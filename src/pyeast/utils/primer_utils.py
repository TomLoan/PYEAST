# Copyright CSIRO 2025. Thomas Loan 
# See LICENSE for full GpLv2 license. 

# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License or 
# (at your option) any later version. 

# This program is distributed in the hope that it will be useful;, 
# but WITHOUT ANY WARRENTY; without even the implied warranty of 
# MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details. 

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 

# ===========================================================================

# src/pyeast/utils/primer_utils.py
"""
Primer utilities for PYEAST. 
"""

# ===========================================================================






import os 
import io
import math
import logging
from typing import Dict, List, Tuple, Generator, Optional
import pandas as pd
from collections import Counter
import primer3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Align 
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from datetime import date, datetime
from PIL import Image
from prompt_toolkit import PromptSession
from prompt_toolkit.completion import WordCompleter
from prompt_toolkit.shortcuts import confirm

def design_screening_primers(sequence, target_start, target_end, primer_length=20, product_size_range=(500,1000)):
    """
    Design screening primers for a given sequence using Primer3.

    Args:
        sequence (str): The relevant DNA sequence to design primers for.
        target_start (int): The start position of the target region within the given sequence.
        target_end (int): The end position of the target region within the given sequence.
        primer_length (int): The desired length of the primers (default: 20).
        product_size_range (tuple): The desired range of product sizes (default: (500,1000)).

    Returns:
        tuple: A tuple containing (forward_primer, reverse_primer, product_size)
    """
    seq_args = {
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_TARGET': [target_start, target_end - target_start],
    }

    global_args = {
        'PRIMER_OPT_SIZE': primer_length,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range],
        'PRIMER_NUM_RETURN': 1,  # We only want one pair of primers
    }

    try:
        result = primer3.design_primers(seq_args, global_args)

        if 'PRIMER_LEFT_0_SEQUENCE' not in result:
            print("Primer3 failed to design primers. No primers returned.")
            return None, None, None

        forward_primer = result['PRIMER_LEFT_0_SEQUENCE']
        reverse_primer = result['PRIMER_RIGHT_0_SEQUENCE']
        product_size = result['PRIMER_PAIR_0_PRODUCT_SIZE']

        return forward_primer, reverse_primer, product_size

    except primer3.PrimerDesignError as e:
        print(f"Primer3 encountered an error: {str(e)}")
        return None, None, None
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        return None, None, None
    

def design_primers(sequence, target_tm, tolerance=3):
    """
    Design forward and reverse primers for a given sequence with a target melting temperature.
    Args:
        sequence (str): The DNA sequence to design primers for.
        target_tm (float): The target melting temperature for the primers.
        tolerance (float): The allowed deviation from the target melting temperature.
    Returns:
        tuple: A tuple containing the forward and reverse primers (Seq objects).
    """
    def adjust_primer(seq, is_forward):
        if is_forward:
            primer = seq
        else:
            primer = seq.reverse_complement()
            
        while abs(mt.Tm_NN(primer) - target_tm) > tolerance:
            current_tm = mt.Tm_NN(primer)
            if current_tm > target_tm:
                # Remove one base from either end
                primer = primer[:-1]
            elif current_tm < target_tm:
                if is_forward:
                    # Add one more base from original sequence
                    new_length = len(primer) + 1
                    primer = seq[:new_length]
                else:
                    # Add one more base from original sequence for reverse primer
                    new_length = len(primer) + 1
                    primer = seq[-new_length:].reverse_complement()
                    
            # Safety check to prevent infinite loop
            if len(primer) < 16 or len(primer) > len(seq):
                break
                
        return primer

    f_primer = adjust_primer(Seq(sequence[:50]), True)
    r_primer = adjust_primer(Seq(sequence[-50:]), False)
    return f_primer, r_primer

def get_primer_locations(primers: Dict[str, Seq], directory: str) -> Tuple[Dict[str, List], Tuple[Dict[str, List]]]:
    """
    Locate primers in IDT spec sheets within the specified directory.

    This function searches for primers in Excel files containing IDT specifications.
    It returns two dictionaries: one for found primers and one for missing primers.

    Args:
        primers (Dict[str, Seq]): A dictionary of primer names and their sequences.
        directory (str): The directory path containing IDT spec sheets.

    Returns:
        Tuple[Dict[str, List], Tuple[Dict[str, List]]]: A tuple containing two dictionaries:
            - primers_found: Mapping of primer names to lists of matching locations.
            - primers_missing: Mapping of primer names to lists of missing primer information.
    """
    primers_found = {}
    primers_missing = {}
    
    for spec_sheet in os.listdir(directory):
        if spec_sheet.endswith('.xlsx'):
            df = pd.read_excel(os.path.join(directory, spec_sheet), header=0)
            if len(df.columns) == 4 and all(col in df.columns for col in ['Plate or Box ID', 'Position', 'Sequence Name', 'Sequence']):
                df['Sequence'] = df['Sequence'].str.replace(" ", "").str.strip()
                
                for name, sequence in primers.items():
                    str_sequence = str(sequence)
                    matching_rows = df[df['Sequence'].str.upper() == str_sequence.upper()]
                    

                    
                    if not matching_rows.empty:
                        if name not in primers_found:
                            primers_found[name] = []
                        # elif name in primers_found:
                        #     if sequence != primers_found[name]: 
                        #         print("""Warning, two primers with the same name and different sequences detected
                        #               This will cause issues in the batching process, run these constructs in different batches""")
                        primers_found[name].append({
                            'Location': matching_rows["Plate or Box ID"].iloc[0],
                            'Position': matching_rows["Position"].iloc[0],
                            'sequence': sequence
                        }) 

                

    for name, sequence in primers.items():
        if name not in primers_found:
            
            primers_missing[name] = []
            primers_missing[name].append({
                'Location' : 'N/A', 
                'Position' : 'N/A',
                'sequence' : sequence
            }) 
                
            
    #print(primers_found) 
    #print("\n") 
    #print(primers_missing)
    return primers_found, primers_missing


def rationalize_primers(primers_found: Dict[str, List],  primers_missing : Dict[str, List]) -> Dict[str, List]:
    """
    Rationalize primer selection to minimize the number of unique plates or Box used.

    This function chooses primers based on the frequency of plate or box usage and
    includes information about missing primers.

    Args:
        primers_found (Dict[str, List]): A dictionary of found primers and their locations.
        primers_missing (Dict[str, List]): A dictionary of missing primers.

    Returns:
        Dict[str, List]: A dictionary of rationalized primer selections, including missing primers.
    """
    plate_count = {}
    #Count how frequently each plate/box appears across all plates/boxes with matching primers
    for primer in primers_found: 
        #print(primer)
        matches = primers_found[primer]
        #print(matches)
        for match in matches: 
            #print(match)
            plate = match['Location']
            #print(plate)
            if plate in plate_count: 
                plate_count[plate] +=1 
            else: 
                plate_count[plate] = 1

            
    #Select primers from the plates that appear most frequently
    selected_primers = {}
    for primer in primers_found: 
        matches = primers_found[primer] 
        valid_plates = [match['Location'] for match in matches]
        
        if valid_plates: 
            best_plate = max(valid_plates, key=lambda x: plate_count[x])
            best_match = next(match for match in matches if match['Location'] == best_plate)
            selected_primers[primer] = best_match
    
    for primer_name, primer_info_list in primers_missing.items():
        selected_primers[primer_name] = primer_info_list[0]

    return selected_primers

def add_circular_overhangs(primers: Dict[str, Seq], parts: List[SeqRecord], overhang_length: int) -> Dict[str, Seq]:
    """
    Add overhangs to primers for yeast assembly.

    Args:
        primers (Dict[str, Seq]): Dictionary of primers keyed by part name + 'F' or 'R'.
        parts (List[SeqRecord]): List of parts in assembly order.
        overhang_length (int): Length of the overhang to add.

    Returns:
        Dict[str, Seq]: Dictionary of primers with overhangs added.
    """
    oh_oligos = {}
    num_parts = len(parts)

    for i, (primer_name, primer_seq) in enumerate(primers.items()):
        part_index = i // 2
        is_forward = primer_name.endswith('F')

        if is_forward:
            prev_part = parts[part_index - 1] if part_index > 0 else parts[-1]
            overhang = prev_part.seq[-overhang_length:]
            new_name = f"{primer_name}-{prev_part.id}"
        else:
            next_part = parts[(part_index + 1) % num_parts]
            overhang = Seq(str(next_part.seq[:overhang_length])).reverse_complement()
            new_name = f"{primer_name}-{next_part.id}"

        oh_oligos[new_name] = overhang + primer_seq

    return oh_oligos

def design_circular_primers(parts: List[SeqRecord], target_tm=50, overhang_length=25):
    """
    Design TAR cloning primers for a set of DNA parts.

    Args:
        parts (List[SeqRecord]): List of parts in assembly order.
        target_tm (float): The target melting temperature for the primers.
        overhang_length (int): Length of the overhang to add for assembly.

    Returns:
        dict: Dictionary of primers with overhangs for TAR cloning.
    """
    primers = {}
    for i, part in enumerate(parts):
        part_name = f"{part.id}_{i}"  # Use a unique identifier for each part
        
        f_primer, r_primer = design_primers(str(part.seq), target_tm)
        primers[f"{part_name}F"] = f_primer
        primers[f"{part_name}R"] = r_primer

    return add_circular_overhangs(primers, parts, overhang_length)

def add_linear_overhangs(primers: Dict[str, Seq], parts: List[SeqRecord], int_site: Tuple[SeqRecord, SeqRecord], overhang_length: int) -> Dict[str, Seq]:
    """
    Add overhangs to primers for linear integration assembly.

    Args:
        primers (Dict[str, Seq]): Dictionary of primers keyed by part name + 'F' or 'R'.
        parts (List[SeqRecord]): List of parts in assembly order (not including integration sites).
        int_site (Tuple[SeqRecord, SeqRecord]): Integration site sequences (up, down).
        overhang_length (int): Length of the overhang to add.

    Returns:
        Dict[str, Seq]: Dictionary of primers with overhangs added where appropriate.
    """
    oh_oligos = {}
    int_up, int_down = int_site
    all_parts = [int_up] + parts + [int_down]
    num_parts = len(all_parts)

    for i, (primer_name, primer_seq) in enumerate(primers.items()):
        part_index = i // 2
        is_forward = primer_name.endswith('F')

        if part_index == 0 and is_forward:
            # First forward primer (int_up_F) - no overhang
            oh_oligos[primer_name] = primer_seq
        elif part_index == num_parts - 1 and not is_forward:
            # Last reverse primer (int_down_R) - no overhang
            oh_oligos[primer_name] = primer_seq
        elif is_forward:
            # All other forward primers
            prev_part = all_parts[part_index - 1]
            overhang = prev_part.seq[-overhang_length:]
            new_name = f"{primer_name}-{prev_part.id}"
            oh_oligos[new_name] = overhang + primer_seq
        else:
            # All other reverse primers
            next_part = all_parts[part_index + 1]
            overhang = Seq(str(next_part.seq[:overhang_length])).reverse_complement()
            new_name = f"{primer_name}-{next_part.id}"
            oh_oligos[new_name] = overhang + primer_seq

    return oh_oligos

def design_linear_primers(components: List[SeqRecord], int_site: Tuple[SeqRecord, SeqRecord], 
                               target_tm: float, homology_length: int) -> Dict[str, Seq]:
    """
    Design primers for integrating components at a specific site.

    Args:
        components (List[SeqRecord]): List of components to be integrated.
        int_site (Tuple[SeqRecord, SeqRecord]): The integration site sequences (up, down).
        target_tm (float): The target melting temperature for the primers.
        homology_length (int): Length of homology for overlaps.

    Returns:
        Dict[str, Seq]: Dictionary of designed primers with overhangs where appropriate.
    """
    primers = {}
    int_up, int_down = int_site
    all_parts = [int_up] + components + [int_down]

    # Design initial primers without overhangs
    for part in all_parts:
        f_primer, r_primer = design_primers(str(part.seq), target_tm)
        primers[f"{part.id}_F"] = f_primer
        primers[f"{part.id}_R"] = r_primer

    # Add overhangs
    primers_with_overhangs = add_linear_overhangs(primers, components, int_site, homology_length)

    return primers_with_overhangs
