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

# src/pyeast/utils/sequence_utils.py
"""
Sequence utilities for PYEAST. 
"""

# ===========================================================================




import os 
import io
import math
import logging
from typing import Dict, List, Tuple, Generator, Optional
import pandas as pd
from collections import Counter
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

def load_sequences(directory: str) -> Dict[str, SeqRecord]:
    """
    Load all sequences from FASTA files in the given directory.
    Uses the FASTA header as the sequence name.

    Args:
        directory (str): Path to the directory containing FASTA files.

    Returns:
        Dict[str, Sequence]: Dictionary of loaded sequences.

    Raises:
        FileNotFoundError: If the directory is not found.
    """
    sequences = {}
    try:
        for filename in os.listdir(directory):
            if filename.endswith(('.fasta', '.fa', 'fsa')):
                file_path = os.path.join(directory, filename)
                for record in SeqIO.parse(file_path, "fasta"): #Change to fasta-person for biopython >= 1.85
                    name = record.id.split()[0]  # Take first word of the header
                    sequences[name] = record
                    #logging.info(f"Loaded sequence: {name} from file: {filename}")
    except FileNotFoundError:
        logging.error(f"Directory not found: {directory}")
        raise
    except Exception as e:
        logging.error(f"Error loading sequences: {str(e)}")
        raise

    return sequences

def get_templates(parts: List[SeqRecord], directory: str) -> Dict[str, List[str]]:
    """
    Find template matches for given parts in the specified directory.

    This function searches for templates in FASTA or GenBank files that contain
    the sequences of the provided parts.

    Args:
        parts (List[SeqRecord]): A list of SeqRecord objects representing the parts.
        directory (str): The directory path containing template files.

    Returns:
        Dict[str, List[str]]: A dictionary mapping part IDs to lists of matching template IDs.
    """
    templates = {}
    templates_used = {}

    for template_file in os.listdir(directory):
        file_path = os.path.join(directory, template_file)
        if file_path.endswith(('.fasta', '.fsa', '.fa', '.gb', '.gbk')):
            format = "fasta" if file_path.endswith(('.fasta', '.fa','.fsa')) else "genbank"
            for record in SeqIO.parse(file_path, format):
                #print(record.name)
                #print(record.annotations)
                templates[record.name] = record.seq

    for part in parts:
        templates_used[part.name] = []
        for template_id, template_seq in templates.items():
            if part.seq.upper() in template_seq.upper() or part.seq.upper() in template_seq.upper().reverse_complement():
                templates_used[part.name].append(template_id)
        
        if not templates_used[part.name]:
            templates_used[part.name] = ["Not found"]
    #print(templates_used)
    return templates_used

def rationalize_templates(template_dict: Dict[str, List[str]]) -> Dict[str, str]:
    """
    Rationalize template selection to minimize the number of unique templates.

    This function chooses templates based on the following priority:
    1. Templates in the preferred list
    2. Highest frequency of use across all parts
    3. Shorter names

    Args:
        template_dict (Dict[str, List[str]]): A dictionary mapping part names to lists of potential templates.

    Returns:
        Dict[str, str]: A dictionary mapping part names to their choseninput template.
    """
    # Define a list of preferred templates
    preferred_templates = ['pUC19', 'pYES2', 'pESC-TRP', ]  # Add more as needed

    # Count the global frequency of each template
    all_templates = [template for templates in template_dict.values() for template in templates if template != "Not found"]
    template_frequency = Counter(all_templates)

    def choose_best_template(templates):
        if not templates or templates[0] == "Not found":
            return "Not found"
        
        # Filter out "Not found" entries
        valid_templates = [t for t in templates if t != "Not found"]
        
        # First, check if any templates are in the preferred list
        preferred_available = [t for t in valid_templates if t in preferred_templates]
        if preferred_available:
            return min(preferred_available, key=lambda x: (preferred_templates.index(x), -template_frequency[x], len(x)))
        
        # If no preferred templates, choose based on frequency and then name length
        return max(valid_templates, key=lambda x: (template_frequency[x], -len(x)))

    # Select the best template for each part
    rationalized_templates = {part: choose_best_template(templates) 
                              for part, templates in template_dict.items()}

    return rationalized_templates

def parse_gb_file(file_path: str) -> Tuple[SeqRecord, List[SeqRecord]]:
    """
    Parse a GenBank file and extract the full plasmid sequence and its parts.

    Args:
    file_path (str): Path to the GenBank file

    Returns:
    Tuple[SeqRecord, List[SeqRecord]]: Full plasmid sequence and a list of parts
    """
    plasmid = SeqIO.read(file_path, "genbank")
    parts = []
    for feature in plasmid.features:
        if feature.type == "misc_feature":
            part_name = feature.qualifiers.get("label", [f"Part_{len(parts)}"])[0]
            part_seq = feature.extract(plasmid.seq)
            parts.append(SeqRecord(part_seq, id=part_name, name=part_name, description=""))
    return plasmid, parts

def write_circular_instructions(rationalized_primers: Dict[str, Dict],
                              rationalized_templates: Dict[str, str],
                              assembly_sequences: List[SeqRecord],
                              homology_length: int) -> List[List[str]]:
    """Generate assembly instructions based on rationalized primers and templates.

    Args:
    rationalized_primers (Dict[str, Dict]): Information about rationalized primers
    rationalized_templates (Dict[str, str]): Mapping of parts to their templates
    assembly_sequences (List[SeqRecord]): Ordered list of parts to be assembled
    homology_length (int): Length of the homology region in primers

    Returns:
    List[List[str]]: List of instruction rows for assembly
    """
    instructions = []
    
    for part in assembly_sequences:
        part_name = part.id
        part_seq = part.seq
        
        # Find matching primers
        f_primer = find_matching_primer(rationalized_primers, part_seq, True, homology_length)
        r_primer = find_matching_primer(rationalized_primers, part_seq, False, homology_length)
        
        # Get template
        template = rationalized_templates.get(part_name, "Not found")
        
        # Calculate amplicon length
        amplicon_length = len(part_seq) + (2 * homology_length)
        
        # Create instruction row
        if f_primer and r_primer:
            instruction = [
                part_name,
                f_primer['Location'], f_primer['name'], f_primer['Position'],
                r_primer['Location'], r_primer['name'], r_primer['Position'],
                template,
                amplicon_length
            ]
            instructions.append(instruction)
    
    return instructions

def write_linear_instructions(rationalized_primers: Dict[str, Dict],
                            rationalized_templates: Dict[str, str],
                            int_site_up: SeqRecord,
                            middle_sequences: List[SeqRecord],
                            int_site_down: SeqRecord,
                            homology_length: int) -> List[List[str]]:
    """Generate assembly instructions for linear integration assembly.

    Args:
    rationalized_primers (Dict[str, Dict]): Information about rationalized primers
    rationalized_templates (Dict[str, str]): Mapping of parts to their templates
    int_site_up (SeqRecord): Upstream integration site sequence
    middle_sequences (List[SeqRecord]): List of sequences to be integrated
    int_site_down (SeqRecord): Downstream integration site sequence
    homology_length (int): Length of the homology region in primers

    Returns:
    List[List[str]]: List of instruction rows for assembly
    """
    instructions = []
    all_parts = [int_site_up] + middle_sequences + [int_site_down]
    
    for part in all_parts:
        part_name = part.id
        part_seq = part.seq
        
        # Find matching primers
        if part_seq == int_site_up.seq:
            f_primer = find_matching_primer(rationalized_primers, part_seq, True, 0)
            r_primer = find_matching_primer(rationalized_primers, part_seq, False, homology_length)
        elif part_seq == int_site_down.seq: 
            f_primer = find_matching_primer(rationalized_primers, part_seq, True, homology_length)
            r_primer = find_matching_primer(rationalized_primers, part_seq, False, 0)
        else:
            f_primer = find_matching_primer(rationalized_primers, part_seq, True, homology_length)
            r_primer = find_matching_primer(rationalized_primers, part_seq, False, homology_length)
        
        # Get template
        template = rationalized_templates.get(part_name, "Not found")
        
        # Calculate amplicon length (accounting for no overhangs on ends)
        if part_seq == int_site_up.seq:
            amplicon_length = len(part_seq) + homology_length  # Only downstream overhang
        elif part_seq == int_site_down.seq:
            amplicon_length = len(part_seq) + homology_length  # Only upstream overhang
        else:
            amplicon_length = len(part_seq) + (2 * homology_length)  # Both overhangs
        
        # Create instruction row
        if f_primer and r_primer:
            instruction = [
                part_name,
                f_primer['Location'], f_primer['name'], f_primer['Position'],
                r_primer['Location'], r_primer['name'], r_primer['Position'],
                template,
                amplicon_length
            ]
            instructions.append(instruction)
    
    return instructions

def find_matching_primer(primers: Dict[str, Dict], part_seq: Seq, is_forward: bool, homology_length: int) -> Optional[Dict]:
    """
    Find a matching primer for a given part sequence.

    Args:
    primers (Dict[str, Dict]): Dictionary of primers and their information
    part_seq (Seq): The sequence of the part to match against
    is_forward (bool): True if searching for a forward primer, False for reverse
    homology_length (int): Length of the homology region to exclude from matching

    Returns:
    Optional[Dict]: Matching primer information or None if no match found
    """
    part_str = str(part_seq)
    for name, info in primers.items():
        primer_seq = info['sequence']
        if not isinstance(primer_seq, Seq):
            primer_seq = Seq(primer_seq)
        
        # Remove homology region from primer
        primer_without_homology = primer_seq[homology_length:]
        check_length = len(primer_without_homology)

        if is_forward:
            if str(primer_without_homology) == part_str[:check_length] and check_length >14:
                return {**info, 'name': name}
        else:
            rev_comp = primer_without_homology.reverse_complement()
            if str(rev_comp) == part_str[-check_length:] and check_length >14:
                return {**info, 'name': name}
    
    return None

def assemble_parts_circular(parts: List[SeqRecord], primers: Dict[str, Seq], homology_length: int) -> SeqRecord:
    """
    Assemble DNA parts using TAR cloning simulation, preserving part information and correctly positioning primers.

    Args:
        parts (List[SeqRecord]): List of parts to be assembled.
        primers (Dict[str, Seq]): Dictionary of primers with overhangs.
        homology_length (int): Length of homology for recombination.

    Returns:
        SeqRecord: Assembled plasmid as a SeqRecord object with features.
    """
    # Assemble the plasmid sequence
    assembled_sequence = "".join(str(part.seq) for part in parts)
    total_length = len(assembled_sequence)
    circular_sequence = assembled_sequence + assembled_sequence[:50]
    
    features = []
    current_position = 0

    # Add part features
    for i, part in enumerate(parts):
        part_length = len(part.seq)
        part_feature = SeqFeature(
            FeatureLocation(current_position, current_position + part_length),
            type="misc_feature",
            qualifiers={"label": f"{part.id}_{i}"}  # Add index to make labels unique
        )
        features.append(part_feature)
        current_position += part_length

    # Add primer features
    for primer_name, primer_seq in primers.items():
        primer_str = str(primer_seq)
        rc_primer_str = str(Seq(primer_str).reverse_complement())
        
        # Find all occurrences of the primer in the circular sequence
        forward_matches = find_all_occurrences(circular_sequence, primer_str)
        reverse_matches = find_all_occurrences(circular_sequence, rc_primer_str)
        
        # Add features for all forward matches
        for match in forward_matches:
            start = match % total_length
            end = (start + len(primer_str)) % total_length
            if start < end:
                primer_feature = SeqFeature(
                    FeatureLocation(start, end, strand=1),
                    type="primer_bind",
                    qualifiers={"label": f"{primer_name}_forward"}
                )
            else:
                primer_feature = SeqFeature(
                    CompoundLocation([
                        FeatureLocation(start, total_length, strand=1),
                        FeatureLocation(0, end, strand=1)
                    ]),
                    type="primer_bind",
                    qualifiers={"label": f"{primer_name}_forward"}
                )
            features.append(primer_feature)
        
        # Add features for all reverse matches
        for match in reverse_matches:
            start = match % total_length
            end = (start + len(primer_str)) % total_length
            if start < end:
                primer_feature = SeqFeature(
                    FeatureLocation(start, end, strand=-1),
                    type="primer_bind",
                    qualifiers={"label": f"{primer_name}_reverse"}
                )
            else:
                primer_feature = SeqFeature(
                    CompoundLocation([
                        FeatureLocation(0, end, strand=-1),
                        FeatureLocation(start, total_length, strand=-1)
                        
                    ]),
                    type="primer_bind",
                    qualifiers={"label": f"{primer_name}_reverse"}
                )
            features.append(primer_feature)
        
        if not forward_matches and not reverse_matches:
            print(f"Warning: Primer {primer_name} not found in the assembled sequence")

    # Check for similar junctions
    junction_length = 100  # 50 bp on either side of the junction
    junctions = []
    
    current_position = 0
    for i, part in enumerate(parts):
        part_length = len(part.seq)
        junction_start = (current_position + part_length - junction_length // 2) % total_length
        junction_end = (junction_start + junction_length) % total_length
        
        if junction_start < junction_end:
            junction_seq = circular_sequence[junction_start:junction_end]
        else:
            junction_seq = circular_sequence[junction_start:]
        
        next_part = parts[(i + 1) % len(parts)]
        junctions.append((junction_seq, f"{part.id}_{i}", f"{next_part.id}_{(i+1)%len(parts)}"))
        
        current_position += part_length
    #print('checking junctions')
    #print(junctions)
    # Perform pairwise alignments and check for similarities
    # print(len(junctions))
    for i in range(len(junctions)):
        for j in range(i + 1, len(junctions)):
            seq1, part1, next_part1 = junctions[i]
            seq2, part2, next_part2 = junctions[j]
            
            alignments = Align.PairwiseAligner().align(seq1, seq2)
            best_alignment = alignments[0]
            max_possible_score = min(len(seq1), len(seq2))
            similarity = best_alignment.score / max_possible_score
            #print(best_alignment)
            if similarity > 0.8:  # You can adjust this threshold
                print(f"Warning: High similarity ({similarity:.2f}) between junctions:")
                print(f"  - {part1} and {next_part1}")
                print(f"  - {part2} and {next_part2}")
                print(f"Alignment:\n{best_alignment}")

    # Create a SeqRecord for the assembled plasmid
    assembled_plasmid = SeqRecord(
        Seq(assembled_sequence),
        id="Assembled_Plasmid",
        name="TAR_Cloned_Plasmid",
        description="Plasmid assembled by TAR cloning simulation",
        features=features
    )

    # Add annotations
    assembled_plasmid.annotations["molecule_type"] = "DNA"
    assembled_plasmid.annotations["topology"] = "circular"
    assembled_plasmid.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()

    return assembled_plasmid

def find_all_occurrences(sequence: str, substring: str) -> Generator[int, None, None]:
    """Find all occurrences of a substring in a sequence."""
    start = 0
    while True:
        start = sequence.find(substring, start)
        if start == -1:  # Substring not found
            return
        yield start
        start += 1  # Move to next possible position


def assemble_parts_linear(parts: List[SeqRecord], primers: Dict[str, Seq]) -> SeqRecord:
    """Assemble parts for integration, add features for components and primers.
    
    Args:
        parts: List of all parts in assembly order (including integration sites)
        primers: Dictionary of primer names and sequences
    
    Returns:
        SeqRecord: Assembled sequence with features and primer annotations
    """
    # Assemble parts
    assembled_sequence = SeqRecord(
        Seq(''.join(str(part.seq) for part in parts)),
        id="Assembled_Integration_Construct",
        name="Integration_Construct",
        description="Assembled sequence for genomic integration"
    )

    # Add features for each part
    current_position = 0
    for i, part in enumerate(parts):
        part_length = len(part.seq)
        if i == 0:
            feature_type = "misc_feature"
            qualifier = {"label": f"{part.id} (upstream)"}
        elif i == len(parts) - 1:
            feature_type = "misc_feature"
            qualifier = {"label": f"{part.id} (downstream)"}
        else:
            feature_type = "misc_feature"
            qualifier = {"label": part.id}
        
        feature = SeqFeature(
            FeatureLocation(current_position, current_position + part_length),
            type=feature_type,
            qualifiers=qualifier
        )
        assembled_sequence.features.append(feature)
        current_position += part_length

    # Add primer annotations
    assembled_seq_str = str(assembled_sequence.seq)
    for primer_name, primer_seq in primers.items():
        primer_seq_str = str(primer_seq)
        forward_pos = assembled_seq_str.find(primer_seq_str)
        reverse_complement = str(primer_seq.reverse_complement())
        reverse_pos = assembled_seq_str.find(reverse_complement)
        
        if forward_pos != -1:
            primer_feature = SeqFeature(
                FeatureLocation(forward_pos, forward_pos + len(primer_seq), strand=1),
                type="primer_bind",
                qualifiers={"label": primer_name}
            )
            assembled_sequence.features.append(primer_feature)
        elif reverse_pos != -1:
            primer_feature = SeqFeature(
                FeatureLocation(reverse_pos, reverse_pos + len(primer_seq), strand=-1),
                type="primer_bind",
                qualifiers={"label": primer_name}
            )
            assembled_sequence.features.append(primer_feature)

    # Set annotations
    assembled_sequence.annotations["molecule_type"] = "DNA"
    assembled_sequence.annotations["topology"] = "linear"
    assembled_sequence.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
    
    return assembled_sequence

