# Golden gate of level 1 plasmids using the MoClo YTK design 



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

# src/pyeast/core/gg_lvl1.py
"""
Level 1 golden gate assemblies using the yeast Moclo standard from Lee et al 2015
"""

# ===========================================================================



from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq
from rich.console import Console
from rich.table import Table
from rich.progress import Progress
import click
from prompt_toolkit import PromptSession 
from prompt_toolkit.completion import WordCompleter
from prompt_toolkit.shortcuts import confirm
import io
from PIL import Image
import dnacauldron as dc

from ..utils.sequence_utils import (
    load_sequences,
    get_templates,  
) 

from ..utils.visualisation import visualise_genbank, save_figure


class gg_lvl1Designer: 
    """A class for designing golden gate level 1 assemblies. 
    
    This class handles the design assemblies according to the 
    parameters in the yeast MoClo took kit laid out in Lee et al
    2015
    """
    #Flanking sequences specific to different part types, should specifiy part type as a string
    #'3a', '2', etc. and sequences as tuple of strings e.g. '1': ('5'seq', '3'seq')
    FLANKING_SEQs = {
        "1": ('CCCT','TTGC'),
        "2": ('AACG','ATAC'),
        "3": ('TATG', 'TAGG'),
        "3a": ('TATG', 'AAGA'),
        "3a'": ('TATG','TGCT'),
        "3b": ('TTCT', 'TAGG'),
        "3b'":('ACGA', 'TAGG'), 
        "4": ('ATCC', 'CGAC'),
        "4a": ('ATCC', 'ACCG'),
        "4b": ('TGGC', 'CGAC'), 
        "5": ('GCTG', 'ATGT'),
        "6": ('TACA', 'CTCA'),
        "7": ('GAGT', 'GGCT'), 
        "8": ('CCGA', 'GGGA'),
        "8a": ('CCGA', 'GTTA'), 
        "8b": ('CAAT', 'GGGA'), 
        }
    

    def __init__(self): 
        """Inintialise a new gg_lvl1 designer. 
        
        Args:
            homology_length: Length of homology regions to be added to primers (default: 25) 
            annealing_temp: Target annealing temperature for primer design (default: 50)
            """
        
        self.console = Console()
        self.session = PromptSession()
        self.multiplex = False
        self.available_sequences = {}     # All loaded sequences
        self.assembly_sequences = []      # Sequences in assembly order may be multiple assemblies for multiplex work (list of lists)
        self.assemblies_names = None
        self.template_dict = {}           # Template information
        self.final_assembled_sequence = None    # Final assembled sequence


    def load_and_get_sequences(self, directory: Path) -> None: 
        """"Loads seqeunces from a directory and store them for assembly. 
        
        Args: 
            directory: Path to a directory containing fasta files containing DNA components 
            as fasta files
            """
        self.available_sequences = load_sequences(directory)
        return self.available_sequences
    
    def print_sequence_grid(self, sequences: dict, title: str = "Available Sequences"):
        """Display available sequences in a formatted table"""
        table = Table(title=title)
        table.add_column("Name", style="cyan")
        table.add_column("Length", justify="right", style="green")
        table.add_column("Description", style="white")
        
        for name, seq in sequences.items():
            table.add_row(
                name,
                f"{len(seq)} bp",
                seq.description[len(name):150+len(name)] + "..." if len(seq.description)-len(name) > 149 else seq.description[len(name):]
            )
        
        self.console.print(table) 

    def get_assembly_order(self, sequences: Dict[str, SeqRecord]) -> List[list[str]]:
        """Get assembly order from user with autocomplete and mutliplex support 

        prompts the user to enter sequences for assembly with autocomplete. 
        Supports multiplex assemblies using '/' syntax to select multiple components 
        in a given position, or /allx to select all parts of a give type x. 
        
        Args: 
            sequences: Dictionary mapping sequence names to SeqRecord objects 

        Returns: 
            List[list[str]]: list of assemblies where each assembly is a list of sequence
                            names. For single assemblies returns [[seq1,seq2 ...]]
                            For multiplex assemblies returns [[seq1, seq2a,..], [seq1, seq2b...]]
        Raises: 
            click.Abort: if user cancels the operation
            ValueError If no valid sequences are available

        Examples:
        Single assembly: "promoter_1 gene_1 terminator_1" -> [["promoter_1", "gene_1", "terminator_1"]]
        Multiplex: "promoter_1 gene_1/gene_2 terminator_1" -> [["promoter_1", "gene_1", "terminator_1"], 
                                                               ["promoter_1", "gene_2", "terminator_1"]]
        All parts: "promoter_1 /all3 terminator_1" -> Multiple assemblies with type 3 parts
        """
        if not sequences: 
            raise ValueError("No sequences available for assembly")
        
        sequence_completer = WordCompleter(sequences.keys(), ignore_case=True)
        
        while True:
            try:
                self.console.print("\n[blue]Enter sequences to assemble (space-separated)[/blue]")
                self.console.print("[dim]Use TAB for autocompletion[/dim]")
                self.console.print("[dim]Remember that golden gate parts need to be assembled in order[/dim]")
                self.console.print("[blue]Separate components with / to multiplex, use /all[red]X[/red] to select all parts of type [red]x[/red]")
                
                user_input = self.session.prompt(
                    "Sequences: ",
                    completer=sequence_completer
                )
                
                selected = user_input.split()
                if not selected: 
                    self.console.print("[yellow]No sequences selected[/yellow]")
                    continue

                #validate selections befor processing
                invalid_selections = []
                for selection in selected: 
                    if selection.startswith("/all"): 
                        part_type = selection[4:].strip()
                        matching_parts = [name for name in sequences.keys() if name.split("_")[0] == part_type]
                        if not matching_parts: 
                            invalid_selections.append(f"{selection} (no parts found with type '{part_type}')")

                    elif "/" in selection: 
                        part_names = selection.split()
                        for part_name in part_names: 
                            if part_name not in sequences: 
                                invalid_selections.append(f"{part_name} (from multplex selection '{selection})")
                    
                    else: 
                        if selection not in sequences: 
                            invalid_selections.append(selection)
                
                #send user back to select if there are invalid entries
                if invalid_selections: 
                    self.console.print(f"[red]Invalid selection(s): {','.join(invalid_selections)}[/red]")
                    continue
                
                # Process valid selections
                assemblies = [[]]
                for selection in selected: 
                    
                    #if the user has input /allx create new assemblies for each part of type x
                    if selection.startswith("/all"): 
                        self.multiplex = True # flag that this is now a multiplex assembly
                        part_type = selection[4:].strip() #part name should follow /all directly

                        multiplex_part_names = []
                        for part in sequences.keys(): 
                            if part.split("_")[0] == part_type: 
                                multiplex_part_names.append(part)
                        
                        print(f"{len(multiplex_part_names)} type {part_type} parts for multiplexing")
                        
                        #validate that the part name exists in the library 
                        if not multiplex_part_names: 
                            self.console.print(f"[red]No parts found with type {part_type}. Available types: {set(p.split('_')[0] for p in sequences.keys())}[/red]")
                      
                        #create new assemblies each with one of the multiplex parts appended
                        new_assemblies = []
                        for assembly in assemblies: 
                            for part in multiplex_part_names: 
                                new_assemblies.append(assembly + [part])
                        assemblies = new_assemblies 



                    elif "/" in selection: 
                        self.multiplex = True
                        part_type = selection.split("_")[0]
                          
                        part_names = selection.split("/")
                        print(f"{len(part_names)} type {part_type} parts for multiplexing")   

                        new_assemblies = []
                        for assembly in assemblies: 
                            for part in part_names: 
                                new_assemblies.append(assembly + [part])
                        assemblies = new_assemblies
                        # print(f"{len(assemblies)} multiplex assemblies with type {part_type} parts")
                                

                    else: #only one option selected, append to all assemblies in place
                        for assembly in assemblies: 
                            assembly.append(selection)

        
                #provide some feed back on how many constructs have been desinged in the combinatorial selection
                if self.multiplex:
                    print(f"{len(assemblies)} constructs for assembly")
                else: 
                    print("one construct for assembly")
                
                self.assemblies_names = assemblies
                return assemblies

            except Exception as e:
                print(f"Error: {str(e)}")
                raise click.Abort()
            
            except KeyboardInterrupt:
                if confirm("\nDo you want to exit?"):
                    raise click.Abort()
                continue
    


    def get_seq_records(self): 
        """"converts self.assembly_sequences to a list of SeqRecord lists with matching names
        
        Converts the list of assembly names lists into SeqRecord objects from self.available_seqeunces
        Returns a list of lists of SeqRecords for both single and multiplex assemblies

        Returns: 
            List[List[SeqRecord]]: List of assemblies, each of which are lists of Seqrecord objects
                    Single assembly: [[Seq1, Seq2, ..]]
                    Multiplex assembly : [[Seq1, Seq2, ..][Seq1, Seq2a, ...]...]
        
        Raises: 
            ValueError: if no sequences have been seelcted for assembly
            KeyError: if any sequence name in assemblies is not found in available sequences
        """
        if not hasattr(self, 'assemblies') or not self.assemblies_names:
            raise ValueError("No assemblies selected. Please run get_assembly_order first.")
        
        if not self.available_sequences:
            raise ValueError("No sequences available. Please load sequences first.")
        
        #get the SeqRecords from the names and populate a list
        for assembly in self.assemblies_names: 
            seq_records = [] 
            for seq_name in assembly: 
                if seq_name not in self.available_sequences: 
                    raise KeyError(f"Sequence {seq_name} not found in availble sequences")
                seq_records.append(self.available_sequences[seq_name])
            self.assembly_sequences.append(seq_records)
        
        return self.assembly_sequences 
    
    def get_get_plasmid_names(self, plasmid_folder): 
        """this will map the sequences in self.assembly_sequences onto plasmids containing those sequences
        in the speficied directory, sets self.plasmid_names = list[str]
        
        Args: 
        
        Returns: 
        
        Raises: 
        
        """
        pass 

    def gg_assembly(self): 
        """Passes the set of plasmids in self.plasmid_names to dc.Type2RestrictionAssembly with 
        expected_constructs = len(self.assemblies_names), if there is warnings or errors it will 
        provide feedback on what went wrong and return the user to the sequence selection stage. 
        Sets self.final_final = AssemblySimulation object. Note this method doesn't know what the
        intention of the user is - that's restricted to the selection stage, but will check that 
        the outcome is roughly as expected. With this method the order of entry won't matter. 
        """
        pass 

    def gg_save_output(self): 
        """Saves the output of the assembly function, creating a new plasmid output file
        single assemblies, or a folder full of .gb files (no images) for multiplex assemblies
        """
        pass
    
    def gg_instructions(self): 
        """Generates human readable and [optionally] machine specific instructions 
        this will just boil down to a table of what plasmids to include and where they are
        I will save this in the appropriate folder - but I'd like to make it possible to 
        run a gg_lvl1 with an option to start here from an existing folder full of plasmids 
        to be assembled - e.g. if you realize you want to use a different machine.
        this might actaully be clearer as a separate command simlar to batch.
        """
        pass


    def find_templates(self, template_folder: Path) -> None:
        """Find template matches for each assembly component.
        
        Uses template functions from sequence_utils to identify potential
        templates for each sequence and rationalize the selections.
        
        Args:
            template_folder: Path to folder containing template files
            
        Raises:
            ValueError: If no sequences have been selected for assembly
        """
        if not self.assembly_sequences:
            raise ValueError("No sequences selected. Please select sequences first.")
        
        if not self.multiplex: 
            self.template_dict = get_templates(self.assembly_sequences, str(template_folder))



    