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

    # flanking_seqs = {
    #     "type1_up" : Seq("CCCT"), 
    #     "type1_down" : Seq("TTGC"),
    #     }
    

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
    
    def get_assembly_order(self, sequences: Dict[str, SeqRecord]) -> List[str]:
        """Get assembly order from user with autocomplete"""
        sequence_completer = WordCompleter(sequences.keys(), ignore_case=True)
        
        while True:
            try:
                self.console.print("\n[blue]Enter sequences to assemble (space-separated)[/blue]")
                self.console.print("[dim]Use TAB for autocompletion[/dim]")
                self.console.print("[dim]\nRemember that golden gate parts need to be assembled in order[/dim]")
                self.console.print("[blue]\nSeparate components with / to multiplex, use /all[red]X[/red] to select all parts of type x")
                
                user_input = self.session.prompt(
                    "Sequences: ",
                    completer=sequence_completer
                )
                
                selected = user_input.split()
                assemblies = [[]]
                for selection in selected: 
                    print(selection)

                    #if the user has input /allx create new assemblies for each part of type x
                    if selection.startswith("/all"): 
                        self.multiplex = True # flag that this is now a multiplex assembly
                        part_type = selection[4:].strip() #part name should follow /all directly

                        multiplex_part_names = []
                        for part in sequences.keys(): 
                            if part.split("_")[0] == part_type: 
                                multiplex_part_names.append(part)
                        
                        #validate that the part name exists in the library 
                        if not multiplex_part_names: 
                            self.console.print(f"[red]No parts found with type {part_type}. Available types: {set(p.split('_')[0] for p in sequence.keys())}[/red]")

                        print(f"{len(multiplex_part_names)} multiplex assemblies for type {part_type} parts")
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
                        print(f"{len(part_names)} multiplex assemblies for {part_type} parts")   

                        new_assemblies = []
                        for assembly in assemblies: 
                            for part in part_names: 
                                new_assemblies.append(assembly + [part])
                        assemblies = new_assemblies
                                

                    else: #only one option selected, append to all assemblies in place
                        for assembly in assemblies: 
                            assembly.append(selection)

        
        
                # print(assemblies)
                print('number of assemblies')
                # print(len(assemblies))
                self.assemblies = assemblies
                return assemblies

            except Exception as e:
                print(f"Error: {str(e)}")
                
                raise click.Abort()

    
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

    def get_seq_records(self): 
        """"converts self.assembly_sequences to a list of seq records with matching names
        
        Uses self.multiplex to determine if there is more than one plasmid to be assembled

        Raises: 
            ValueError: if no sequences have been seelcted for assembly
        """
        pass 

    def gg_assembly(self): 
        """should maybe be a helper function? or specfic to yeast moclo level 1 assemblies? 
        """
        pass 

    def gg_save_output(self): 
        """Saves the output of the assembly function, creating a new plasmid output file and image for 
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



    