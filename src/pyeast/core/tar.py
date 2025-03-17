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

# src/pyeast/core/tar.py
"""
Transformation Assisted Recombinaiton (TAR) toolkit for plasmid assembly in S. cerevisiae.
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
from prompt_toolkit import PromptSession 
from prompt_toolkit.completion import WordCompleter
from prompt_toolkit.shortcuts import confirm
import io
from PIL import Image

from ..utils.sequence_utils import (
    load_sequences,
    get_templates, 
    rationalize_templates, 
    write_circular_instructions, 
    assemble_parts_circular
) 
                                    
from ..utils.primer_utils import (
    design_circular_primers, 
    get_primer_locations, 
    rationalize_primers, 

)

from ..utils.visualisation import visualise_genbank, save_figure


class TARDesigner: 
    """A class for designing TAR experiments . 
    
    This class handles the design of primers and assembly 
    strategies for Transformation Assisted Recombinaiton (TAR)
    in Sacharomyces cerevisiae. 
    """

    def __init__(self, homology_length: int = 25, annealing_temp: float = 50): 
        """Inintialise a new TARDesigner. 
        
        Args:
            homology_length: Length of homology regions to be added to primers (default: 25) 
            annealing_temp: Target annealing temperature for primer design (default: 50)
            """
        
        self.homology_length = homology_length 
        self.annealing_temp = annealing_temp 
        self.console = Console()
        self.session = PromptSession()

        self.available_sequences = {}     # All loaded sequences
        self.assembly_sequences = []      # Sequences in assembly order
        self.primers = {}                 # Designed primers
        self.primers_found = {}           # Primers found in plates
        self.missing_primers = {}         # Primers that need ordering
        self.template_dict = {}           # Template information
        self.rationalized_templates = {}  # Final template selections
        self.rationalized_primers = {}    # Final primer selections
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
                
                user_input = self.session.prompt(
                    "Sequences: ",
                    completer=sequence_completer
                )
                
                selected = user_input.split()
                if not selected:
                    self.console.prinpt("[yellow]No sequences selected[/yellow]")
                    continue

                # Validate selections
                invalid = [name for name in selected if name not in sequences]
                if invalid:
                    self.console.print(f"[red]Invalid sequence(s): {', '.join(invalid)}[/red]")
                    continue

                # Show selection and confirm
                self.console.print("\n[green]Selected sequences:[/green]")
                for i, name in enumerate(selected, 1):
                    self.console.print(f"{i}. {name}")

                if confirm("\nProceed with these sequences?"):
                    return selected
                
            except KeyboardInterrupt:
                if confirm("\nDo you want to exit?"):
                    return None
                continue

    def display_instructions(self, instructions: List[List[str]]):
        """Display assembly instructions in a formatted table"""
        table = Table(title="Assembly Instructions")
        table.add_column("Part Name", style="bold cyan")
        table.add_column("Fwd Primer", style="green")
        table.add_column("Fwd Plate", style = "bold cyan")
        table.add_column("Well", style="bold blue")
        table.add_column("Rev Primer", style="green")
        table.add_column("Rev Plate", style = "bold cyan")
        table.add_column("Well", style="bold blue")
        table.add_column("Template", style="magenta")
        table.add_column("Size", justify="right")
        
        for line in instructions:
            table.add_row(
                str(line[0]),  # Part name
                str(line[2]),   # Fwd primer name
                str(line[1] if line[1] != "N/A" else "[bold red]Not found[/bold red]"), # Fwd primer plate
                str(line[3] if line[3] != "N/A" else "[bold red]Not found[/bold red]"),  # Fwdc well
                str(line[5]),  # Rev primer name
                str(line[4] if line[4] != "N/A" else "[bold red]Not found[/bold red]"),   # Rev primer Plate
                str(line[6] if line[6] != "N/A" else "[bold red]Not found[/bold red]"),  # Rev well
                str(line[7] if line[7] != "Not found" else "[bold red]Not found[/bold red]"), #Template
                str(line[8])  # Amplicon Size
            )
        
        self.console.print(table)
    
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
                seq.description[:150] + "..." if len(seq.description) > 149 else seq.description
            )
        
        self.console.print(table)

    def set_assembly_order(self, selected_names: List[str])-> None: 
        """Set the assmebly order from user selected seqeunce names"""

        self.assembly_sequences = [
            self.available_sequences[name] for name in selected_names 
        ]

    def design_tar_primers(self) -> Dict[str, Seq]: 
        """Design TAR primers with the necesary overhangs
        
        uses the design_circular_primers function from primer_utils to design primers 
        for the selected assembly sequences. 
        
        Returns: 
            Dictionary mapping primer names to their sequence
            
        Raises: 
            ValueError: If no sequences have been selected for assembly. 
        """

        if not self.assembly_sequences: 
            raise ValueError("No sequences selected for assembly. Please select sequences first")
        
        self.primers = design_circular_primers(
            self.assembly_sequences, 
            target_tm = self.annealing_temp, 
            overhang_length = self.homology_length
        )

        return self.primers

    def check_primer_locations(self, primer_folder: Path) -> None: 
        """Check for exisiting primers stored in plates
        
        Uses primer locations functions from primer_utils to find exisiting primers in 
        IDT plate maps  and identify which primers need to be ordered. 
        
        Args: 
            primer_folder: Path to folder containing primer Excel files

        Raises: 
            ValueError: if primers haven't been desinged yet
        """ 
        if not self.primers: 
            raise ValueError("No primers to look for, please design primers first")
        
        self.primers_found, self.missing_primers = get_primer_locations(
            self.primers, 
            str(primer_folder)
        )
        

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

        self.template_dict = get_templates(self.assembly_sequences, str(template_folder))


    def rationalize_selections(self) -> Tuple[Dict[str, Dict], Dict[str, str]]:
        """Rationalize primer and template selections to minimize plate usage.
        
        Uses rationalization functions from primer_utils to optimize primer plate
        usage and template selections.
        
        Returns:
            Tuple containing:
                - Dictionary of rationalized primer selections
                - Dictionary of rationalized template selections
                
        Raises:
            ValueError: If primer locations or templates haven't been checked
        """
        if not self.primers_found and not self.missing_primers:
            raise ValueError("No primer location data. Please check primer locations first.")
        
        if not self.template_dict:
            raise ValueError("No template data. Please find templates first.")
            
        self.rationalized_primers = rationalize_primers(
            self.primers_found,
            self.missing_primers
        )
        
        self.rationalized_templates = rationalize_templates(self.template_dict)

        return self.rationalized_primers, self.rationalized_templates
    
    def write_instructions(self) -> List[List[str]]:
        """Generate assembly instructions for the TAR cloning experiment.
        
        Returns:
            List of instruction rows containing primer and template details
            
        Raises:
            ValueError: If primers and templates haven't been rationalized
        """
        if not self.rationalized_primers:
            raise ValueError("No rationalized primer data. Please rationalize selections first.")
            
        if not self.rationalized_templates:
            raise ValueError("No rationalized template data. Please rationalize selections first.")
            
        return write_circular_instructions(
            self.rationalized_primers,
            self.rationalized_templates,
            self.assembly_sequences,
            self.homology_length
        )
    
    def create_assembly(self) -> SeqRecord:
        """Create the assembled sequence with all parts and primers.
        
        Returns:
            SeqRecord object representing the assembled construct with features
            
        Raises:
            ValueError: If no primers have been designed
        """
        if not self.primers:
            raise ValueError("No primers available. Please design primers first.")
            
        if not self.assembly_sequences:
            raise ValueError("No sequences selected. Please select sequences first.")
            
        self.final_assembled_plasmid = assemble_parts_circular(
            self.assembly_sequences,
            self.primers,
            self.homology_length
        )

        return self.final_assembled_plasmid