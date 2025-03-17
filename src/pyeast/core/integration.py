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

# src/pyeast/core/integration.py
"""
Integration designer for integrations of DNA sequences into the genome of S. cerevisiae.
"""

# ===========================================================================



from pathlib import Path
from typing import Dict, List, Tuple
from Bio.SeqRecord import SeqRecord
from rich.console import Console
from rich.table import Table
from prompt_toolkit import PromptSession
from prompt_toolkit.completion import WordCompleter
from prompt_toolkit.shortcuts import confirm
from Bio.SeqRecord import SeqRecord 
from Bio import SeqIO
from Bio.Seq import Seq

from ..utils.sequence_utils import (
    load_sequences,
    get_templates,
    rationalize_templates,
    assemble_parts_linear,
    write_linear_instructions
)
from ..utils.primer_utils import (
    design_linear_primers,
    get_primer_locations,
    rationalize_primers
)
from ..utils.visualisation import visualise_genbank, save_figure

class IntegrationDesigner:
    def __init__(self, homology_length: int = 25):
        self.homology_length = homology_length
        self.console = Console()
        self.session = PromptSession()
        
        # State storage
        self.components = {}
        self.int_sites = {} #all opitoins for integration sites
        self.assembly_sequences = []
        self.int_site = None #chosen integration site
        self.primers = {}
        self.primers_found = {}
        self.missing_primers = {}
        self.template_dict = {}
        self.rationalized_templates = {}
        self.rationalized_primers = {}
        self.final_assembly = None

    def load_sequences(self, components_dir: Path) -> None:
        """Load component sequences and integration sites."""
        # Load components
        self.components = load_sequences(components_dir)
        
        # Load integration sites from standard location
        int_sites_dir = Path("data/integration sites")
        self.int_sites = self._load_int_sites(int_sites_dir)
        
        if not self.components:
            raise ValueError("No component sequences found")
        if not self.int_sites:
            raise ValueError("No integration sites found")


    def get_assembly_selections(self) -> None:
        """Get user selections and return ordered list of sequences for assembly."""
        component_completer = WordCompleter(list(self.components.keys()), ignore_case=True)
        int_site_completer = WordCompleter(list(self.int_sites.keys()), ignore_case=True)

        while True:
            try:
                # Print available sequences first
                self.print_sequence_grid(self.components, "Available Components")
                self.print_integration_sites(self.int_sites,"\n[bold cyan]Available Integration Sites:[/bold cyan]")
                # for site in self.int_sites.keys():
                #     self.console.print(f"  {site}")

                # Ask for components
                user_input = self.session.prompt(
                    "\nEnter the names of the components you want to assemble, in order (space-separated): ",
                    completer=component_completer
                )
                
                self.console.print(f"Got user input: {user_input}")  # Debug print
                
                # Get component sequences
                component_seqs = []
                for name in user_input.split():
                    matches = [comp for comp in self.components.keys() 
                            if comp.lower() == name.lower()]
                    if matches:
                        component_seqs.append(self.components[matches[0]])
                    else:
                        raise ValueError(f"Invalid component name: {name}")

                self.console.print(f"Found {len(component_seqs)} component sequences")  # Debug print

                # Get integration site
                while True:
                    int_site = self.session.prompt(
                        "\nEnter the name of the integration site: ",
                        completer=int_site_completer
                    )
                    
                    matches = [site for site in self.int_sites.keys() 
                            if site.lower() == int_site.lower()]
                    if matches:
                        upstream_seq, downstream_seq = self.int_sites[matches[0]]
                        self.console.print(f"Found integration site: {matches[0]}")  # Debug print
                        break
                    else:
                        self.console.print(f"[red]Invalid integration site name: {int_site}[/red]")

                # Show selection and confirm
                self.console.print("\n[green]Selected assembly:[/green]")
                self.console.print("Components:")
                for i, seq in enumerate(component_seqs, 1):
                    self.console.print(f"{i}. {seq.id}")
                self.console.print(f"Integration site: {int_site}")

                if confirm("Is this correct?"):
                    assembled_sequences = [upstream_seq] + component_seqs + [downstream_seq]
                    #self.console.print(f"Returning {len(assembled_sequences)} sequences")  # Debug print
                    self.assembly_sequences = assembled_sequences
                    return
                else:
                    self.console.print("Let's try again.")
                    continue

            except ValueError as e:
                self.console.print(f"[red]Error: {e}[/red]")
            except KeyboardInterrupt:
                if confirm("\nCancel selection?"):
                    raise KeyboardInterrupt
                
    def _load_int_sites(self, directory: Path = None) -> Dict[str, Tuple[SeqRecord, SeqRecord]]:
        """Load integration sites from data/integration sites directory.
        
        Each site should have two sequences: upstream and downstream homology.
        """
        if directory is None:
            directory = Path("data/integration sites")
        
        int_sites = {}
        if not directory.exists():
            self.console.print(f"[red]Integration sites directory not found: {directory}[/red]")
            raise FileNotFoundError(f"Integration sites directory not found: {directory}")
            
        for file_path in directory.glob("*.fasta"):
            records = list(SeqIO.parse(file_path, "fasta"))
            if len(records) == 2:
                site_name = file_path.stem
                int_sites[site_name] = (records[0], records[1])
            else:
                self.console.print(f"[yellow]Warning: Skipping {file_path.name} - expected 2 sequences, found {len(records)}[/yellow]")
            
        return int_sites
    
    def design_integration_primers(self) -> Dict[str, Seq]: 
        """Design integration primers with the necesary overhangs
        
        uses the design_linear_primers function from primer_utils to design primers 
        for the selected assembly sequences. 
        
        Returns: 
            Dictionary mapping primer names to their sequence
            
        Raises: 
            ValueError: If no sequences have been selected for assembly. 
        """

        if not self.assembly_sequences: 
            raise ValueError("No sequences selected for assembly. Please select sequences first")
        
        
        self.primers = design_linear_primers(
                self.assembly_sequences[1:-1],  # middle components
                (self.assembly_sequences[0], self.assembly_sequences[-1]),  # int sites
                target_tm=50,
                homology_length=self.homology_length
            )
        return self.primers
                
    def print_sequence_grid(self, sequences: Dict[str, SeqRecord], title: str = "Available Sequences"):
        """Print the available sequences in a grid format.
        
        Args:
            sequences: Dictionary of sequences to display
            title: Title for the sequence table
        """
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

    def print_integration_sites(self, sites: Dict[str, Tuple[SeqRecord, SeqRecord]], title: str = "Available Integration Sites"):
        """Display integration sites in a formatted table.
        
        Args:
            sites: Dictionary mapping site names to (upstream, downstream) sequence tuples
            title: Title for the table
        """
        table = Table(title=title)
        table.add_column("Name", style="cyan")
        table.add_column("Up", justify="right", style="green")
        table.add_column("Down", justify="right", style="green")
        table.add_column("Description", style="white")
        
        for name, (up_seq, down_seq) in sites.items():
            table.add_row(
                name,
                f"{len(up_seq)} bp",
                f"{len(down_seq)} bp",
                up_seq.description[:150] + "..." if len(up_seq.description) > 149 else up_seq.description
            )
        
        self.console.print(table)
    
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
        
        return write_linear_instructions(
            self.rationalized_primers,
            self.rationalized_templates,
            self.assembly_sequences[0],
            self.assembly_sequences[1:-1],
            self.assembly_sequences[-1],
            self.homology_length
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
        #self.console.print(self.template_dict)

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


    def create_linear_assembly(self) -> SeqRecord:
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
            
        self.final_assembly = assemble_parts_linear(
            self.assembly_sequences,
            self.primers,
        )

        return self.final_assembly
