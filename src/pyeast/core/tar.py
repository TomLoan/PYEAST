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



from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyeast.utils.primer_utils import design_circular_primers, get_primer_locations, rationalize_primers
from pyeast.utils.sequence_utils import assemble_parts_circular, get_templates, load_sequences, rationalize_templates, write_circular_instructions


@dataclass
class TARResult:
    """Result of a TAR cloning design run."""
    name: str
    assembly: SeqRecord
    primers: dict
    instructions: list
    missing_primers: dict = field(default_factory=dict)

    def save(self, output_prefix: Path) -> None:
        """Save GenBank file and TSV outputs.

        Args:
            output_prefix: Path stem for output files, e.g. Path("output/my_plasmid/my_plasmid").
                Files are written as <output_prefix>.gb, <output_prefix>_instructions.tsv, etc.
        """
        output_prefix = Path(output_prefix)
        SeqIO.write(self.assembly, f"{output_prefix}.gb", "genbank")

        with open(f"{output_prefix}_all_primers.tsv", "w") as f:
            f.write("Name\tSequence\n")
            for name, seq in self.primers.items():
                f.write(f"{name}\t{seq}\n")

        if self.missing_primers:
            with open(f"{output_prefix}_missing_primers.tsv", "w") as f:
                f.write("Name\tSequence\n")
                for name, info in self.missing_primers.items():
                    f.write(f"{name}\t{info[0]['sequence']}\n")

        with open(f"{output_prefix}_instructions.tsv", "w") as f:
            f.write("Part\tF_Plate\tF_Name\tF_Well\tR_Plate\tR_Name\tR_Well\tTemplate\tSize\n")
            for row in self.instructions:
                f.write("\t".join(map(str, row)) + "\n")


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

    def set_assembly_order(self, selected_names: list[str])-> None:
        """Set the assmebly order from user selected seqeunce names"""

        self.assembly_sequences = [
            self.available_sequences[name] for name in selected_names
        ]

    def design_tar_primers(self) -> dict[str, Seq]:
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


    def rationalize_selections(self) -> tuple[dict[str, dict], dict[str, str]]:
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

    def write_instructions(self) -> list[list[str]]:
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

    def design(
        self,
        library_path: Path,
        assembly_order: list[str],
        primer_folder: Path,
        template_folder: Path,
        name: str = "assembly",
    ) -> TARResult:
        """Design a TAR cloning experiment programmatically.

        Runs the full design pipeline - load sequences, design primers, find templates,
        rationalize selections, generate instructions, and assemble.

        Args:
            library_path: Path to the component library directory (FASTA files).
            assembly_order: Ordered list of component names to assemble.
            primer_folder: Path to primer plate library folder.
            template_folder: Path to template folder.
            name: Name for the assembled construct (default: "assembly").

        Returns:
            TARResult containing the assembled construct, primers, instructions,
            and any missing primers. Call result.save(output_prefix) to write files.

        Example:
            designer = TARDesigner(homology_length=25)
            result = designer.design(
                library_path=Path("data/component_libraries/YeastToolKit"),
                assembly_order=["part1", "part2", "part3"],
                primer_folder=Path("data/primers"),
                template_folder=Path("data/templates"),
                name="my_plasmid",
            )
            result.save(Path("output/my_plasmid/my_plasmid"))
        """
        self.load_and_get_sequences(library_path)
        self.set_assembly_order(assembly_order)
        self.design_tar_primers()
        self.check_primer_locations(primer_folder)
        self.find_templates(template_folder)
        self.rationalize_selections()
        instructions = self.write_instructions()
        assembly = self.create_assembly()
        assembly.name = name

        return TARResult(
            name=name,
            assembly=assembly,
            primers=self.primers,
            instructions=instructions,
            missing_primers=self.missing_primers,
        )
