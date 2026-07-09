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



from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyeast.utils.path_utils import get_integration_sites_path
from pyeast.utils.primer_utils import design_linear_primers, get_primer_locations, rationalize_primers
from pyeast.utils.pydna_utils import (
    annotate_linear_product,
    assemble_in_vivo,
    build_amplicons,
    cloning_history_json,
    diagnose_assembly,
    load_template_index,
    pcr_specificity,
)
from pyeast.utils.sequence_utils import get_templates, load_sequences, rationalize_templates, write_linear_instructions


@dataclass
class IntegrationResult:
    """Result of an integration cloning design run."""
    name: str
    assembly: SeqRecord
    primers: dict
    instructions: list
    missing_primers: dict = field(default_factory=dict)
    diagnostics: list = field(default_factory=list)
    specificity: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    provenance: object = None
    history: str = None

    def save(self, output_prefix: Path) -> None:
        """Save GenBank file and TSV outputs.

        Args:
            output_prefix: Path stem for output files, e.g. Path("output/my_construct/my_construct").
                Files are written as <output_prefix>.gb, <output_prefix>_instructions.tsv, etc.
        """
        output_prefix = Path(output_prefix)
        SeqIO.write(self.assembly, f"{output_prefix}.gb", "genbank")

        # Cloning history (OpenCloning CloningStrategy: primers + templates + how they
        # assembled) from pydna's assembly source DAG - additive.
        if self.history:
            with open(f"{output_prefix}_history.json", "w") as f:
                f.write(self.history)

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


class IntegrationDesigner:
    def __init__(self, homology_length: int = 25):
        self.homology_length = homology_length

        # State storage
        self.components = {}
        self.int_sites = {}  # all options for integration sites
        self.assembly_sequences = []
        self.int_site = None  # chosen integration site name
        self.primers = {}
        self.primers_found = {}
        self.missing_primers = {}
        self.template_dict = {}
        self.template_index = {}          # {template_name: SeqRecord} for pydna PCR
        self.rationalized_templates = {}
        self.rationalized_primers = {}
        self.final_assembly = None
        self.diagnostics = []             # Ambiguity diagnostic findings
        self.specificity = []             # PCR-specificity warnings
        self.assembly_notes = []          # Amplicon-build + assembly notes
        self.provenance = None            # pydna assembly source (provenance DAG)
        self.history = None               # OpenCloning CloningStrategy JSON
        self._amplicons = []              # Built pydna amplicons (internal)

    def load_sequences(self, components_dir: Path) -> None:
        """Load component sequences and integration sites."""
        self.components = load_sequences(components_dir)

        int_sites_dir = get_integration_sites_path()
        self.int_sites = self._load_int_sites(int_sites_dir)

        if not self.components:
            raise ValueError("No component sequences found")
        if not self.int_sites:
            raise ValueError("No integration sites found")

    def set_assembly_order(self, assembly_order: list[str], integration_site_name: str) -> None:
        """Set assembly sequences from component names and an integration site name.

        Args:
            assembly_order: Ordered list of component names (must be present in self.components).
            integration_site_name: Name of the integration site (must be present in self.int_sites).

        Raises:
            ValueError: If any name is not found in the loaded sequences.
        """
        invalid_components = [n for n in assembly_order if n not in self.components]
        if invalid_components:
            raise ValueError(f"Unknown component(s): {', '.join(invalid_components)}")

        matches = [s for s in self.int_sites if s.lower() == integration_site_name.lower()]
        if not matches:
            raise ValueError(f"Unknown integration site: {integration_site_name}")

        self.int_site = matches[0]
        upstream_seq, downstream_seq = self.int_sites[self.int_site]
        component_seqs = [self.components[n] for n in assembly_order]
        self.assembly_sequences = [upstream_seq] + component_seqs + [downstream_seq]

    def _load_int_sites(self, directory: Path = None) -> dict[str, tuple[SeqRecord, SeqRecord]]:
        """Load integration sites from both public and private directories.

        Each site should have two sequences: upstream and downstream homology.
        """
        if directory is None:
            directory = get_integration_sites_path()

        private_dir = get_integration_sites_path(private=True)

        int_sites = {}

        for search_dir in [directory, private_dir]:
            if not search_dir.exists():
                continue

            for file_path in search_dir.glob("*.fasta"):
                records = list(SeqIO.parse(file_path, "fasta"))
                if len(records) == 2:
                    site_name = file_path.stem
                    int_sites[site_name] = (records[0], records[1])
                else:
                    import warnings
                    warnings.warn(
                        f"Skipping {file_path.name} - expected 2 sequences, found {len(records)}"
                    )

        if not int_sites:
            raise FileNotFoundError(f"Integration sites directory not found: {directory}")

        return int_sites

    def design_integration_primers(self) -> dict[str, Seq]:
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

    def write_instructions(self) -> list[list[str]]:
        """Generate assembly instructions for the integration experiment.

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
        self.template_index = load_template_index(str(template_folder))

    def diagnose(self) -> list[dict]:
        """Run the pre-assembly ambiguity diagnostic and PCR-specificity analysis.

        Builds pydna amplicons for the full linear construct ([int_up] + components +
        [int_down]) and stores the ambiguity findings, specificity warnings, and
        amplicon-build notes on the designer.

        Returns:
            The ambiguity findings (empty list means a clean, unambiguous design).
        """
        if not self.primers:
            raise ValueError("No primers available. Please design primers first.")
        if not self.assembly_sequences:
            raise ValueError("No sequences selected. Please select sequences first.")

        self._amplicons, self.assembly_notes = build_amplicons(
            self.assembly_sequences,
            self.primers,
            self.rationalized_templates,
            self.template_index,
            self.homology_length,
        )
        self.diagnostics = diagnose_assembly(
            self._amplicons, limit=self.homology_length, overhang=self.homology_length
        )
        self.specificity = pcr_specificity(
            self.assembly_sequences,
            self.primers,
            self.rationalized_templates,
            self.template_index,
            self.homology_length,
        )
        return self.diagnostics

    def create_linear_assembly(self, proceed_on_ambiguity: bool = True) -> SeqRecord:
        """Create the assembled linear construct via pydna in-vivo-assembly.

        Simulates homologous recombination of the amplicons into a single linear molecule
        (the flanking integration-site primers carry no outer overhang, so the product stays
        linear) and adds component annotations. Degrades gracefully if the graph explodes.

        Args:
            proceed_on_ambiguity: When False, an over-dense (ambiguous) design re-raises
                instead of falling back to a constructed product.

        Returns:
            SeqRecord (pydna Dseqrecord) representing the assembled construct with features.

        Raises:
            ValueError: If no primers have been designed or no sequences selected.
        """
        if not self.primers:
            raise ValueError("No primers available. Please design primers first.")

        if not self.assembly_sequences:
            raise ValueError("No sequences selected. Please select sequences first.")

        if not getattr(self, "_amplicons", None):
            self.diagnose()

        product, assembly_notes = assemble_in_vivo(
            self._amplicons,
            limit=self.homology_length,
            circular_only=False,
            proceed_on_ambiguity=proceed_on_ambiguity,
        )
        self.assembly_notes = self.assembly_notes + assembly_notes
        product.name = "assembly"
        self.provenance = getattr(product, "source", None)
        self.history = cloning_history_json(product)
        #print(self.history)

        self.final_assembly = annotate_linear_product(
            product,
            self.assembly_sequences,
        )

        return self.final_assembly

    def design(
        self,
        components_dir: Path,
        assembly_order: list[str],
        integration_site_name: str,
        primer_folder: Path,
        template_folder: Path,
        name: str = "assembly",
        proceed_on_ambiguity: bool = True,
    ) -> IntegrationResult:
        """Design a chromosomal integration experiment programmatically.

        Runs the full design pipeline - load sequences, set assembly order, design primers,
        find templates, rationalize selections, generate instructions, and assemble.

        Args:
            components_dir: Path to the component library directory (FASTA files).
            assembly_order: Ordered list of component names to assemble (not including flanks).
            integration_site_name: Name of the integration site (must be in data/integration_sites).
            primer_folder: Path to primer plate library folder.
            template_folder: Path to template folder.
            name: Name for the assembled construct (default: "assembly").

        Returns:
            IntegrationResult containing the assembled construct, primers, instructions,
            and any missing primers. Call result.save(output_prefix) to write files.

        Example:
            designer = IntegrationDesigner(homology_length=25)
            result = designer.design(
                components_dir=Path("data/component_libraries/Saccharomyces_cerevisiae"),
                assembly_order=["GFP", "CYC1tt"],
                integration_site_name="HO",
                primer_folder=Path("data/primers"),
                template_folder=Path("data/templates"),
                name="GFP_integration",
            )
            result.save(Path("output/GFP_integration/GFP_integration"))
        """
        self.load_sequences(components_dir)
        self.set_assembly_order(assembly_order, integration_site_name)
        self.design_integration_primers()
        self.check_primer_locations(primer_folder)
        self.find_templates(template_folder)
        self.rationalize_selections()
        instructions = self.write_instructions()
        self.diagnose()
        assembly = self.create_linear_assembly(proceed_on_ambiguity=proceed_on_ambiguity)
        assembly.name = name

        return IntegrationResult(
            name=name,
            assembly=assembly,
            primers=self.primers,
            instructions=instructions,
            missing_primers=self.missing_primers,
            diagnostics=self.diagnostics,
            specificity=self.specificity,
            notes=self.assembly_notes,
            provenance=self.provenance,
            history=self.history,
        )
