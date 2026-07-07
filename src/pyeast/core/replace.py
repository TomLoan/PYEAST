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

# src/pyeast/core/replace.py

"""
Replace Designer for S. cerevisiae genome modifications.

This module provides tools for designing pop-in/pop-out replacements in S. cerevisiae.
"""
# ===========================================================================




from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from pyeast.utils.path_utils import get_component_libraries_path, get_templates_path
from pyeast.utils.primer_utils import design_screening_primers


@dataclass
class ReplaceResult:
    """Result of a replacement cassette design run."""
    name: str
    cassette: SeqRecord
    forward_primer: str
    reverse_primer: str
    product_sizes: dict
    genome_location: tuple  # (chrom_id, start, end, orientation)

    def save(self, output_prefix: Path) -> None:
        """Save GenBank, FASTA, and screening primers TSV.

        Args:
            output_prefix: Path stem, e.g. Path("output/my_replace/my_replace").
        """
        output_prefix = Path(output_prefix)
        SeqIO.write(self.cassette, f"{output_prefix}.gb", "genbank")
        SeqIO.write(self.cassette, f"{output_prefix}.fasta", "fasta")
        with open(f"{output_prefix}_screening_primers.tsv", "w") as f:
            f.write(f"{self.name}_ScreenF\t{self.forward_primer}\n")
            f.write(f"{self.name}_ScreenR\t{self.reverse_primer}")


class ReplaceDesigner:
    """Designer class for creating pop-in/pop-out replacement cassettes."""

    def __init__(self,
                 upstream_homology_len: int = 200,
                 downstream_homology_len: int = 200,
                 repeat_length: int = 160,
                 genome_file: Path | None = None,
                 ura3_file: Path | None = None):
        """Initialize the ReplaceDesigner.

        Args:
            upstream_homology_len: Length of upstream homology for recombination (default: 200)
            downstream_homology_len: Length of downstream homology for recombination (default: 200)
            repeat_length: Length of repeat sequence for marker removal (default: 160)
            genome_file: Path to the genome file (default: BY4741_Toronto_2012.fsa from data directory)
            ura3_file: Path to the URA3 marker file (default: URA3.fasta from data directory)
        """
        self.upstream_homology_len = upstream_homology_len
        self.downstream_homology_len = downstream_homology_len
        self.repeat_length = repeat_length

        self.genome_file = genome_file if genome_file is not None else get_templates_path() / "BY4741_Toronto_2012.fsa"
        self.ura3_file = ura3_file if ura3_file is not None else get_component_libraries_path() / "Saccharomyces_cerevisiae" / "URA3.fasta"

        # State storage
        self.target_sequence = None
        self.replacement_sequence = None
        self.genome_location = None  # (chrom_id, start, end, orientation)
        self.replacement_cassette = None
        self.screening_primers = None
        self.product_sizes = None
        self.marker_position = "upstream"  # Default position

    def find_target_sequence(self, genome_file: Path, target_seq: str) -> tuple[str, int, int, str] | None:
        """Locate a target sequence in the genome.

        Args:
            genome_file: Path to the genome file
            target_seq: The DNA sequence to find

        Returns:
            Tuple containing (chromosome_id, start, end, orientation) or None if not found
        """
        target_seq = Seq(target_seq.upper())
        for record in SeqIO.parse(genome_file, "fasta"):
            start = record.seq.find(target_seq)
            if start != -1:
                end = start + len(target_seq)
                return (record.id, start, end, "forward")

            start = record.seq.find(target_seq.reverse_complement())
            if start != -1:
                end = start + len(target_seq)
                return (record.id, start, end, "reverse")

        return None

    def extract_sequences(self, genome_seq: Seq, start: int, end: int, orientation: str) -> tuple[Seq, Seq, Seq]:
        """Extract required sequences based on marker position and orientation.

        Args:
            genome_seq: Full genome sequence
            start: Start position of target
            end: End position of target
            orientation: Either "forward" or "reverse"

        Returns:
            Tuple of (upstream_homology, downstream_homology, repeat)
        """
        if self.marker_position == "upstream" and orientation == "forward":
            upstream_homology = genome_seq[start :start + self.upstream_homology_len]
            downstream_homology = genome_seq[end:end + self.downstream_homology_len]
            repeat = genome_seq[start - self.repeat_length:start]

        elif self.marker_position == "upstream" and orientation == "reverse":
            upstream_homology = genome_seq[end-self.upstream_homology_len:end].reverse_complement()
            downstream_homology = genome_seq[start - self.downstream_homology_len:start].reverse_complement()
            repeat = genome_seq[end:end + self.repeat_length].reverse_complement()

        elif self.marker_position == "downstream" and orientation == "forward":
            upstream_homology = genome_seq[start - self.upstream_homology_len:start]
            downstream_homology = genome_seq[end - self.downstream_homology_len : end]
            repeat = genome_seq[end:end + self.repeat_length]

        else:  # downstream and reverse
            upstream_homology = genome_seq[end:end + self.upstream_homology_len].reverse_complement()
            downstream_homology = genome_seq[start : start + self.downstream_homology_len].reverse_complement()
            repeat = genome_seq[start - self.repeat_length:start].reverse_complement()

        return upstream_homology, downstream_homology, repeat

    def make_replacement_cassette(self, genome_file: Path, ura3_file: Path) -> SeqRecord:
        """Create the replacement cassette with URA3 marker.

        Ensures the replacement sequence is correctly oriented relative to the target
        sequence before assembly.

        Cassette structure with upstream marker:
        [upstream_homology]-[URA3]-[repeat]-[replacement]-[downstream_homology]

        Cassette structure with downstream marker:
        [upstream_homology]-[replacement]-[repeat]-[URA3]-[downstream_homology]
        """
        if not self.genome_location:
            raise ValueError("No target sequence location found")

        if not self.replacement_sequence:
            raise ValueError("No replacement sequence selected")

        chrom_id, start, end, orientation = self.genome_location

        # Get genome sequence
        for record in SeqIO.parse(genome_file, "fasta"):
            if record.id == chrom_id:
                genome = record
                break

        # Extract sequences with correct orientation
        upstream_homology, downstream_homology, repeat = self.extract_sequences(
            genome.seq, start, end, orientation
        )

        # Get URA3 marker
        ura3_marker = SeqIO.read(ura3_file, "fasta").seq

        # Orient replacement sequence correctly
        replacement_seq = self.replacement_sequence.seq
        if orientation == "reverse":
            replacement_seq = replacement_seq.reverse_complement()

        # Construct cassette based on marker position
        if self.marker_position == "upstream":
            cassette_seq = (upstream_homology + ura3_marker + repeat +
                        replacement_seq + downstream_homology)
        else:  # downstream
            cassette_seq = (upstream_homology + replacement_seq + repeat +
                        ura3_marker + downstream_homology)

        # Create SeqRecord
        cassette = SeqRecord(
            cassette_seq,
            id="replacement_cassette",
            name="replacement_cassette",
            description=f"Replacement cassette for {chrom_id}:{start}-{end}",
            annotations={
                "molecule_type": "DNA",
                "topology": "linear"
            }
        )

        # Add features with correct positions
        features = []
        current_pos = 0

        # Upstream homology
        features.append(SeqFeature(
            FeatureLocation(0, self.upstream_homology_len),
            type="misc_feature",
            qualifiers={"label": "upstream homology"}
        ))
        current_pos += self.upstream_homology_len

        if self.marker_position == "upstream":
            # URA3 marker
            features.append(SeqFeature(
                FeatureLocation(current_pos, current_pos + len(ura3_marker)),
                type="gene",
                qualifiers={"label": "URA3"}
            ))
            current_pos += len(ura3_marker)

            # Repeat
            features.append(SeqFeature(
                FeatureLocation(current_pos, current_pos + self.repeat_length),
                type="misc_feature",
                qualifiers={"label": "repeat"}
            ))
            current_pos += self.repeat_length

            # Replacement sequence
            features.append(SeqFeature(
                FeatureLocation(current_pos, current_pos + len(self.replacement_sequence.seq)),
                type="misc_feature",
                qualifiers={"label": "replacement"}
            ))
            current_pos += len(self.replacement_sequence.seq)

        else:  # downstream
            # Replacement sequence
            features.append(SeqFeature(
                FeatureLocation(current_pos, current_pos + len(self.replacement_sequence.seq)),
                type="misc_feature",
                qualifiers={"label": "replacement"}
            ))
            current_pos += len(self.replacement_sequence.seq)

            # Repeat
            features.append(SeqFeature(
                FeatureLocation(current_pos, current_pos + self.repeat_length),
                type="misc_feature",
                qualifiers={"label": "repeat"}
            ))
            current_pos += self.repeat_length

            # URA3 marker
            features.append(SeqFeature(
                FeatureLocation(current_pos, current_pos + len(ura3_marker)),
                type="gene",
                qualifiers={"label": "URA3"}
            ))
            current_pos += len(ura3_marker)

        # Downstream homology
        features.append(SeqFeature(
            FeatureLocation(current_pos, current_pos + self.downstream_homology_len),
            type="misc_feature",
            qualifiers={"label": "downstream homology"}
        ))

        cassette.features = features
        self.replacement_cassette = cassette
        return cassette

    def design_screening_strategy(self, genome_file: Path) -> tuple[str, str, dict[str, int]]:
        """Design screening primers and calculate expected product sizes.

        Args:
            genome_file: Path to the genome file

        Returns:
            Tuple containing (forward_primer, reverse_primer, product_sizes)

        Raises:
            ValueError: If no target sequence has been found
        """
        if not self.genome_location:
            raise ValueError("No target sequence location found")

        chrom_id, start, end, orientation = self.genome_location

        # Get relevant genomic sequence
        for record in SeqIO.parse(genome_file, "fasta"):
            if record.id == chrom_id:
                primer_start = max(0, start - 550)
                primer_end = min(len(record.seq), end + 550)
                relevant_sequence = str(record.seq[primer_start:primer_end])
                break

        # Calculate relative positions
        relative_start = 250
        relative_end = relative_start + (end - start) + 550

        # Design primers
        forward_primer, reverse_primer, product_size = design_screening_primers(
            relevant_sequence,
            relative_start,
            relative_end,
            product_size_range=(len(relevant_sequence) - 500, len(relevant_sequence))
        )

        if not all([forward_primer, reverse_primer, product_size]):
            raise ValueError("Failed to design screening primers")

        # Calculate expected product sizes
        replacement_size = len(self.replacement_sequence.seq)
        ura3_size = len(SeqIO.read(self.ura3_file, "fasta").seq)

        product_sizes = {
            'parent': product_size,
            'replacement_intermediate': (product_size - len(self.target_sequence) +
                                      replacement_size + self.repeat_length + ura3_size),
            'final_replacement': product_size - len(self.target_sequence) + replacement_size
        }

        self.screening_primers = (forward_primer, reverse_primer)
        self.product_sizes = product_sizes

        return forward_primer, reverse_primer, product_sizes

    def design(
        self,
        target_sequence: str,
        replacement_sequence: SeqRecord,
        marker_position: str = "upstream",
        name: str = "replacement_cassette",
        genome_location: tuple | None = None,
    ) -> ReplaceResult:
        """Design a replacement cassette programmatically.

        Args:
            target_sequence: DNA sequence to replace (A/T/G/C only).
            replacement_sequence: SeqRecord to insert in place of the target.
            marker_position: Where to place the URA3 marker - "upstream" or "downstream"
                of the replacement sequence (default: "upstream").
            name: Name for the output construct (default: "replacement_cassette").
            genome_location: Optional pre-computed (chrom_id, start, end, orientation) tuple.
                If provided, skips the genome search step (useful when displaying
                the location interactively before committing to the design).

        Returns:
            ReplaceResult containing the cassette, screening primers, product sizes,
            and genome location. Call result.save(output_prefix) to write files.

        Raises:
            ValueError: If target_sequence is not found in the genome, or marker_position
                is not "upstream" or "downstream".

        Example:
            from Bio import SeqIO
            replacement = SeqIO.read("pTEF1.fasta", "fasta")

            designer = ReplaceDesigner(upstream_homology_len=200)
            result = designer.design(
                target_sequence="ATGCATGC...",
                replacement_sequence=replacement,
                marker_position="upstream",
                name="YFG1_pTEF1_replace",
            )
            result.save(Path("output/YFG1_pTEF1_replace/YFG1_pTEF1_replace"))
        """
        if marker_position not in ("upstream", "downstream"):
            raise ValueError(f"marker_position must be 'upstream' or 'downstream', got: {marker_position!r}")

        self.target_sequence = target_sequence.upper()
        self.replacement_sequence = replacement_sequence
        self.marker_position = marker_position

        if genome_location is not None:
            self.genome_location = genome_location
        else:
            location = self.find_target_sequence(self.genome_file, self.target_sequence)
            if location is None:
                raise ValueError(f"Target sequence not found in genome: {self.genome_file}")
            self.genome_location = location

        cassette = self.make_replacement_cassette(self.genome_file, self.ura3_file)
        cassette.name = name

        fwd, rev, product_sizes = self.design_screening_strategy(self.genome_file)

        return ReplaceResult(
            name=name,
            cassette=cassette,
            forward_primer=fwd,
            reverse_primer=rev,
            product_sizes=product_sizes,
            genome_location=self.genome_location,
        )
