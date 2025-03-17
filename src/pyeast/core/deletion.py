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

# src/pyeast/core/deletion.py
"""
Deletion Designer for S. cerevisiae genome modifications.

This module provides tools for designing scarless deletion cassettes in S. cerevisiae.
"""

# ===========================================================================





from pathlib import Path
from typing import Optional, Tuple, Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from rich.console import Console
from rich.table import Table
from prompt_toolkit import PromptSession

from ..utils.primer_utils import design_screening_primers
from ..utils.visualisation import visualise_genbank, save_figure

class DeletionDesigner:
    """Designer class for creating scarless deletion cassettes."""
    
    def __init__(self, 
                 upstream_homology_len: int = 300,
                 downstream_homology_len: int = 200,
                 repeat_length: int = 160,
                 genome_file: Path = Path("data/templates/BY4741_Toronto_2012.fsa"),
                 ura3_file: Path = Path("data/component libraries/Saccharomyces cerevisiae/URA3.fasta")):
        """Initialize the DeletionDesigner.
        
        Args:
            upstream_homology_len: Length of upstream homology for recombination (default: 300)
            downstream_homology_len: Length of downstream homology for recombination (default: 200)
            repeat_length: Length of repeat sequence for marker removal (default: 160)
            genome_file: Path to the genome file (default: BY4741_Toronto_2012.fsa)
            ura3_file: Path to the URA3 marker file
        """
        self.upstream_homology_len = upstream_homology_len  # keeping variable name for compatibility
        self.downstream_homology_len = downstream_homology_len
        self.repeat_length = repeat_length
        self.genome_file = genome_file
        self.ura3_file = ura3_file
        self.console = Console()
        self.session = PromptSession()
        
        # State storage
        self.target_sequence = None
        self.genome_location = None  # (chrom_id, start, end, orientation)
        self.deletion_cassette = None
        self.screening_primers = None
        self.product_sizes = None
        
    def find_target_sequence(self, genome_file: Path, target_seq: str) -> Optional[Tuple[str, int, int, str]]:
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
    
    def extract_sequences(self, genome_seq: Seq, start: int, end: int, orientation: str) -> Tuple[Seq, Seq, Seq]:
        """Extract required sequences with correct orientation.
        
        Args:
            genome_seq: Full genome sequence
            start: Start position of target
            end: End position of target
            orientation: Either "forward" or "reverse"
            
        Returns:
            Tuple of (upstream_homology, downstream_homology, repeat)
        """
        if orientation == "forward":
            upstream_homology = genome_seq[start-self.upstream_homology_len:start]
            downstream_homology = genome_seq[end-self.downstream_homology_len:end]
            repeat = genome_seq[end :end + self.repeat_length]
        else:  # reverse
            upstream_homology = genome_seq[end:end+self.upstream_homology_len].reverse_complement()
            downstream_homology = genome_seq[start : start + self.downstream_homology_len].reverse_complement()
            repeat = genome_seq[start - self.repeat_length:start].reverse_complement()
        
        return upstream_homology, downstream_homology, repeat

        
    def make_deletion_cassette(self, genome_file: Path, ura3_file: Path) -> SeqRecord:
        """Create the deletion cassette with URA3 marker.
        
        Args:
            genome_file: Path to the genome file
            ura3_file: Path to the URA3 marker file
            
        Returns:
            SeqRecord object representing the deletion cassette
            
        Raises:
            ValueError: If no target sequence has been found
        """
        if not self.genome_location:
            raise ValueError("No target sequence location found. Run find_target_sequence first.")
            
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
        
        # Construct cassette
        cassette_seq = upstream_homology + repeat + ura3_marker + downstream_homology
        
        # Create SeqRecord
        cassette = SeqRecord(
            cassette_seq,
            id="deletion_cassette",
            name="deletion_cassette",
            description=f"Deletion cassette for {chrom_id}:{start}-{end}",
            annotations={
                "molecule_type": "DNA",
                "topology": "linear"
            }
        )
        
        # Add features
        features = [
            SeqFeature(
                FeatureLocation(0, self.upstream_homology_len),
                type="misc_feature",
                qualifiers={"label": "upstream homology"}
            ),
            SeqFeature(
                FeatureLocation(self.upstream_homology_len, 
                            self.upstream_homology_len + self.repeat_length),
                type="misc_feature",
                qualifiers={"label": "repeat"}
            ),
            SeqFeature(
                FeatureLocation(self.upstream_homology_len + self.repeat_length,
                            self.upstream_homology_len + self.repeat_length + len(ura3_marker)),
                type="gene",
                qualifiers={"label": "URA3"}
            ),
            SeqFeature(
                FeatureLocation(self.upstream_homology_len + self.repeat_length + len(ura3_marker),
                            len(cassette_seq)),
                type="misc_feature",
                qualifiers={"label": "downstream homology"}
            )
        ]
        
        cassette.features = features
        self.deletion_cassette = cassette
        return cassette
        
    def design_screening_strategy(self, genome_file: Path) -> Tuple[str, str, Dict[str, int]]:
        """Design screening primers and calculate expected product sizes.
        
        Args:
            genome_file: Path to the genome file
            
        Returns:
            Tuple containing (forward_primer, reverse_primer, product_sizes)
            
        Raises:
            ValueError: If no target sequence has been found
        """
        if not self.genome_location:
            raise ValueError("No target sequence location found. Run find_target_sequence first.")
            
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
        product_sizes = {
            'parent': product_size,
            'deletion_intermediate': product_size - len(self.target_sequence) + self.repeat_length  + len(SeqIO.read(self.ura3_file, "fasta").seq),
            'final_deletion': product_size - len(self.target_sequence)
        }
        
        self.screening_primers = (forward_primer, reverse_primer)
        self.product_sizes = product_sizes
        
        return forward_primer, reverse_primer, product_sizes
        