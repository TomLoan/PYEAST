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




from pathlib import Path
from typing import Optional, Tuple, Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from rich.console import Console
from rich.table import Table
from prompt_toolkit import PromptSession
from prompt_toolkit.completion import WordCompleter

from ..utils.primer_utils import design_screening_primers
from ..utils.sequence_utils import load_sequences
from ..utils.visualisation import visualise_genbank, save_figure

class ReplaceDesigner:
    """Designer class for creating pop-in/pop-out replacement cassettes."""
    
    def __init__(self, 
                 upstream_homology_len: int = 200,
                 downstream_homology_len: int = 200,
                 repeat_length: int = 160,
                 genome_file: Path = Path("data/templates/BY4741_Toronto_2012.fsa"),
                 ura3_file: Path = Path("data/component libraries/Saccharomyces cerevisiae/URA3.fasta")):
        """Initialize the ReplaceDesigner.
        
        Args:
            upstream_homology_len: Length of upstream homology for recombination (default: 200)
            downstream_homology_len: Length of downstream homology for recombination (default: 200)
            repeat_length: Length of repeat sequence for marker removal (default: 160)
            genome_file: Path to the genome file
            ura3_file: Path to the URA3 marker file
        """
        self.upstream_homology_len = upstream_homology_len
        self.downstream_homology_len = downstream_homology_len
        self.repeat_length = repeat_length
        self.genome_file = genome_file
        self.ura3_file = ura3_file
        self.console = Console()
        self.session = PromptSession()
        
        # State storage
        self.target_sequence = None
        self.replacement_sequence = None
        self.genome_location = None  # (chrom_id, start, end, orientation)
        self.replacement_cassette = None
        self.screening_primers = None
        self.product_sizes = None
        self.marker_position = "upstream"  # Default position
        
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

    def print_sequence_grid(self, sequences: Dict[str, SeqRecord], title: str = "Available Sequences"):
        """Display available sequences in a formatted table."""
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

    def get_replacement_selection(self, components_dir: Path) -> Optional[SeqRecord]:
        """Get user selection for replacement sequence from a library.
        
        Args:
            components_dir: Directory containing sequence libraries
            
        Returns:
            Selected sequence as SeqRecord or None if selection cancelled
        """
        # Load sequences
        sequences = load_sequences(components_dir)
        if not sequences:
            self.console.print("[red]No sequences found in directory[/red]")
            return None
            
        # Display sequences and create completer
        self.print_sequence_grid(sequences)
        sequence_completer = WordCompleter(list(sequences.keys()), ignore_case=True)
        
        while True:
            try:
                self.console.print("\n[blue]Select a sequence to insert:[/blue]")
                user_input = self.session.prompt(
                    "Sequence name: ",
                    completer=sequence_completer
                ).strip()
                
                # Find matching sequence
                matches = [seq for name, seq in sequences.items() 
                         if name.lower() == user_input.lower()]
                
                if matches:
                    return matches[0]
                else:
                    self.console.print("[red]Invalid sequence name[/red]")
                    
            except KeyboardInterrupt:
                return None
            
    def extract_sequences(self, genome_seq: Seq, start: int, end: int, orientation: str) -> Tuple[Seq, Seq, Seq]:
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
            repeat = genome_seq[start - self.repeat_length:
                            start]
                            
        elif self.marker_position == "upstream" and orientation == "reverse":
            upstream_homology = genome_seq[end-self.upstream_homology_len:end].reverse_complement()
            downstream_homology = genome_seq[start - self.downstream_homology_len:start].reverse_complement()
            repeat = genome_seq[end :
                            end + self.repeat_length].reverse_complement()
                            
        elif self.marker_position == "downstream" and orientation == "forward":
            upstream_homology = genome_seq[start - self.upstream_homology_len:start]
            downstream_homology = genome_seq[end - self.downstream_homology_len : end]
            repeat = genome_seq[end:
                            end + self.repeat_length]
                            
        else:  # downstream and reverse
            upstream_homology = genome_seq[end:end + self.upstream_homology_len].reverse_complement()
            downstream_homology = genome_seq[start : start + self.downstream_homology_len].reverse_complement()
            repeat = genome_seq[start - self.repeat_length:
                            start].reverse_complement()
        
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