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

import logging
import os
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .path_utils import get_private_equivalent

# Marker recorded in a feature's /note qualifier to identify a PYEAST component part.
# The feature type itself is the standard "misc_feature" so the annotation survives
# round-trips through GenBank tools (SnapGene / Geneious / OpenCloning) that coerce
# non-standard feature keys. Constructs written before the changeover used a custom
# type="PYEAST_component"; that is still recognised for backward compatibility.
PYEAST_COMPONENT_NOTE = "PYEAST_component"


def is_pyeast_component(feature) -> bool:
    """True if a SeqFeature marks a PYEAST component part.

    Recognises the current form (type="misc_feature" carrying a /note="PYEAST_component"
    marker) and the legacy form (custom type="PYEAST_component").
    """
    if feature.type == "PYEAST_component":
        return True
    return PYEAST_COMPONENT_NOTE in feature.qualifiers.get("note", [])


def get_component_features(record) -> list:
    """Return a record's PYEAST component features, most-preferred detection first.

    Prefers features carrying the PYEAST marker (new note-marked misc_feature, or the
    legacy custom type). Only when a record has no marked components does it fall back to
    treating every bare ``misc_feature`` as a component - the legacy path for externally
    supplied GenBank files that predate the marker. This keeps template-derived
    misc_features (which now flow through into outputs) from being mistaken for parts.
    """
    marked = [f for f in record.features if is_pyeast_component(f)]
    if marked:
        return marked
    return [f for f in record.features if f.type == "misc_feature"]


def load_sequences(directory: str) -> dict[str, SeqRecord]:
    """
    Load all sequences from FASTA files in both public and private directories.
    Private sequences override public ones if they have the same name.
    Can handle directories that exist only in public, only in private, or both.

    Args:
        directory (str): Path to the directory containing FASTA files.

    Returns:
        Dict[str, Sequence]: Dictionary of loaded sequences.

    Raises:
        FileNotFoundError: If neither public nor private directory is found.
    """
    sequences = {}

    # Convert to Path for easier manipulation
    public_dir = Path(directory)

    # Construct private directory path (if path is within data directory)
    try:
        private_dir = get_private_equivalent(public_dir)
    except ValueError:
        # Path is not within data directory, so no private equivalent exists
        private_dir = None

    # Track if we found at least one directory
    found_any = False

    # Load from public directory if it exists
    if public_dir.exists():
        found_any = True
        try:
            for filename in os.listdir(public_dir):
                if filename.endswith(('.fasta', '.fa', 'fsa')):
                    file_path = public_dir / filename
                    for record in SeqIO.parse(file_path, "fasta"):
                        name = record.id.split()[0]
                        sequences[name] = record
        except Exception as e:
            logging.error(f"Error loading sequences from {public_dir}: {str(e)}")
            raise

    # Load from private directory if it exists (overwrites public if same name)
    if private_dir and private_dir.exists():
        found_any = True
        try:
            for filename in os.listdir(private_dir):
                if filename.endswith(('.fasta', '.fa', 'fsa')):
                    file_path = private_dir / filename
                    for record in SeqIO.parse(file_path, "fasta"):
                        name = record.id.split()[0]
                        sequences[name] = record
        except Exception as e:
            logging.error(f"Error loading sequences from private directory: {str(e)}")
            raise

    # If neither directory exists, raise error
    if not found_any:
        if private_dir:
            logging.error(f"Directory not found: {public_dir} (also checked {private_dir})")
        else:
            logging.error(f"Directory not found: {public_dir}")
        raise FileNotFoundError(f"Directory not found: {public_dir}")

    return sequences

def get_templates(parts: list[SeqRecord], directory: str) -> dict[str, list[str]]:
    """
    Find template matches for given parts in both public and private directories.
    Can handle directories that exist only in public, only in private, or both.

    This function searches for templates in FASTA or GenBank files that contain
    the sequences of the provided parts.

    Args:
        parts (List[SeqRecord]): A list of SeqRecord objects representing the parts.
        directory (str): The directory path containing template files.

    Returns:
        Dict[str, List[str]]: A dictionary mapping part IDs to lists of matching template IDs.
    """
    templates = {}

    # Convert to Path for easier manipulation
    public_dir = Path(directory)

    # Construct private directory path (if path is within data directory)
    try:
        private_dir = get_private_equivalent(public_dir)
    except ValueError:
        # Path is not within data directory, so no private equivalent exists
        private_dir = None

    # Load templates from both public and private directories
    search_dirs = [public_dir] + ([private_dir] if private_dir else [])
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue

        for template_file in os.listdir(search_dir):
            file_path = search_dir / template_file
            if str(file_path).endswith(('.fasta', '.fsa', '.fa', '.gb', '.gbk')):
                format_type = "fasta" if str(file_path).endswith(('.fasta', '.fa', '.fsa')) else "genbank"
                for record in SeqIO.parse(file_path, format_type):
                    templates[record.name] = record.seq

    # Find matching templates for each part
    templates_used = {}
    for part in parts:
        templates_used[part.name] = []
        for template_id, template_seq in templates.items():
            if part.seq.upper() in template_seq.upper() or part.seq.upper() in template_seq.upper().reverse_complement():
                templates_used[part.name].append(template_id)

        if not templates_used[part.name]:
            templates_used[part.name] = ["Not found"]

    return templates_used

def rationalize_templates(template_dict: dict[str, list[str]],
                          preferred_templates: list[str] | None = None) -> dict[str, str]:
    """
    Rationalize template selection to minimize the number of unique templates.

    This function chooses templates based on the following priority:
    1. Templates in the preferred list
    2. Highest frequency of use across all parts
    3. Shorter names

    Args:
        template_dict (Dict[str, List[str]]): A dictionary mapping part names to lists of potential templates.
        preferred_templates (Optional[List[str]]): Templates to prefer when a part matches
            several. Matched against the template record name. When None, the list is read
            from the user config (empty by default, so the shortest-match sorter decides).

    Returns:
        Dict[str, str]: A dictionary mapping part names to their choseninput template.
    """
    # Load the preferred list from config unless one was passed explicitly.
    if preferred_templates is None:
        from pyeast.config import get_config
        preferred_templates = get_config().preferred_templates

    # Count the global frequency of each template

    all_templates = [template for templates in template_dict.values() for template in templates if template != "Not found"]
    template_frequency = Counter(all_templates)
    # print(all_templates)
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

def parse_gb_file(file_path: str) -> tuple[SeqRecord, list[SeqRecord]]:
    """
    Parse a GenBank file and extract the full plasmid sequence and its parts.

    Args:
    file_path (str): Path to the GenBank file

    Returns:
    Tuple[SeqRecord, List[SeqRecord]]: Full plasmid sequence and a list of parts
    """
    plasmid = SeqIO.read(file_path, "genbank")
    parts = []
    for feature in get_component_features(plasmid):
        part_name = feature.qualifiers.get("label", [f"Part_{len(parts)}"])[0]
        part_seq = feature.extract(plasmid.seq)
        parts.append(SeqRecord(part_seq, id=part_name, name=part_name, description=""))
    return plasmid, parts

def write_circular_instructions(rationalized_primers: dict[str, dict],
                              rationalized_templates: dict[str, str],
                              assembly_sequences: list[SeqRecord],
                              homology_length: int) -> list[list[str]]:
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
        else:
            print(f"primer(s) missing for {part_name}")

    return instructions

def write_linear_instructions(rationalized_primers: dict[str, dict],
                            rationalized_templates: dict[str, str],
                            int_site_up: SeqRecord,
                            middle_sequences: list[SeqRecord],
                            int_site_down: SeqRecord,
                            homology_length: int) -> list[list[str]]:
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

def find_matching_primer(primers: dict[str, dict], part_seq: Seq, is_forward: bool, homology_length: int) -> dict | None:
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

# NB: assemble_parts_circular / assemble_parts_linear (blind concatenation + primer
# annotation) and their find_all_occurrences helper have been removed. Assembly is now a
# real recombination simulation via pydna - see pyeast.utils.pydna_utils and the
# create_assembly / create_linear_assembly methods in core/tar.py and core/integration.py.

