"""Tool schemas and handlers for the PYEAST agent."""

import json
import shutil
import tempfile
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path

from pyeast.utils.path_utils import (
    ensure_output_dir_exists,
    get_component_libraries_path,
    get_integration_sites_path,
    get_output_path,
    get_primers_path,
    get_templates_path,
)

_SEQ_EXTENSIONS = {".fasta", ".fa", ".gb", ".genbank", ".gbk"}

# ---------------------------------------------------------------------------
# Tool schemas
# ---------------------------------------------------------------------------

TOOL_SCHEMAS = [
    {
        "name": "list_components",
        "description": "List available PYEAST libraries, parts, integration sites, or templates.",
        "input_schema": {
            "type": "object",
            "properties": {
                "component_type": {
                    "type": "string",
                    "enum": ["libraries", "components", "integration_sites", "templates"],
                    "description": (
                        "'libraries': library folder names. "
                        "'components': parts in a library (requires library_name). "
                        "'integration_sites': genomic insertion loci. "
                        "'templates': PCR template files."
                    ),
                },
                "library_name": {
                    "type": "string",
                    "description": "Library name; required when component_type='components'.",
                },
            },
            "required": ["component_type"],
        },
    },
    {
        "name": "lookup_gene_sequence",
        "description": (
            "Fetch a S. cerevisiae gene sequence from SGD by name. "
            "Call before design_deletion or design_replacement when the user names a gene."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene_name": {
                    "type": "string",
                    "description": "Standard or systematic gene name (e.g. 'PEP4', 'YPL154C').",
                },
                "sequence_type": {
                    "type": "string",
                    "enum": ["orf", "genomic", "upstream_1kb"],
                    "description": (
                        "'orf': coding sequence — use for deletions. "
                        "'genomic': ORF + UTRs. "
                        "'upstream_1kb': 1 kb upstream region — use for promoter replacements."
                    ),
                },
            },
            "required": ["gene_name"],
        },
    },
    {
        "name": "read_component",
        "description": (
            "Read the DNA sequence of an existing library part. "
            "Use this when you need to concatenate multiple parts into a single replacement sequence, "
            "or to inspect what a part contains before using it in a design."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {"type": "string", "description": "Part name (without file extension)."},
                "library_name": {"type": "string", "description": "Library folder name."},
            },
            "required": ["name", "library_name"],
        },
    },
    {
        "name": "save_component",
        "description": (
            "Save a user-provided DNA sequence as a FASTA part in a component library, "
            "making it immediately available for design_tar and design_integration."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {
                    "type": "string",
                    "description": "Part name (no spaces; e.g. 'HFBI_codon_opt'). Becomes the filename.",
                },
                "sequence": {
                    "type": "string",
                    "description": "DNA sequence — raw A/T/G/C string, or FASTA format with a >header line.",
                },
                "library_name": {
                    "type": "string",
                    "description": "Library folder to save into. Must match library_name used in the design call.",
                },
                "description": {
                    "type": "string",
                    "description": "Short description written into the FASTA header (optional).",
                },
                "overwrite": {
                    "type": "boolean",
                    "description": "If true, replace an existing part with the same name. Use to update a description or sequence.",
                },
            },
            "required": ["name", "sequence", "library_name"],
        },
    },
    {
        "name": "list_outputs",
        "description": "List constructs already saved to the output directory. Call before run_batch.",
        "input_schema": {
            "type": "object",
            "properties": {},
        },
    },
    {
        "name": "design_tar",
        "description": (
            "Assemble library parts into a circular plasmid via TAR cloning. "
            "Parts from all listed libraries (public and private) are merged automatically."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {"type": "string", "description": "Construct name."},
                "assembly_order": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Ordered part names, 5' to 3'.",
                },
                "library_names": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "One or more library folder names. Parts from all are available in the design.",
                },
            },
            "required": ["name", "assembly_order", "library_names"],
        },
    },
    {
        "name": "design_integration",
        "description": (
            "Insert a linear cassette at a chromosomal locus via homologous recombination. "
            "Parts from all listed libraries (public and private) are merged automatically."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {"type": "string", "description": "Construct name."},
                "assembly_order": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Part names to integrate, 5' to 3' (flanks not included).",
                },
                "library_names": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "One or more library folder names. Parts from all are available in the design.",
                },
                "integration_site_name": {
                    "type": "string",
                    "description": "Integration locus name (e.g. 'HO', 'HIS3', 'LEU2').",
                },
            },
            "required": ["name", "assembly_order", "library_names", "integration_site_name"],
        },
    },
    {
        "name": "design_replacement",
        "description": (
            "Design a pop-in/pop-out gene replacement cassette using URA3 counter-selection. "
            "Supply replacement via library (replacement_library_name + replacement_component_name) "
            "or raw DNA string (replacement_sequence)."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {"type": "string", "description": "Construct name."},
                "target_sequence": {"type": "string", "description": "Genomic DNA to replace."},
                "replacement_library_name": {"type": "string", "description": "Library for replacement part."},
                "replacement_component_name": {"type": "string", "description": "Part name in library."},
                "replacement_sequence": {"type": "string", "description": "Raw DNA to insert (A/T/G/C)."},
                "marker_position": {
                    "type": "string",
                    "enum": ["upstream", "downstream"],
                    "description": "URA3 marker position relative to replacement.",
                },
            },
            "required": ["name", "target_sequence", "marker_position"],
        },
    },
    {
        "name": "design_deletion",
        "description": "Design a scarless deletion cassette using URA3 marker recycling.",
        "input_schema": {
            "type": "object",
            "properties": {
                "target_sequence": {"type": "string", "description": "Genomic sequence to delete."},
                "name": {"type": "string", "description": "Construct name."},
            },
            "required": ["target_sequence", "name"],
        },
    },
    {
        "name": "run_batch",
        "description": (
            "Consolidate TAR/integration designs into optimized PCR batches with instructions. "
            "Does not apply to deletion/replacement cassettes."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "construct_names": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Constructs to batch; omit to batch all available.",
                },
                "liquid_handler": {
                    "type": "string",
                    "enum": ["Human", "epMotion", "Janus", "Hamilton"],
                    "description": "Target instrument (default: Human).",
                },
            },
        },
    },
]

# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------

_SGD_API = "https://www.yeastgenome.org/backend/locus/{gene}/sequence_details"
_SGD_KEY = {"orf": "coding_dna", "genomic": "genomic_dna", "upstream_1kb": "1kb"}


def _handle_lookup_gene_sequence(gene_name: str, sequence_type: str = "orf") -> dict:
    url = _SGD_API.format(gene=gene_name.upper())
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PYEAST-agent/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return {
                "success": False,
                "error": (
                    f"Gene '{gene_name}' not found in SGD. "
                    "Check spelling or try the systematic name (e.g. YPL154C)."
                ),
            }
        return {"success": False, "error": f"SGD API error {e.code}: {e.reason}"}
    except (urllib.error.URLError, TimeoutError) as e:
        return {
            "success": False,
            "error": f"Could not reach SGD (network error: {e}). Provide the sequence manually.",
        }

    api_key = _SGD_KEY.get(sequence_type, "coding_dna")
    entries = data.get(api_key, [])
    if not entries:
        return {"success": False, "error": f"No {sequence_type} sequence found for '{gene_name}'."}

    # Prefer S288C / Reference strain; fall back to first entry
    preferred = next(
        (e for e in entries if "S288C" in e.get("strain", {}).get("display_name", "")
         or "Reference" in e.get("strain", {}).get("display_name", "")),
        entries[0],
    )
    sequence = preferred.get("residues", "")
    if not sequence:
        return {"success": False, "error": f"Empty sequence returned for '{gene_name}'."}

    locus = preferred.get("locus", {})
    return {
        "success": True,
        "gene_name": gene_name.upper(),
        "systematic_name": locus.get("format_name", locus.get("systematic_name", "")),
        "sequence_type": sequence_type,
        "sequence": sequence,
        "length_bp": len(sequence),
        "strand": preferred.get("strand", ""),
        "description": locus.get("headline", locus.get("description", "")),
        "source": "SGD S288C reference",
        "sgd_url": f"https://www.yeastgenome.org/locus/{gene_name.upper()}",
        "note": "Sequence is from S288C, compatible with BY4741 background.",
    }


def _handle_list_components(component_type: str, library_name: str | None = None) -> dict:
    pub_base = get_component_libraries_path()
    priv_base = get_component_libraries_path(private=True)

    if component_type == "libraries":
        if not pub_base.exists():
            return {
                "error": f"Component libraries not found at {pub_base}. "
                "Run 'pyeast init' to configure your data directory."
            }
        pub_libs = {d.name for d in pub_base.iterdir() if d.is_dir()}
        priv_libs = {d.name for d in priv_base.iterdir() if d.is_dir()} if priv_base.exists() else set()
        all_libs = sorted(pub_libs | priv_libs)
        return {
            "libraries": all_libs,
            "private_libraries": sorted(priv_libs),
            "count": len(all_libs),
        }

    if component_type == "components":
        if not library_name:
            return {"error": "library_name is required when component_type='components'"}
        parts: set[str] = set()
        found_any = False
        for base in [pub_base, priv_base]:
            lib = base / library_name
            if lib.exists():
                found_any = True
                parts |= {f.stem for f in lib.iterdir() if f.is_file() and f.suffix.lower() in _SEQ_EXTENSIONS}
        if not found_any:
            all_libs = sorted({d.name for d in pub_base.iterdir() if d.is_dir()} |
                              ({d.name for d in priv_base.iterdir() if d.is_dir()} if priv_base.exists() else set()))
            return {"error": f"Library '{library_name}' not found.", "available_libraries": all_libs}
        return {"library": library_name, "components": sorted(parts), "count": len(parts)}

    if component_type == "integration_sites":
        path = get_integration_sites_path()
        if not path.exists():
            return {"error": f"Integration sites directory not found: {path}"}
        sites = sorted([
            f.stem for f in path.iterdir()
            if f.is_file() and f.suffix.lower() in _SEQ_EXTENSIONS
        ])
        return {"integration_sites": sites, "count": len(sites)}

    if component_type == "templates":
        path = get_templates_path()
        if not path.exists():
            return {"error": f"Templates directory not found: {path}"}
        templates = sorted([f.name for f in path.iterdir() if f.is_file()])
        return {"templates": templates, "count": len(templates)}

    return {"error": f"Unknown component_type: '{component_type}'"}


def _handle_list_outputs() -> dict:
    output_dir = get_output_path()
    if not output_dir.exists():
        return {"constructs": [], "count": 0, "note": "Output directory does not exist yet."}

    constructs = []
    # Subfolders (new organized structure)
    for subdir in sorted(output_dir.iterdir()):
        if not subdir.is_dir():
            continue
        for gb in sorted(subdir.glob("*.gb")):
            constructs.append({"name": gb.stem, "folder": subdir.name})
    # Root level (legacy)
    for gb in sorted(output_dir.glob("*.gb")):
        constructs.append({"name": gb.stem, "folder": "output/"})

    return {"constructs": constructs, "count": len(constructs)}


def _handle_read_component(name: str, library_name: str) -> dict:
    from Bio import SeqIO

    pub_base = get_component_libraries_path()
    priv_base = get_component_libraries_path(private=True)

    for base in [priv_base, pub_base]:  # private takes precedence
        for ext in (".fasta", ".fa", ".gb", ".genbank", ".gbk"):
            file_path = base / library_name / f"{name}{ext}"
            if file_path.exists():
                fmt = "genbank" if ext in (".gb", ".genbank", ".gbk") else "fasta"
                try:
                    record = SeqIO.read(str(file_path), fmt)
                    return {
                        "success": True,
                        "name": name,
                        "library": library_name,
                        "sequence": str(record.seq),
                        "length_bp": len(record.seq),
                        "description": record.description,
                    }
                except Exception as e:
                    return {"success": False, "error": f"Failed to parse {file_path.name}: {e}"}

    return {
        "success": False,
        "error": f"Part '{name}' not found in library '{library_name}' (public or private). "
        "Use list_components to see available parts.",
    }


def _handle_save_component(
    name: str,
    sequence: str,
    library_name: str,
    description: str = "",
    overwrite: bool = False,
) -> dict:
    # Parse FASTA if a header line is present
    lines = sequence.strip().splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    cleaned = "".join(lines).strip().upper().replace(" ", "").replace("\t", "")

    if not cleaned:
        return {"success": False, "error": "Sequence is empty after stripping whitespace."}

    invalid = set(cleaned) - set("ACGTRYSWKMBDHVN")
    if invalid:
        return {
            "success": False,
            "error": f"Sequence contains non-DNA characters: {sorted(invalid)}. Expected A/T/G/C.",
        }

    library_path = get_component_libraries_path() / library_name
    library_path.mkdir(parents=True, exist_ok=True)

    file_path = library_path / f"{name}.fasta"
    if file_path.exists() and not overwrite:
        return {
            "success": False,
            "error": (
                f"A component named '{name}' already exists in library '{library_name}'. "
                "Pass overwrite=true to replace it (e.g. to update the description)."
            ),
        }

    header = f">{name} {description}".strip()
    file_path.write_text(f"{header}\n{cleaned}\n")

    part_count = sum(
        1 for f in library_path.iterdir()
        if f.is_file() and f.suffix.lower() in _SEQ_EXTENSIONS
    )

    return {
        "success": True,
        "name": name,
        "library": library_name,
        "length_bp": len(cleaned),
        "file_path": str(file_path),
        "library_part_count": part_count,
        "note": (
            f"'{name}' is now available for use in design_tar / design_integration "
            f"with library_name='{library_name}'."
        ),
    }


def _merge_libraries_to_tempdir(library_names: list[str]) -> Path:
    """Copy parts from named libraries (public + private) into a fresh temp directory.
    Private files overwrite public files with the same name. Caller must rmtree when done."""
    tmpdir = Path(tempfile.mkdtemp())
    pub_base = get_component_libraries_path()
    priv_base = get_component_libraries_path(private=True)
    for lib_name in library_names:
        for base in [pub_base, priv_base]:
            lib_path = base / lib_name
            if lib_path.exists():
                for f in lib_path.iterdir():
                    if f.is_file() and f.suffix.lower() in _SEQ_EXTENSIONS:
                        shutil.copy2(f, tmpdir / f.name)
    return tmpdir


def _find_missing_parts(library_names: list[str], assembly_order: list[str]) -> list[str]:
    """Return part names from assembly_order not found in any of the named libraries (public or private)."""
    pub_base = get_component_libraries_path()
    priv_base = get_component_libraries_path(private=True)
    available: set[str] = set()
    for lib_name in library_names:
        for base in [pub_base, priv_base]:
            lib_path = base / lib_name
            if lib_path.exists():
                available |= {f.stem for f in lib_path.iterdir() if f.is_file() and f.suffix.lower() in _SEQ_EXTENSIONS}
    return [p for p in assembly_order if p not in available]


def _handle_design_tar(name: str, assembly_order: list[str], library_names: list[str]) -> dict:
    from pyeast.core.tar import TARDesigner

    missing = _find_missing_parts(library_names, assembly_order)
    if missing:
        return {
            "success": False,
            "missing_parts": missing,
            "error": (
                f"Cannot design: {len(missing)} part(s) not found in libraries {library_names}: "
                f"{missing}. "
                "Ask the user to provide the DNA sequence(s) for these parts before proceeding."
            ),
        }

    output_dir = ensure_output_dir_exists(name)
    output_prefix = output_dir / name

    tmpdir = _merge_libraries_to_tempdir(library_names)
    try:
        designer = TARDesigner()
        result = designer.design(
            library_path=tmpdir,
            assembly_order=assembly_order,
            primer_folder=get_primers_path(),
            template_folder=get_templates_path(),
            name=name,
        )
        result.save(output_prefix)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    output_files = [
        str(output_prefix) + ".gb",
        str(output_prefix) + "_all_primers.tsv",
        str(output_prefix) + "_instructions.tsv",
    ]
    if result.missing_primers:
        output_files.append(str(output_prefix) + "_missing_primers.tsv")

    return {
        "success": True,
        "type": "TAR",
        "construct_name": name,
        "assembly_length_bp": len(result.assembly.seq),
        "parts_assembled": len(assembly_order),
        "primers_designed": len(result.primers),
        "missing_primers": len(result.missing_primers),
        "output_files": output_files,
        "primer_note": (
            f"{len(result.missing_primers)} primers need ordering — see _missing_primers.tsv."
            if result.missing_primers
            else "All primers found in plate inventory."
        ),
        "batch_eligible": True,
    }


def _handle_design_integration(
    name: str,
    assembly_order: list[str],
    library_names: list[str],
    integration_site_name: str,
) -> dict:
    from pyeast.core.integration import IntegrationDesigner

    missing = _find_missing_parts(library_names, assembly_order)
    if missing:
        return {
            "success": False,
            "missing_parts": missing,
            "error": (
                f"Cannot design: {len(missing)} part(s) not found in libraries {library_names}: "
                f"{missing}. "
                "Ask the user to provide the DNA sequence(s) for these parts before proceeding."
            ),
        }

    output_dir = ensure_output_dir_exists(name)
    output_prefix = output_dir / name

    tmpdir = _merge_libraries_to_tempdir(library_names)
    try:
        designer = IntegrationDesigner()
        result = designer.design(
            components_dir=tmpdir,
            assembly_order=assembly_order,
            integration_site_name=integration_site_name,
            primer_folder=get_primers_path(),
            template_folder=get_templates_path(),
            name=name,
        )
        result.save(output_prefix)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    output_files = [
        str(output_prefix) + ".gb",
        str(output_prefix) + "_all_primers.tsv",
        str(output_prefix) + "_instructions.tsv",
    ]
    if result.missing_primers:
        output_files.append(str(output_prefix) + "_missing_primers.tsv")

    return {
        "success": True,
        "type": "Integration",
        "construct_name": name,
        "integration_site": integration_site_name,
        "assembly_length_bp": len(result.assembly.seq),
        "parts_integrated": len(assembly_order),
        "primers_designed": len(result.primers),
        "missing_primers": len(result.missing_primers),
        "output_files": output_files,
        "primer_note": (
            f"{len(result.missing_primers)} primers need ordering — see _missing_primers.tsv."
            if result.missing_primers
            else "All primers found in plate inventory."
        ),
        "batch_eligible": True,
    }


def _handle_design_replacement(
    name: str,
    target_sequence: str,
    marker_position: str,
    replacement_library_name: str | None = None,
    replacement_component_name: str | None = None,
    replacement_sequence: str | None = None,
) -> dict:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from pyeast.core.replace import ReplaceDesigner

    # Resolve the replacement SeqRecord
    if replacement_library_name and replacement_component_name:
        lib_path = get_component_libraries_path() / replacement_library_name
        if not lib_path.exists():
            return {"success": False, "error": f"Library '{replacement_library_name}' not found"}
        matches = [
            f for f in lib_path.iterdir()
            if f.is_file() and f.stem == replacement_component_name
            and f.suffix.lower() in _SEQ_EXTENSIONS
        ]
        if not matches:
            return {
                "success": False,
                "error": f"Component '{replacement_component_name}' not found in library '{replacement_library_name}'",
            }
        from Bio import SeqIO
        replacement_record = SeqIO.read(str(matches[0]), matches[0].suffix.lstrip(".") or "fasta")
    elif replacement_sequence:
        replacement_record = SeqRecord(
            Seq(replacement_sequence.upper()),
            id=name,
            name=name,
            description="replacement sequence",
        )
    else:
        return {
            "success": False,
            "error": "Provide either (replacement_library_name + replacement_component_name) or replacement_sequence.",
        }

    output_dir = ensure_output_dir_exists(name)
    output_prefix = output_dir / name

    designer = ReplaceDesigner()
    result = designer.design(
        target_sequence=target_sequence,
        replacement_sequence=replacement_record,
        marker_position=marker_position,
        name=name,
    )
    result.save(output_prefix)

    return {
        "success": True,
        "type": "Replacement",
        "construct_name": name,
        "cassette_length_bp": len(result.cassette.seq),
        "genome_location": {
            "chromosome": result.genome_location[0],
            "start": result.genome_location[1],
            "end": result.genome_location[2],
            "orientation": result.genome_location[3],
        },
        "marker_position": marker_position,
        "screening_primers": {
            "forward": result.forward_primer,
            "reverse": result.reverse_primer,
        },
        "expected_pcr_products_bp": result.product_sizes,
        "output_files": [
            str(output_prefix) + ".gb",
            str(output_prefix) + ".fasta",
            str(output_prefix) + "_screening_primers.tsv",
        ],
        "batch_eligible": False,
        "note": "Replacement cassettes are not batched. Use the screening primers TSV for colony PCR.",
    }


def _handle_design_deletion(target_sequence: str, name: str) -> dict:
    from pyeast.core.deletion import DeletionDesigner

    output_dir = ensure_output_dir_exists(name)
    output_prefix = output_dir / name

    designer = DeletionDesigner()
    result = designer.design(target_sequence=target_sequence, name=name)
    result.save(output_prefix)

    return {
        "success": True,
        "type": "Deletion",
        "construct_name": name,
        "cassette_length_bp": len(result.cassette.seq),
        "genome_location": {
            "chromosome": result.genome_location[0],
            "start": result.genome_location[1],
            "end": result.genome_location[2],
            "orientation": result.genome_location[3],
        },
        "screening_primers": {
            "forward": result.forward_primer,
            "reverse": result.reverse_primer,
        },
        "expected_pcr_products_bp": result.product_sizes,
        "output_files": [
            str(output_prefix) + ".gb",
            str(output_prefix) + ".fasta",
            str(output_prefix) + "_screening_primers.tsv",
        ],
        "batch_eligible": False,
        "note": "Deletion cassettes are not batched. Use the screening primers TSV for colony PCR.",
    }


def _handle_run_batch(
    construct_names: list[str] | None = None,
    liquid_handler: str = "Human",
) -> dict:
    from pyeast.core.batch import BatchDesigner

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = get_output_path()
    output_prefix = str(output_dir / f"batch_{timestamp}")

    designer = BatchDesigner()
    designer.load_constructs()

    if not designer.available_constructs:
        return {
            "success": False,
            "error": "No TAR or integration constructs found in output directory. "
            "Design some constructs first using design_tar or design_integration.",
        }

    # Select constructs
    if construct_names:
        missing = [n for n in construct_names if n not in designer.available_constructs]
        if missing:
            return {
                "success": False,
                "error": f"Constructs not found: {missing}. "
                f"Available: {list(designer.available_constructs.keys())}",
            }
        designer.selected_constructs = {n: designer.available_constructs[n] for n in construct_names}
    else:
        designer.selected_constructs = dict(designer.available_constructs)

    designer.validate_constructs()
    designer.process_selected_constructs()
    designer.find_primers_and_templates()
    designer.organize_pcr_batches()

    designer.generate_human_instructions(output_prefix)
    designer.save_input_record(output_prefix)

    output_files = [
        f"{output_prefix}_instructions.txt",
        f"{output_prefix}_inputs.txt",
    ]

    if liquid_handler != "Human":
        designer.generate_assembly_groups(output_prefix)
        machine_file = designer.generate_machine_assembly_instructions(
            output_prefix, liquid_handler, timestamp
        )
        if machine_file:
            output_files.append(machine_file)

    return {
        "success": True,
        "constructs_batched": list(designer.selected_constructs.keys()),
        "total_constructs": len(designer.selected_constructs),
        "liquid_handler": liquid_handler,
        "output_files": output_files,
        "output_prefix": output_prefix,
    }


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def dispatch_tool(tool_name: str, tool_input: dict) -> str:
    """Execute a tool call and return a JSON string for Claude."""
    try:
        if tool_name == "read_component":
            result = _handle_read_component(
                name=tool_input["name"],
                library_name=tool_input["library_name"],
            )
        elif tool_name == "save_component":
            result = _handle_save_component(
                name=tool_input["name"],
                sequence=tool_input["sequence"],
                library_name=tool_input["library_name"],
                description=tool_input.get("description", ""),
                overwrite=tool_input.get("overwrite", False),
            )
        elif tool_name == "lookup_gene_sequence":
            result = _handle_lookup_gene_sequence(
                gene_name=tool_input["gene_name"],
                sequence_type=tool_input.get("sequence_type", "orf"),
            )
        elif tool_name == "list_components":
            result = _handle_list_components(
                component_type=tool_input["component_type"],
                library_name=tool_input.get("library_name"),
            )
        elif tool_name == "list_outputs":
            result = _handle_list_outputs()
        elif tool_name == "design_tar":
            result = _handle_design_tar(
                name=tool_input["name"],
                assembly_order=tool_input["assembly_order"],
                library_names=tool_input["library_names"],
            )
        elif tool_name == "design_integration":
            result = _handle_design_integration(
                name=tool_input["name"],
                assembly_order=tool_input["assembly_order"],
                library_names=tool_input["library_names"],
                integration_site_name=tool_input["integration_site_name"],
            )
        elif tool_name == "design_replacement":
            result = _handle_design_replacement(
                name=tool_input["name"],
                target_sequence=tool_input["target_sequence"],
                marker_position=tool_input["marker_position"],
                replacement_library_name=tool_input.get("replacement_library_name"),
                replacement_component_name=tool_input.get("replacement_component_name"),
                replacement_sequence=tool_input.get("replacement_sequence"),
            )
        elif tool_name == "design_deletion":
            result = _handle_design_deletion(
                target_sequence=tool_input["target_sequence"],
                name=tool_input["name"],
            )
        elif tool_name == "run_batch":
            result = _handle_run_batch(
                construct_names=tool_input.get("construct_names") or None,
                liquid_handler=tool_input.get("liquid_handler", "Human"),
            )
        else:
            result = {"error": f"Unknown tool: '{tool_name}'"}
    except Exception as e:
        result = {"error": str(e), "tool": tool_name}

    return json.dumps(result, default=str)
