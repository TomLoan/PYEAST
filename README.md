[![DOI](https://zenodo.org/badge/946400824.svg)](https://doi.org/10.5281/zenodo.15393309)
[![Installation Test](https://github.com/TomLoan/PYEAST/actions/workflows/install-test.yml/badge.svg)](https://github.com/TomLoan/PYEAST/actions/workflows/install-test.yml)

# PYEAST

PYEAST is a command-line toolkit and Python library that automates the design of DNA cloning experiments in *Saccharomyces cerevisiae* (baker's yeast). Given a set of genetic parts, it designs the PCR primers and liquid-handling instructions for common yeast genetic engineering techniques — reducing manual design work and minimising errors.

If you are new to yeast cloning, PYEAST is designed to complement standard wet-lab protocols: it handles the computational design steps so you can focus on the biology.

For full methodological details, see our [publication](https://doi.org/10.1038/s41540-026-00712-4).

## Prerequisites

- **Python** (3.12 or later) — [python.org](https://www.python.org/downloads/)
- **Git** — [git-scm.com](https://git-scm.com/downloads)
- **uv** - [docs.astral.sh](https://docs.astral.sh/uv/getting-started/installation/)

## Installation

```bash
uv tool install git+https://github.com/TomLoan/PYEAST.git

# Download the sequence data automatically
pyeast init
```

`pyeast init` clones the [data repository](https://github.com/TomLoan/PYEAST_data) to `~/.pyeast/data-repo/` and configures PYEAST automatically. If you have already downloaded the data elsewhere, point PYEAST at it instead:

```bash
pyeast init --data-dir /path/to/PYEAST_data
```

## Commands

Run `pyeast --help` to see all available commands:

| Command | Description |
|---|---|
| `pyeast tar` | Design primers for Transformation-Assisted Recombination (TAR) cloning |
| `pyeast integrate` | Design primers for chromosomal integration |
| `pyeast delete` | Design deletion cassettes (scarless marker-recycling method) |
| `pyeast replace` | Design replacement cassettes (scarless marker-recycling method)|
| `pyeast gg` | Design Golden Gate / MoClo assemblies |
| `pyeast batch` | Regenerate instruction files for previously designed assemblies |
| `pyeast agent` | Interactive LLM-powered experiment design session |
| `pyeast init` | Configure the data directory |

For help with any command: `pyeast COMMAND --help`

## TAR Cloning

Transformation-Assisted Recombination (TAR) exploits *S. cerevisiae*'s natural homologous recombination to assemble PCR products into new plasmids in a single yeast transformation. PYEAST designs all primers, including the homology overhangs, for each part in your assembly.

## Chromosomal Integration

Design primers for inserting a DNA construct into a specific locus in the yeast genome. PYEAST calculates the flanking homology sequences required for efficient integration.

## Gene Deletion and Replacement

Design cassettes for deleting or replacing genes using the scarless marker-recycling method described by [Akada et al., 2006](https://doi.org/10.1002/yea.1365).

## Beta feature: Golden Gate / MoClo Assembly

== Beta feature, please double check outputs ==

Design Golden Gate assemblies using the MoClo standard. The following part libraries are included out of the box:

| Kit | Reference |
|---|---|
| Yeast ToolKit (YTK) | [Lee et al., 2015](http://doi.org/10.1021/sb500366v) — [Addgene](https://www.addgene.org/kits/moclo-ytk/) |
| Yeast Secretion and Display Toolkit | [O'Riordan et al., 2023](https://doi.org/10.1021/acssynbio.3c00743) — [Addgene](https://www.addgene.org/kits/young-moclo-ysd/) |
| OPENPichia MoClo Kit | [Claes et al., 2024](https://www.nature.com/articles/s41564-023-01574-w) — [BCCM](https://bccm.belspo.be/GeneCorner-OPENPichia) |

To add parts to an existing kit, save a FASTA file to `component_libraries/<kit_name>/` in your data directory. For liquid-handling support, also save the matching plasmid `.gb` file to `component_libraries/<kit_name>/plasmids/` and add a well entry to `templates/TemPlates.xlsx`.

To add a new Golden Gate kit, create a new subfolder in `component_libraries/` and follow the same structure.

## LLM Agent

`pyeast agent` starts an interactive session where you describe your experiment in plain language. The agent discovers available components, looks up gene sequences from SGD, designs TAR cloning, chromosomal integrations, gene deletions and replacements, and generates consolidated PCR batch instructions — calling the same underlying functions as the CLI commands.

> **Note:** The agent is an interface layer over the standard design tools. Double-check outputs against direct CLI command results if you have any doubts.

Three LLM backends are supported. Copy `.env.example` to `.env` and fill in the relevant keys.

**Anthropic (default)** — requires `ANTHROPIC_API_KEY`.
```bash
pyeast agent
```

**Ollama** — requires a locally running [Ollama](https://ollama.com) server with a tool-calling model.
```bash
ollama serve && ollama pull qwen2.5
pyeast agent --provider ollama --model qwen2.5
```

**OpenAI-compatible servers (e.g. LM Studio)** — requires a tool-calling model loaded and the server started in LM Studio's Developer tab.
```bash
pyeast agent --provider openai --model <model-key>
# or, if your server is not on the default port (1234):
pyeast agent --provider openai --model <model-key> --base-url http://localhost:5678/v1
```

For full setup instructions see [docs/](docs/).

## Private Data

You can store proprietary sequences, primers, and templates in a `private/` subdirectory of your data folder. PYEAST searches both public and private locations automatically, and the private folder is excluded from version control.

```
<data directory>/
├── component_libraries/       ← public parts
├── integration_sites/
├── primers/
├── templates/
└── private/
    ├── component_libraries/   ← your private parts (gitignored)
    ├── integration_sites/
    ├── primers/
    └── templates/
```

When you run `pyeast init`, the `private/` folder structure is created automatically.

## Configuration

PYEAST looks for data in this order:

1. Environment variable: `PYEAST_DATA_DIR`
2. Config file: `~/.pyeast/config.yaml`
3. Default: `~/PYEAST/data/`

To view or update your current configuration, run `pyeast init`.

## Troubleshooting

**"Data directory not found"**

Run `pyeast init` to check what path is configured. To reconfigure:

```bash
pyeast init --data-dir /correct/path/to/PYEAST_data
```

**Commands fail with missing files**

Verify that your data directory contains the expected folder structure and that FASTA/Excel files are present in the right locations.

**Import errors on startup**

Reinstall to ensure all dependencies are present:

```bash
pip install git+https://github.com/TomLoan/PYEAST.git
```

## Python API

All designer classes can be imported and used directly — no interactive prompts, no CLI required. Useful for scripting, Jupyter notebooks, or integrating PYEAST into larger pipelines.

Each designer follows the same pattern: create a designer, call `design()` with your inputs, get back a result object, and call `result.save()` to write files.

```python
from pathlib import Path
from pyeast.core.tar import TARDesigner
from pyeast.utils.path_utils import get_component_libraries_path, get_primers_path, get_templates_path

designer = TARDesigner(homology_length=25)
result = designer.design(
    library_path=get_component_libraries_path() / "YeastToolKit",
    assembly_order=["pTEF1", "YeGFP", "tCYC1", "Ura3", "AmpR_ColE1", "2Micron"],
    primer_folder=get_primers_path(),
    template_folder=get_templates_path(),
    name="my_plasmid",
)
result.save(Path("output/my_plasmid/my_plasmid"))
```

The same `designer.design() → result.save()` pattern applies to `IntegrationDesigner`, `DeletionDesigner`, and `ReplaceDesigner`. See [docs/](docs/) or run `help(ClassName)` for constructor options and full parameter details.

## For Developers

```bash
git clone https://github.com/TomLoan/PYEAST.git
cd PYEAST
uv sync --group dev
uv run pytest tests/ -v
```

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## Citation

If you use PYEAST in your research, please cite:

> Madika, A., Suri, A., Purohit, A. et al. PYEAST – A Computational Toolkit for Saccharomyces cerevisiae Genetic Engineering. npj Syst Biol Appl 12, 83 (2026). https://doi.org/10.1038/s41540-026-00712-4
