[![DOI](https://zenodo.org/badge/946400824.svg)](https://doi.org/10.5281/zenodo.15393309)
[![Installation Test](https://github.com/TomLoan/PYEAST/actions/workflows/install-test.yml/badge.svg)](https://github.com/TomLoan/PYEAST/actions/workflows/install-test.yml)

# Set of tools for yeast cloning, built on the CSIRO uv-python template. 

## Quick start:

### Prerequisites

UV, Git

#### Install uv

To install uv, see the instructions: https://docs.astral.sh/uv/getting-started/installation/, but in short:

In bash on (most) Linux systems and Mac:
```bash
curl -LsSf https://astral.sh/uv/0.4.6/install.sh | sh
# restart your shell and make sure `uv --version` works
```

In Windows, from any shell run:
```shell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex" # windows
# restart your shell and make sure `uv --version` works
```

# Clone the repo
```shel
git clone https://github.com/TomLoan/PYEAST.git
```

### Using PYEAST 

Run PYEAST in Command line
```shell
uv run pyeast
```
uv will handle the package management and create athe required virtual environment in the local directory. 
Once this process is complete a range of commands will be printed to the terminal, use uv run pyeast command --help for more information on running each command

Pyeast provides  useful functions for  your Saccharomyces cerevisiae cloning needs: 
tar: Transformation assisted recombination is a cloning method that relies on S. cerevisiaes native capcity for homologous recombination to create new plasmids out of PCR products using homology added to the ends of PCR primers. 
insert: Similar to tar this script designs PCR primers for insertion of new DNA into the chromosomes of S. cerevisiae 
del: This script designs DNA fragments that can be used to delete regions of S. cerevisiae gDNA using the scarless method details by Akada et al 2006 (Yeast 23(5):399-405). 
replace: Similar to del this script designs DNA fragments that can be used to replace regions of S. cerevisiae gDNA with a method based on that described by Akada et al 2006 (Yeast 23(5):399-405)
batch: regenate instructions files for previously designed tar and integrate command outputs stored in the output file. You might need to do this for example after ordering primers or adding templates sequences to ./data/templates. 

For more details see our [pre-print](https://doi.org/10.1101/2025.05.19.655004) on BioRxiv 


## Version 1.1 Adds support for other liquid handlering robots and for MoClo Golden gate assemblies. 

The following Kits come preloaded: 
- Yeast ToolKit: [Addgene](https://www.addgene.org/kits/moclo-ytk/) from [Lee at al, 2015](http://doi.org/10.1021/sb500366v)

- Yeast Secretion and Display Toolkit: [Addgene](https://www.addgene.org/kits/young-moclo-ysd/) from [O'Riordan et al, 2023](https://doi.org/10.1021/acssynbio.3c00743)

- OPENPichia MoClo Kit: [BCCM](https://bccm.belspo.be/GeneCorner-OPENPichia) from [Claes et al, 2024](https://www.nature.com/articles/s41564-023-01574-w)

New components can be added to any oy these kits by saving a fasta file to .data/component libraries/*kit_name* and a gb file of the MoClo complatable lvl 0 plasmid to .data/components libraries/*kit_name*/plasmids. For liquid handling the plasmid name should be added to a well in .data/templates/TemPlates.xlsx.

To add new golden gate kits save fasta files describing the parts to a new folder in data/component libraries and the lvl 0 plasmids saved in a plasmids subfolder in the same directory. For liquid handling it is often simpliest to save a new 96 plate in .data/templates/TemPlates.xlsx as a new sheet. Be sure the positions of these plasmids are accurate. 

Private Data Support (New in v1.1)
PYEAST now supports private data directories for proprietary sequences, primers, and templates. Create data/private/component_libraries/ and add your files - PYEAST automatically searches both public and private locations. See data/README.md for details.
