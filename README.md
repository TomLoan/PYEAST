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
git clone https://github.com/TomLoan/PYEAST.git

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

For more details see our BioRxiv 
Link
