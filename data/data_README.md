# PYEAST Data Directory

This directory contains all the genetic components, primers, templates, and integration sites used by PYEAST.

## Directory Structure

PYEAST organizes data into **public** and **private** directories:

```
data/
├── component_libraries/     # Shared genetic parts (in Git)
├── integration_sites/       # Shared integration sites (in Git)
├── primers/                 # Shared primer plates (in Git)
├── templates/               # Shared template plates (in Git)
└── private/                 # YOUR private data (NOT in Git)
    ├── component_libraries/
    ├── integration_sites/
    ├── primers/
    └── templates/
```

## Public vs Private

### Public Data (Committed to Git)
- Standard genetic parts (promoters, terminators, markers, etc.)
- Common integration sites
- Lab-wide primer plates
- Lab-wide template plates
- **Everyone using PYEAST sees this data**

### Private Data (Local Only, NOT in Git)
- Proprietary sequences
- Unpublished constructs
- Client work
- Personal projects
- Your own primer plates
- Your own template plates
- **Only you see this data**

## How It Works

PYEAST automatically searches both locations and combines the results:

1. **Adding to public libraries**: Drop files in `private/component_libraries/YTK/` to add your own parts to the public YTK library
2. **Creating private libraries**: Create `private/component_libraries/MyProject/` for completely private work
3. **Private overrides public**: If you have a private sequence with the same name as a public one, the private version is used

**You'll never see any indication of whether data is public or private** - PYEAST just shows you everything you have access to.

## Data Location Configuration

As of version 1.1.1, PYEAST supports flexible data directory locations. You can install PYEAST separately from its data.

### Data Location Priority

PYEAST looks for data in this order:

1. **Environment variable**: `PYEAST_DATA_DIR=/path/to/data`
2. **Config file**: `~/.pyeast/config.yaml`
3. **Dev mode**: `./data/` (if running from git checkout)
4. **Default**: `~/PYEAST/data/`

### Using Data in Different Locations

If you installed PYEAST via `pip install` (not from the git repository), you can configure where data lives:

#### Option 1: Environment Variable
```bash
export PYEAST_DATA_DIR=/path/to/PYEAST/data
pyeast tar
```

#### Option 2: Config File
Create `~/.pyeast/config.yaml`:
```yaml
data_dir: /path/to/PYEAST/data
output_dir: /path/to/output  # optional
```

#### Option 3: Init Command
```bash
# Point PYEAST to existing data directory
pyeast init --data-dir /path/to/PYEAST/data
```

#### Option 4: Git Checkout (Dev Mode)
When running from a git checkout, PYEAST automatically uses `./data/`:
```bash
git clone https://github.com/TomLoan/PYEAST.git
cd PYEAST
pip install -e .
pyeast tar  # Automatically uses ./data/
```

This flexibility means you can:
- Install PYEAST once, use it with multiple data directories
- Keep code updated via pip while maintaining stable data
- Share a single PYEAST installation across a team with different data sets
- Separate proprietary data from the PYEAST codebase

## Getting Started

### Using PYEAST with Public Data Only
Nothing special needed - just run PYEAST commands normally:
```bash
uv run pyeast integrate
```

### Adding Private Data

#### 1. Create Private Directories

Simply use the file editor on your machine to make new olders and files. Or in a shell:

```bash
# For sequences/components
mkdir -p data/private/component_libraries/MyProject

# For primers
mkdir -p data/private/primers

# For templates
mkdir -p data/private/templates

# For integration sites
mkdir -p data/private/integration_sites
```

#### 2. Add Your Files
```bash
# Add sequences (FASTA format)
cp my_proprietary_promoter.fasta data/private/component_libraries/MyProject/

# Add primer plates (Excel format - IDT spec sheets)
cp MyPlate_IDT.xlsx data/private/primers/

# Add template plates (Excel format)
cp MyTemPlates.xlsx data/private/templates/
```

#### 3. Use PYEAST Normally
PYEAST will automatically find and use both public and private data.

## Common Workflows

### Workflow 1: Add Private Parts to a Public Kit
You're using YTK but have some proprietary promoters:

```bash
mkdir -p data/private/component_libraries/YTK
cp my_secret_promoters.fasta data/private/component_libraries/YTK/
```

When you select YTK in PYEAST, you'll see both the public YTK parts AND your private additions.

### Workflow 2: Private-Only Project
You're working on confidential client work:

```bash
mkdir -p data/private/component_libraries/ClientX
cp client_sequences.fasta data/private/component_libraries/ClientX/
```

"ClientX" will appear in the library menu alongside public libraries.

### Workflow 3: Private Primer Plates
You ordered primers that only you use:

```bash
cp Plate2024-11-27_IDT.xlsx data/private/primers/
```

PYEAST will find primers in both public and private plates and choose the best option.

### Workflow 4: Golden Gate with Private Parts
For Golden Gate assemblies, you need both FASTA sequences and GenBank plasmids:

```bash
mkdir -p data/private/component_libraries/MyMoClo/plasmids
cp *.fasta data/private/component_libraries/MyMoClo/
cp *.gb data/private/component_libraries/MyMoClo/plasmids/
```

## File Formats

### Component Libraries (Sequences)
- **Format**: FASTA (.fasta, .fa, .fsa)
- **Location**: `component_libraries/LibraryName/` or `private/component_libraries/LibraryName/`
- **Content**: DNA sequences of genetic parts
- **For Golden Gate**: Also need GenBank files in `plasmids/` subdirectory

### Integration Sites  
- **Format**: FASTA (.fasta)
- **Location**: `integration_sites/` or `private/integration_sites/`
- **Content**: Two sequences per file (upstream and downstream homology)
- **Naming**: File name becomes the integration site name

### Primers
- **Format**: Excel (.xlsx) - IDT spec sheet format
- **Location**: `primers/` or `private/primers/`
- **Columns**: 'Plate or Box ID', 'Position', 'Sequence Name', 'Sequence'
- **Content**: Primer sequences and their plate positions

### Templates
- **Format**: Excel (.xlsx) - TemPlates.xlsx format
- **Location**: `templates/` or `private/templates/`
- **Layout**: 96-well plate format (rows A-H, columns 1-12)
- **Content**: Template names in each well position

### Genome Mappings
- **Format**: TSV (.tsv)
- **Location**: `templates/genome_well_mapping.tsv` or `private/templates/genome_well_mapping.tsv`
- **Columns**: Plate, Genome_Name, Well_Position, Contig_Names
- **Content**: Maps chromosome/contig names to plate positions

## Privacy & Git

The `private/` directory is in `.gitignore`, which means:

- ✅ Your private data stays on your computer
- ✅ You can share PYEAST code and public data
- ✅ No risk of accidentally committing proprietary sequences
- ✅ Each person can have different private data

**Important**: Even though the private directory exists, Git will never track or commit anything inside it.

## Tips

### Organization
- Use descriptive library names: `ClientABC`, `Project2024`, `Proprietary_Strains`
- Mirror the public structure: if public has `YTK/`, private can have `YTK/` too
- Keep related files together in the same library

### Backups
Private data is NOT in Git, so:
- Back up your `data/private/` directory separately
- Consider using your lab's file server or cloud storage
- Don't rely on Git for backup of private data

### Sharing
To share PYEAST with public data only:
```bash
git clone <repository>
# Your private/ directory doesn't exist in the clone
# The tool works fine with just public data
```

To share some private data with a collaborator:
```bash
# On your machine
zip -r client_data.zip data/private/component_libraries/ClientProject/

# Send client_data.zip to collaborator, they unzip into their data/private/
```

## Troubleshooting

### "Directory not found" error
Make sure you created the directory structure correctly:
```bash
# Check if directories exist
ls data/private/component_libraries/
```

### PYEAST doesn't see my private data
1. Check the directory structure matches exactly: `data/private/component_libraries/LibraryName/`
2. Check file extensions are correct (.fasta, .xlsx, etc.)
3. Try with public data first to verify PYEAST is working

### Private data accidentally in Git
This shouldn't happen (private/ is gitignored), but if it does:
```bash
# Remove from Git but keep local files
git rm -r --cached data/private/
git commit -m "Remove private data"
```

### Golden Gate assembly can't find plasmids
For Golden Gate, make sure your private library has BOTH:
1. FASTA files: `data/private/component_libraries/MyKit/*.fasta`
2. Plasmid GenBank files: `data/private/component_libraries/MyKit/plasmids/*.gb`

## Questions?

For more details, see the PYEAST documentation or open an issue on GitHub.

Remember: Public data is shared, private data is yours alone. PYEAST handles both simultaineously.
