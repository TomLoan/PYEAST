# Installing PYEAST

PYEAST can be installed in several ways depending on your use case.

## Quick Start (Recommended for Most Users)

### Option 1: Development Installation

Clone the repository and install in editable mode:

```bash
# Clone the repository
git clone https://github.com/TomLoan/PYEAST.git
cd PYEAST

# Install with uv (recommended)
uv pip install -e .

# Or install with pip
pip install -e .
```

Data will automatically be found in `./data/` when running from the repository directory.

### Option 2: Install from GitHub

Install directly from GitHub without cloning:

```bash
pip install git+https://github.com/TomLoan/PYEAST.git
```

After installation, set up data directory:

```bash
# Clone repo separately for data
git clone https://github.com/TomLoan/PYEAST.git ~/PYEAST-data

# Configure PYEAST to use this data
pyeast init --data-dir ~/PYEAST-data/data
```

## Configuration

PYEAST looks for data in this priority order:

1. **Environment variable:** `PYEAST_DATA_DIR`
   ```bash
   export PYEAST_DATA_DIR=/path/to/data
   ```

2. **Config file:** `~/.pyeast/config.yaml`
   ```yaml
   data_dir: /path/to/data
   output_dir: /path/to/output  # optional
   ```

3. **Dev mode:** `./data/` (if running from git checkout)

4. **Default:** `~/PYEAST/data/`

### Using the Init Command

The `pyeast init` command helps configure data paths:

```bash
# Point to existing data directory
pyeast init --data-dir /path/to/PYEAST/data

# Check current configuration
pyeast init
```

## Data Directory Structure

PYEAST expects data to be organized as follows:

```
data/
├── component libraries/      # Genetic part libraries
│   ├── OpenPichia MoClo lvl1/
│   ├── Saccharomyces cerevisiae/
│   ├── Yeast Moclo lvl 0/
│   └── Yeast MoClo lvl 1/
├── integration sites/        # Integration target sequences
├── primers/                  # Primer plates (Excel format)
├── templates/                # Template plates and genomes
│   ├── TemPlates.xlsx
│   ├── genome_well_mapping.tsv
│   └── BY4741_Toronto_2012.fsa
└── private/                  # Optional private data
    ├── component libraries/
    ├── integration sites/
    ├── primers/
    └── templates/
```

## Private Data

PYEAST supports private genetic parts and data:

- Public data: `$PYEAST_DATA_DIR/component libraries/`
- Private data: `$PYEAST_DATA_DIR/private/component libraries/`

Both directories are searched automatically. Private entries override public ones with the same name.

To add private data, create a `private/` subdirectory in your data folder:

```bash
mkdir -p ~/my-pyeast-data/private/component\ libraries/MyPrivateParts
# Add your private FASTA files...
```

## Verification

Test your installation:

```bash
# Check help works
pyeast --help

# Verify configuration
pyeast init

# Test a command (will prompt for input)
pyeast tar
```

## Troubleshooting

### "Data directory not found" errors

1. Check your configuration:
   ```bash
   pyeast init
   ```

2. Verify data directory exists:
   ```bash
   ls $PYEAST_DATA_DIR
   # Or check config file
   cat ~/.pyeast/config.yaml
   ```

3. Reconfigure if needed:
   ```bash
   pyeast init --data-dir /correct/path/to/data
   ```

### Commands fail with path errors

- Ensure data directory contains the expected structure
- Check file permissions
- Verify FASTA/Excel files are present in expected locations

### Import errors

Make sure all dependencies are installed:

```bash
pip install -e .
# Or
uv pip install -e .
```

## Advanced: Multiple Data Locations

You can switch between different data sets using environment variables:

```bash
# Use production data
export PYEAST_DATA_DIR=~/PYEAST/data
pyeast tar

# Use test data
export PYEAST_DATA_DIR=~/PYEAST-test/data
pyeast tar
```

Or create different config files and symlink as needed:

```bash
# Create configs
cat > ~/.pyeast/config-prod.yaml <<EOF
data_dir: /path/to/prod/data
EOF

cat > ~/.pyeast/config-test.yaml <<EOF
data_dir: /path/to/test/data
EOF

# Switch configs
ln -sf ~/.pyeast/config-prod.yaml ~/.pyeast/config.yaml
```

## Development Setup

For development work:

```bash
# Clone and enter repository
git clone https://github.com/TomLoan/PYEAST.git
cd PYEAST

# Install in editable mode with dev dependencies
uv sync

# Install pre-commit hooks
uv run pre-commit install

# Run tests
uv run pytest tests/ -v

# Run with coverage
uv run pytest tests/ --cov=src/pyeast --cov-report=html
```

Data will automatically be found in `./data/` when running from the repository.

## Getting Help

- Documentation: [README.md](README.md)
- Issues: https://github.com/TomLoan/PYEAST/issues
- Command help: `pyeast --help` or `pyeast COMMAND --help`
