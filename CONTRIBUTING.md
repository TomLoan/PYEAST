# Contributing to PYEAST

Thanks for your interest in contributing to PYEAST! This project is developed primarily for the yeast synthetic biology community, and contributions that help make it more useful are very welcome.

## Two Repositories

PYEAST is split across two repositories:

- **[PYEAST](https://github.com/TomLoan/PYEAST)** (this repo) — the Python package, CLI, and tests
- **[PYEAST_data](https://github.com/TomLoan/PYEAST_data)** — the sequence data: component libraries, integration sites, primers, and templates

If you want to add or fix **sequences** (promoters, terminators, integration sites, etc.), contribute to **PYEAST_data**. If you want to fix a bug or add a feature to the tool itself, contribute to this repo.

## Ways to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub with:
- A clear description of the problem
- Steps to reproduce it
- Your operating system and Python version

### Suggesting Features

Feature requests are welcome! Please open an issue describing:
- What you'd like PYEAST to do
- Why this would be useful for your workflow
- Any ideas on how it might work

### Submitting Code

1. **Fork the repository** and create a branch from `develop`
2. **Make your changes** - try to keep them focused on a single issue
3. **Test locally** - at minimum, verify `uv run pyeast --help` and any commands you've modified
4. **Submit a pull request** to the `develop` branch

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR-USERNAME/PYEAST.git
cd PYEAST

# Install with dev dependencies
uv sync --group dev

# Download the sequence data
uv run pyeast init

# Verify it works
uv run pyeast --help
```

`pyeast init` clones the PYEAST_data repository to `~/.pyeast/data-repo/`. If you already have the data elsewhere, point PYEAST at it instead:

```bash
uv run pyeast init --data-dir /path/to/PYEAST_data
```

Don't worry too much about perfect style - things can be tidied up during review.

## Adding Sequences

To add new sequences, contribute to the [PYEAST_data](https://github.com/TomLoan/PYEAST_data) repository:

- **Promoters, terminators, genes**: Add `.fasta` files to `component_libraries/Saccharomyces_cerevisiae/`, or create a new folder in `component_libraries/` for a different organism or kit.
- **Integration sites**: Add upstream/downstream flanking sequences in a single `.fasta` file to `integration_sites/`.
- **Templates**: Add `.gb` or `.fasta` files to `templates/`.

If you have a collection that would be useful to others, consider submitting a PR to PYEAST_data!

## Questions?

Feel free to open an issue if you have questions about contributing. Happy to help newcomers get started.
