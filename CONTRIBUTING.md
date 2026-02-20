# Contributing to PYEAST

Thanks for your interest in contributing to PYEAST! This project is developed primarily for the yeast synthetic biology community, and contributions that help make it more useful are very welcome.

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
2. **Make your changes** – try to keep them focused on a single issue
3. **Test locally** – at minimum, verify `uv run pyeast --help` and any commands you've modified
4. **Submit a pull request** to the `develop` branch

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR-USERNAME/PYEAST.git
cd PYEAST

# Install with dev dependencies
uv sync --group dev

# Verify it works
uv run pyeast --help
```

Don't worry too much about perfect style – things can be tidied up during review.

## Adding Components

One of the easiest ways to contribute is adding new components to the default libraries:

- **Promoters, terminators, genes**: Add `.fasta` files to `component_libraries/Saccharomyces_cerevisiae/` in the [PYEAST_data](https://github.com/TomLoan/PYEAST_data) repository, or create a new folder in `component_libraries/` and populate it with `.fasta` files.
- **Integration sites**: Add upstream/downstream flanking sequences in a single `.fasta` file to `integration_sites/`.
- **Templates**: Add `.gb` or `.fasta` files to `templates/`.

If you have a collection that would be useful to others, consider submitting a PR!

## Questions?

Feel free to open an issue if you have questions about contributing. Happy to help newcomers get started.
