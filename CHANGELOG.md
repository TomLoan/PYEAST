# Changelog

All notable changes to PYEAST are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2026-07-30 
pyeast gg now only allows folders containing a correctly configured collection of golden 
gate plasmids to be selected.  

## [2.0.0] - 2026-07-07

Version 2.0.0 is a major release. The assembly engine, the packaged data, and the
supported Python version have all changed - see **Breaking changes** below for what to
do when upgrading from 1.x.

### Breaking changes

- **Sequence data now lives in a separate repository.** The `data/` directory has been
  removed from the code repo; sequences (component libraries, integration sites,
  primers, templates) are maintained in
  [PYEAST_data](https://github.com/TomLoan/PYEAST_data). Run `pyeast init` after
  installing to clone the data and write the config, or point PYEAST at an existing copy
  with `pyeast init --data-dir /path/to/PYEAST_data`.
- **Python 3.12 or later is now required** (previously 3.10+).
- **Component annotations use a standard feature key.** Assembled constructs now mark
  each part as a GenBank `misc_feature` carrying a `/note="PYEAST_component"` qualifier,
  instead of the custom `PYEAST_component` feature type. This survives round-trips
  through GenBank tools (SnapGene, Geneious, OpenCloning) that rewrite non-standard
  feature keys. Files written by earlier versions (custom type) are still read
  correctly.

### Added

- **pydna-based assembly engine for TAR and integration.** The assembly step now runs a
  real homologous-recombination simulation (pydna `in_vivo_assembly`) rather than a
  blind concatenation, including:
  - a pre-assembly ambiguity diagnostic that flags parts with repeated/shared homology
    before enumeration,
  - a PCR-specificity analysis that judges candidate priming sites by the melting
    temperature of their exact 3' footprint,
  - graceful degradation to a deterministically constructed product when the
    recombination graph is too dense to enumerate.
- **OpenCloning cloning-history export.** Each TAR/integration design writes a
  `*_history.json` (OpenCloning `CloningStrategy`) capturing the templates, primers, PCR
  products, and how they assemble - ready to upload to OpenCloning. Component parts are
  included in the history.
- **Complete Python library API.** The designers are decoupled from the CLI and can be
  driven programmatically: `TARDesigner`/`IntegrationDesigner` expose a `design(...)`
  method returning a result object (`TARResult`, `IntegrationResult`) with a
  `.save(output_prefix)` for writing GenBank, instruction, primer, map, and history
  files.
- **LLM-assisted design agent (`pyeast agent`).** An interactive session that helps
  design experiments, with support for Anthropic (default), OpenAI, and Ollama (local)
  providers. Install with the optional extra: `pip install "pyeast[agent]"`.
- **Index or name selection in `gg` and `batch`.** Parts and instruction files can be
  selected by list index as well as by name.

### Changed

- Plasmid-map colour palette updated.
- Dependencies modernised to current versions.

### Removed

- Bundled `data/` directory (now the separate PYEAST_data repository - see Breaking
  changes).
- Unused runtime dependencies `loguru` and `pyfastcopy`.
- `python-dotenv` moved from the core dependencies into the optional `agent` extra.

[2.0.0]: https://github.com/TomLoan/PYEAST/releases/tag/v2.0.0
