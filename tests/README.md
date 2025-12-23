# PYEAST Test Suite

## Golden Master Tests

These tests preserve PYEAST behavior during the PyDNA refactoring. They verify that:
- Final assembled DNA sequences remain identical
- Correct number of primers are designed
- Output file formats remain stable
- Core functionality doesn't break

**Location:** `tests/fixtures/golden_*/`

### Test Cases

1. **golden_tar**: 6-part circular plasmid assembly using TAR
   - Assembly: pTEF1 → YeRFP → tDIT1 → Ura3 → AmpR_ColE1 → 2Micron
   - Expected output: 5722 bp circular plasmid, 12 primers
   - Validates: Primer design, circular assembly logic

2. **golden_integrate**: Chromosomal integration cassette
   - Components: pTEF1 → YeGFP → tDIT1
   - Integration site: Ura3MX (upstream/downstream arms)
   - Expected output: 3630 bp linear construct, 8 primers
   - Validates: Linear assembly, integration site handling

## Running Tests

```bash
# All tests
uv run pytest tests/ -v

# Just golden master tests
uv run pytest tests/ -k golden -v

# Specific test
uv run pytest tests/test_real_world_cases.py::test_tar_designer_real_world_case -v

# With coverage
uv run pytest tests/ --cov=src/pyeast --cov-report=html
```

## Test Structure

```
tests/
├── README.md                    (this file)
├── test_real_world_cases.py     (golden master tests)
├── fixtures/
│   ├── golden_tar/             (TAR test data)
│   │   ├── input/              (6 FASTA files)
│   │   ├── expected/           (primers.tsv, instructions.tsv, etc.)
│   │   └── README.md
│   └── golden_integrate/       (Integration test data)
│       ├── input/              (4 FASTA files)
│       ├── expected/           (primers.tsv, instructions.tsv, etc.)
│       └── README.md
└── unit/                       (placeholder for future unit tests)
```

## During PyDNA Refactoring

**Run golden master tests frequently:**
```bash
uv run pytest tests/ -k golden -v
```

**Interpreting Results:**

- ✅ **All pass:** Refactor preserved core behavior - safe to continue
- ⚠️ **Format tests fail:** Output structure changed (update tests if intentional)
- ⚠️ **Primer count changed:** Number of primers differs (investigate - might be valid optimization)
- ❌ **Assembly sequence differs:** DNA sequence changed - STOP and investigate! This means the refactor broke assembly logic.

## What the Tests Check

### Core Behavior (MUST NOT change):
- **Assembled DNA sequence**: Exact nucleotide sequence of final assembly
- **Assembly topology**: Circular for TAR, linear for Integration
- **Assembly length**: 5722 bp for TAR, 3630 bp for Integration
- **Primer count**: 12 for TAR (6 parts × 2), 8 for Integration (3 components × 2 + 2 sites)

### Implementation Details (CAN change):
- Individual primer sequences (PyDNA may use different algorithms)
- Primer naming conventions
- Exact primer binding positions
- Melting temperatures or other primer properties

The key principle: **If the assembled DNA sequence is correct, the primers work correctly.**

## After Refactoring

- Verify all golden master tests still pass
- Update fixtures if intentional changes were made
- Keep golden masters as regression tests
