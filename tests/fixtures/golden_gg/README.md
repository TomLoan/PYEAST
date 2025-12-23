# Golden GG (Golden Gate) Test Case

**Assembly Type:** Golden Gate (Type IIS restriction-based cloning)
**Parts:** 4 components (using multi-part fragments)
**Assembly Order:** 1_ConL1, 234_Spacer, 5_ConR2, 678_AmpR-ColE1
**Expected Output:** 2278 bp circular plasmid
**Description:** Minimal MoClo level 1 assembly using multi-part fragments

## Assembly Details

This test uses a simplified 4-part Golden Gate assembly to minimize complexity while still validating core GG functionality.

**Component List:**
- `1_ConL1` - Left connector (MoClo part type 1)
- `234_Spacer` - Multi-part fragment combining promoter/CDS/terminator (types 2-3-4)
- `5_ConR2` - Right connector (MoClo part type 5)
- `678_AmpR-ColE1` - Multi-part fragment combining marker/origin/resistance (types 6-7-8)

**Plasmid Sources (from Yeast MoClo lvl 0 library):**
- pYTK003 (1_ConL1) - Well H11
- pYTK048 (234_Spacer) - Well A3
- pYTK068 (5_ConR2) - Well D12
- pYTK095 (678_AmpR-ColE1) - Well F8

## Expected Outputs

- **golden_gg.gb** - 2278 bp circular assembled construct
- **golden_gg_worklist.csv** - Janus liquid handler instructions (4 plasmid transfers)
- **golden_gg_summary.csv** - Assembly summary information

**Note:** GG assembly does not generate screening primers (unlike TAR/Integration).

## Test Validation

This test verifies that Golden Gate assembly remains stable during PyDNA refactoring by checking:
1. Final assembled DNA sequence is correct (2278 bp)
2. Assembly topology is circular
3. Part-to-plasmid mapping succeeds for all 4 components
4. Assembly produces exactly 1 construct

## Dependencies

This test uses existing component libraries:
- Component sequences: `data/component libraries/Yeast Moclo lvl 0/`
- Plasmid files: `data/component libraries/Yeast Moclo lvl 0/plasmids/`

The test does NOT replicate these files in the fixtures directory to avoid duplication.
