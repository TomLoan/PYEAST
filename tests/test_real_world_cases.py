# tests/test_real_world_cases.py

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from Bio import SeqIO

from pyeast.core.deletion import DeletionDesigner
from pyeast.core.gg import ggDesigner
from pyeast.core.integration import IntegrationDesigner
from pyeast.core.replace import ReplaceDesigner
from pyeast.core.tar import TARDesigner
from pyeast.utils.path_utils import get_component_libraries_path, get_templates_path
from pyeast.utils.sequence_utils import load_sequences


#Helper functions to load test data
def load_expected_primers(test_case_dir):
    '''Load expected primers from a tsv file'''
    primers_file = Path('tests/fixtures')/test_case_dir/'expected'/'primers.tsv'
    if not primers_file.exists():
        print("didn't find some primers")
        return {}

    df = pd.read_csv(primers_file, sep = '\t')
    return dict(zip(df['Name'], df['Sequence']))

def load_expected_instructions(test_case_dir):
    '''Load expected instructions from TSV file'''
    instructions_file = Path('tests/fixtures')/test_case_dir/'expected'/'instructions.tsv'
    if not instructions_file.exists():
        return []
    df = pd.read_csv(instructions_file, sep = '\t')
    return df.to_numpy()

def load_expected_assembly(test_case_dir, assembly_type='tar'):
    """Load expected assembled sequence from GenBank file."""
    if assembly_type == 'tar':
        gb_file = Path('tests/fixtures')/test_case_dir/'expected'/'tar_test.gb'
    else:  # integration
        gb_file = Path('tests/fixtures')/test_case_dir/'expected'/'Integrate_test.gb'

    if not gb_file.exists():
        return None

    record = SeqIO.read(gb_file, 'genbank')
    return record

@pytest.fixture
def mock_console_interaction():
    """Mock click.confirm to avoid interactive prompts during tests.

    TAR, Integration, Deletion, Replace, Batch, and GG designers no longer use
    interactive components directly - those live in cli/main.py.
    """
    with patch('click.confirm') as mock_click_confirm:
        mock_click_confirm.return_value = True

        yield {
            'click_confirm': mock_click_confirm
        }

def test_tar_designer_real_world_case(mock_console_interaction):
    '''Test TARDesigner against a known real world case.'''

    test_case_dir = 'golden_tar'

    #skip if test data doesn't exist yet
    input_dir = Path('tests/fixtures')/test_case_dir/'input'
    if not input_dir.exists() or not any(input_dir.glob('*.fasta')):
        pytest.skip('Test data not yet set up')

    #Load test sequences
    sequences = load_sequences(input_dir)

    designer = TARDesigner(homology_length=25)
    designer.available_sequences = sequences

    #defined assembly order matches the order in the example case
    assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']

    designer.set_assembly_order(assembly_order)
    designer.design_tar_primers()

    # Test 1: Check primer count (should have 2 primers per part for circular assembly)
    expected_primer_count = len(assembly_order) * 2  # 6 parts * 2 = 12 primers
    assert len(designer.primers) == expected_primer_count, \
        f"Expected {expected_primer_count} primers, got {len(designer.primers)}"

    # Test 2: Validate assembled sequence matches expected
    actual_assembly = designer.create_assembly()
    expected_assembly = load_expected_assembly('golden_tar', 'tar')

    assert actual_assembly is not None, "Failed to create assembly"

    if expected_assembly:
        # Compare sequences case-insensitively (case may vary in FASTA files)
        assert str(actual_assembly.seq).upper() == str(expected_assembly.seq).upper(), \
            "Assembled DNA sequence doesn't match expected golden master"
        assert actual_assembly.annotations.get('topology') == 'circular', \
            "TAR assembly should be circular"
        assert len(actual_assembly.seq) == 5722, \
            f"Expected 5722 bp assembly, got {len(actual_assembly.seq)}"

# def test_tar_designer_full_workflow(mock_console_interaction):
#     '''Test the complete TAR workflow including instruction generation'''
#     test_case_dir = 'known_good_tar_case'
#     input_dir = Path('tests/data')/test_case_dir/'input'

#     # Skip if test data doesn't exist
#     if not input_dir.exists() or not any(input_dir.glob('*.fasta')):
#         pytest.skip('Test data not yet set up')

#     assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']

#     # Load sequences and run complete workflow
#     sequences = load_sequences(input_dir)
#     designer = TARDesigner(homology_length=25)
#     designer.available_sequences = sequences
#     designer.set_assembly_order(assembly_order)
#     designer.design_tar_primers()


#     # Mock template and primer locations for testing
#     designer.check_primer_locations(Path('tests/data')/test_case_dir/'input')
#     designer.template_dict = {seq.name: ["TestTemplate"] for seq in sequences.values()}


#     # Run rationalization and instruction generation
#     designer.rationalize_selections()
#     instructions = designer.write_instructions()

#     # Basic validation of instructions format
#     assert len(instructions) == len(sequences), "Should have instructions for each sequence"

#     # Check that instructions have the correct structure
#     for instruction in instructions:
#         assert len(instruction) == 9, 'Each instruction should have 9 elements'
#         assert instruction[0] in [seq.name for seq in sequences.values()], 'Part name should match sequence name'

# def test_assembly_creation(mock_console_interaction):
#     '''Test that the final assembly is correct'''
#     test_case_dir = 'known_good_tar_case'
#     input_dir = Path('tests/data')/test_case_dir/'input'
#     assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']

#     sequences = load_sequences(input_dir)

#     #load_seqeunces and run complete work flow

#     sequences = load_sequences(input_dir)
#     designer = TARDesigner(homology_length=25)
#     designer.available_sequences = sequences
#     designer.set_assembly_order(assembly_order)
#     designer.design_tar_primers()

#     #create the final assembly
#     assembly = designer.create_assembly()

#     #basic validation
#     assert assembly is not None, 'Assembly aught to be created'
#     assert len(assembly.seq)>0, 'Assembly should have a sequence'
#     assert assembly.annotations.get('topology') =='circular', ' TAR assembly aught to be circular'

#     #chack that all input sequences are represeented in the assmebly
#     assembly_str = str(assembly.seq)
#     for seq in sequences.values():
#         assert str(seq.seq) in assembly_str, f'Sequence {seq.name} should be in the assembly'

# def test_tar_designer_without_interactive_methods(mock_console_interaction):
#     '''Test TARDesigner using only non-interactive methods'''
#     # This test doesn't need mocking since it doesn't call interactive methods

#     test_case_dir = 'known_good_tar_case'
#     input_dir = Path('tests/data')/test_case_dir/'input'

#     # Skip if test data doesn't exist
#     if not input_dir.exists() or not any(input_dir.glob('*.fasta')):
#         pytest.skip('Test data not yet set up')

#     # Load test sequences
#     sequences = load_sequences(input_dir)

#     # Create designer without any interactive calls
#     designer = TARDesigner(homology_length=25)

#     # Test basic functionality
#     assert designer.homology_length == 25
#     assert designer.annealing_temp == 50  # default value

#     # Test sequence loading
#     designer.available_sequences = sequences
#     assert len(designer.available_sequences) > 0

#     # Test assembly order setting
#     assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']
#     designer.set_assembly_order(assembly_order)
#     assert len(designer.assembly_sequences) == len(assembly_order)

#     # Test primer design
#     designer.design_tar_primers()
#     assert len(designer.primers) > 0

#     # Verify primer naming convention
#     for primer_name in designer.primers.keys():
#         assert '_' in primer_name, f"Primer {primer_name} should contain underscore separator"

# def test_primer_validation(mock_console_interaction):
#     '''Test primer sequence validation and properties'''
#     test_case_dir = 'known_good_tar_case'
#     input_dir = Path('tests/data')/test_case_dir/'input'

#     # Skip if test data doesn't exist
#     if not input_dir.exists() or not any(input_dir.glob('*.fasta')):
#         pytest.skip('Test data not yet set up')

#     sequences = load_sequences(input_dir)
#     designer = TARDesigner(homology_length=25)
#     designer.available_sequences = sequences

#     assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']
#     designer.set_assembly_order(assembly_order)
#     designer.design_tar_primers()

#     # Test primer properties
#     for primer_name, primer_seq in designer.primers.items():
#         # Check primer length is reasonable (typically 40-80 bp for TAR)
#         assert 30 <= len(primer_seq) <= 100, f"Primer {primer_name} length {len(primer_seq)} is outside expected range"

#         # Check primer contains only valid DNA bases
#         valid_bases = set('ATCG')
#         primer_bases = set(str(primer_seq).upper())
#         assert primer_bases.issubset(valid_bases), f"Primer {primer_name} contains invalid bases: {primer_bases - valid_bases}"

#         # Check primer name format (should contain sequence names and direction)
#         assert any(seq_name in primer_name for seq_name in assembly_order), f"Primer {primer_name} should reference assembly sequences"

### Claude tests for integrate follow - re-written by hand.
def test_integration_designer_real_world_case(mock_console_interaction):
    """Test IntegrationDesigner against a known real world case."""
    test_case_dir = 'golden_integrate'


    # Skip if test data doesn't exist yet
    input_dir = Path('tests/fixtures')/test_case_dir/'input'
    if not input_dir.exists() or not any(input_dir.glob('*.fasta')):
        pytest.skip('Test data not yet set up')

    # Load test sequences
    sequences = load_sequences(input_dir)

    # Create designer
    designer = IntegrationDesigner(homology_length=25)

    # Set up sequences directly (bypass interactive methods)
    # Components: pTEF1, YeGFP, tDIT1 (biologically sensible integration)
    component_names = ['pTEF1', 'YeGFP', 'tDIT1']

    # Set components
    designer.components = {}
    for name in component_names:
        if name in sequences:
            designer.components[name] = sequences[name]

    # Set up Ura3MX integration site (assuming it exists in sequences)
    # Look for Ura3MX upstream and downstream sequences
    ura3mx_upstream = None
    ura3mx_downstream = None

    for seq_name, seq_record in sequences.items():
        if 'ura3mx' in seq_name.lower():
            if 'upstream' in seq_name.lower() or 'up' in seq_name.lower():
                ura3mx_upstream = seq_record
            elif 'downstream' in seq_name.lower() or 'down' in seq_name.lower():
                ura3mx_downstream = seq_record

    # Set up integration site
    if ura3mx_upstream and ura3mx_downstream:
        designer.int_sites = {'Ura3MX': (ura3mx_upstream, ura3mx_downstream)}
        designer.int_site = 'Ura3MX'


        # Set assembly sequences: upstream + components + downstream
        component_seqs = [designer.components[name] for name in component_names if name in designer.components]
        designer.assembly_sequences = [ura3mx_upstream] + component_seqs + [ura3mx_downstream]

        # Design integration primers (this should be non-interactive)
        designer.design_integration_primers()

        # Test 1: Check primer count
        # Integration creates: 2 primers for each component + 2 primers for each integration site arm
        # = (3 components × 2) + (2 integration arms × 2) = 6 + 4 = 10 primers
        expected_primer_count = 10
        assert len(designer.primers) == expected_primer_count, \
            f"Expected {expected_primer_count} primers, got {len(designer.primers)}"

        # Test 2: Validate assembled sequence matches expected
        actual_assembly = designer.create_linear_assembly()
        expected_assembly = load_expected_assembly('golden_integrate', 'integration')

        if expected_assembly:
            # Compare sequences case-insensitively (case may vary in FASTA files)
            assert str(actual_assembly.seq).upper() == str(expected_assembly.seq).upper(), \
                "Assembled DNA sequence doesn't match expected golden master"
            assert actual_assembly.annotations.get('topology') == 'linear', \
                "Integration assembly should be linear"
            assert len(actual_assembly.seq) == 3630, \
                f"Expected 3630 bp assembly, got {len(actual_assembly.seq)}"
    else:
        pytest.skip('Ura3MX integration site sequences not found in test data')


def test_primers_tsv_format():
    """Validate that primer TSV files have correct structure."""
    # Test both TAR and Integration primer outputs
    test_cases = ['golden_tar', 'golden_integrate']

    for test_case in test_cases:
        primers_file = Path('tests/fixtures')/test_case/'expected'/'primers.tsv'

        # Skip if file doesn't exist
        if not primers_file.exists():
            pytest.skip(f'Primers file not found for {test_case}')

        # Load TSV
        df = pd.read_csv(primers_file, sep='\t')

        # Required columns check
        assert 'Name' in df.columns, f"{test_case}: Missing 'Name' column"
        assert 'Sequence' in df.columns, f"{test_case}: Missing 'Sequence' column"

        # No duplicate primer names
        assert len(df['Name']) == len(df['Name'].unique()), \
            f"{test_case}: Found duplicate primer names"

        # Valid DNA sequences (ATGC only)
        for idx, seq in enumerate(df['Sequence']):
            assert len(seq) > 0, f"{test_case}: Empty sequence at row {idx}"
            assert all(base in 'ATGC' for base in str(seq).upper()), \
                f"{test_case}: Invalid DNA sequence at row {idx}: {seq}"


def test_instructions_tsv_format():
    """Validate that instructions TSV files have correct structure."""
    test_cases = ['golden_tar', 'golden_integrate']

    for test_case in test_cases:
        instructions_file = Path('tests/fixtures')/test_case/'expected'/'instructions.tsv'

        # Skip if file doesn't exist
        if not instructions_file.exists():
            pytest.skip(f'Instructions file not found for {test_case}')

        # Load TSV
        df = pd.read_csv(instructions_file, sep='\t')

        # Required columns (9 total from TARDesigner/IntegrationDesigner)
        expected_columns = ['Part', 'F_Plate', 'F_Name', 'F_Well',
                           'R_Plate', 'R_Name', 'R_Well', 'Template', 'Size']
        for col in expected_columns:
            assert col in df.columns, f"{test_case}: Missing column: {col}"

        # Validate data types
        assert df['Size'].dtype in [int, 'int64', 'float64'], \
            f"{test_case}: Size column should be numeric"
        assert all(df['Part'].notna()), \
            f"{test_case}: Part names should not be empty"


@pytest.mark.realdata
def test_deletion_designer_golden_case(mock_console_interaction):
    """Test DeletionDesigner against a known deletion case."""
    test_case_dir = 'golden_delete'

    # Skip if test data doesn't exist yet
    input_file = Path('tests/fixtures')/test_case_dir/'input'/'input.fasta'
    if not input_file.exists():
        pytest.skip('Test data not yet set up')

    # Load input sequence (the target gene to delete)
    input_seq = SeqIO.read(input_file, 'fasta')

    # Create designer with parameters matching golden case
    designer = DeletionDesigner(
        upstream_homology_len=300,
        downstream_homology_len=100,
        repeat_length=80
    )

    # Set the target sequence
    designer.target_sequence = str(input_seq.seq)

    # Find target in genome (using default genome file in templates/)
    genome_file = get_templates_path() / 'BY4741_Toronto_2012.fsa'
    if not genome_file.exists():
        pytest.skip('Genome file not found - required for deletion tests')

    location = designer.find_target_sequence(genome_file, designer.target_sequence)
    assert location is not None, "Target sequence should be found in genome"
    designer.genome_location = location

    # Create deletion cassette
    ura3_file = get_component_libraries_path() / 'Saccharomyces_cerevisiae' / 'URA3.fasta'
    if not ura3_file.exists():
        pytest.skip('URA3 file not found - required for deletion tests')

    cassette = designer.make_deletion_cassette(genome_file, ura3_file)
    assert cassette is not None, "Failed to create deletion cassette"

    # Test 1: Validate cassette sequence matches expected
    expected_cassette = SeqIO.read(
        Path('tests/fixtures')/test_case_dir/'expected'/'golden_delete.gb',
        'genbank'
    )

    # Compare sequences case-insensitively (case may vary in FASTA files)
    assert str(cassette.seq).upper() == str(expected_cassette.seq).upper(), \
        "Deletion cassette sequence doesn't match golden master"
    assert len(cassette.seq) == 1533, \
        f"Expected 1533 bp cassette, got {len(cassette.seq)}"

    # Test 2: Validate cassette has correct feature structure
    feature_labels = [f.qualifiers.get('label', [''])[0] if isinstance(f.qualifiers.get('label'), list)
                     else f.qualifiers.get('label', '') for f in cassette.features]
    assert any('upstream' in str(label).lower() for label in feature_labels), \
        f"Cassette should have upstream homology feature. Found labels: {feature_labels}"
    assert any('repeat' in str(label).lower() for label in feature_labels), \
        f"Cassette should have repeat feature. Found labels: {feature_labels}"
    assert any('ura3' in str(label).lower() for label in feature_labels), \
        f"Cassette should have URA3 marker feature. Found labels: {feature_labels}"
    assert any('downstream' in str(label).lower() for label in feature_labels), \
        f"Cassette should have downstream homology feature. Found labels: {feature_labels}"

    # Test 3: Design screening primers
    forward_primer, reverse_primer, product_sizes = designer.design_screening_strategy(genome_file)

    # Check that 2 primers were designed (not checking specific sequences)
    assert forward_primer is not None, "Forward primer should be designed"
    assert reverse_primer is not None, "Reverse primer should be designed"
    assert len(forward_primer) > 0, "Forward primer should not be empty"
    assert len(reverse_primer) > 0, "Reverse primer should not be empty"


@pytest.mark.realdata
def test_replace_designer_golden_case(mock_console_interaction):
    """Test ReplaceDesigner against a known replacement case."""
    test_case_dir = 'golden_replace'

    # Skip if test data doesn't exist yet
    input_dir = Path('tests/fixtures')/test_case_dir/'input'
    target_file = input_dir/'input.fasta'
    replacement_file = input_dir/'pTEF1.fasta'

    if not target_file.exists() or not replacement_file.exists():
        pytest.skip('Test data not yet set up')

    # Load input sequences
    target_seq = SeqIO.read(target_file, 'fasta')
    replacement_seq = SeqIO.read(replacement_file, 'fasta')

    # Create designer with parameters matching golden case
    # From README: marker position is "up" (upstream)
    designer = ReplaceDesigner(
        upstream_homology_len=200,
        downstream_homology_len=200,
        repeat_length=80
    )

    # Set the target sequence and replacement sequence
    designer.target_sequence = str(target_seq.seq)
    designer.replacement_sequence = replacement_seq
    designer.marker_position = "upstream"  # From README

    # Find target in genome (using default genome file in templates/)
    genome_file = get_templates_path() / 'BY4741_Toronto_2012.fsa'
    if not genome_file.exists():
        pytest.skip('Genome file not found - required for replace tests')

    location = designer.find_target_sequence(genome_file, designer.target_sequence)
    assert location is not None, "Target sequence should be found in genome"
    designer.genome_location = location

    # Create replacement cassette
    ura3_file = get_component_libraries_path() / 'Saccharomyces_cerevisiae' / 'URA3.fasta'
    if not ura3_file.exists():
        pytest.skip('URA3 file not found - required for replace tests')

    cassette = designer.make_replacement_cassette(genome_file, ura3_file)
    assert cassette is not None, "Failed to create replacement cassette"

    # Test 1: Validate cassette sequence matches expected
    expected_cassette = SeqIO.read(
        Path('tests/fixtures')/test_case_dir/'expected'/'golden_replace.gb',
        'genbank'
    )

    # Compare sequences case-insensitively (case may vary in FASTA files)
    assert str(cassette.seq).upper() == str(expected_cassette.seq).upper(), \
        "Replacement cassette sequence doesn't match golden master"
    assert len(cassette.seq) == 1945, \
        f"Expected 1945 bp cassette, got {len(cassette.seq)}"

    # Test 2: Validate cassette has correct feature structure
    # Expected features: upstream homology, URA3, repeat, replacement, downstream homology
    feature_labels = [f.qualifiers.get('label', [''])[0] if isinstance(f.qualifiers.get('label'), list)
                     else f.qualifiers.get('label', '') for f in cassette.features]

    assert any('upstream' in str(label).lower() for label in feature_labels), \
        f"Cassette should have upstream homology feature. Found labels: {feature_labels}"
    assert any('ura3' in str(label).lower() for label in feature_labels), \
        f"Cassette should have URA3 marker feature. Found labels: {feature_labels}"
    assert any('repeat' in str(label).lower() for label in feature_labels), \
        f"Cassette should have repeat feature. Found labels: {feature_labels}"
    assert any('replacement' in str(label).lower() for label in feature_labels), \
        f"Cassette should have replacement sequence feature. Found labels: {feature_labels}"
    assert any('downstream' in str(label).lower() for label in feature_labels), \
        f"Cassette should have downstream homology feature. Found labels: {feature_labels}"

    # Test 3: Design screening primers
    forward_primer, reverse_primer, product_sizes = designer.design_screening_strategy(genome_file)

    # Check that 2 primers were designed (not checking specific sequences)
    assert forward_primer is not None, "Forward primer should be designed"
    assert reverse_primer is not None, "Reverse primer should be designed"
    assert len(forward_primer) > 0, "Forward primer should not be empty"
    assert len(reverse_primer) > 0, "Reverse primer should not be empty"


@pytest.mark.realdata
def test_gg_designer_golden_case(mock_console_interaction):
    """Test ggDesigner (Golden Gate) against a known assembly case."""
    test_case_dir = 'golden_gg'

    # Skip if test data doesn't exist yet
    input_file = Path('tests/fixtures')/test_case_dir/'input'/'input.txt'
    if not input_file.exists():
        pytest.skip('Test data not yet set up')

    # Load component sequences from existing MoClo library
    component_lib = get_component_libraries_path() / 'Yeast_Moclo_lvl_0'
    if not component_lib.exists():
        pytest.skip('MoClo component library not found - required for GG tests')

    sequences = load_sequences(component_lib)

    # Create designer
    designer = ggDesigner(is_library=False)
    designer.available_sequences = sequences

    # Set plasmid directory (critical for GG!)
    designer.gg_plasmids = component_lib / 'plasmids'
    if not designer.gg_plasmids.exists():
        pytest.skip('MoClo plasmid directory not found - required for GG tests')

    # Set assembly order (from input.txt)
    # Assembly: 1_ConL1, 234_Spacer, 5_ConR2, 678_AmpR-ColE1
    assembly_order = ['1_ConL1', '234_Spacer', '5_ConR2', '678_AmpR-ColE1']
    designer.assemblies_names = [assembly_order]  # List of assemblies (single assembly in this case)

    # Get plasmid names (maps parts to plasmids)
    try:
        plasmid_names = designer.get_plasmid_names()
        assert plasmid_names is not None, "Failed to get plasmid names"
        assert len(plasmid_names) == 4, f"Expected 4 plasmids, got {len(plasmid_names)}"
    except Exception as e:
        pytest.fail(f"Part-to-plasmid mapping failed: {e}")

    # Run GG assembly simulation
    try:
        assembly = designer.gg_assembly('golden_gg')
        assert assembly is not None, "Failed to create assembly"
    except Exception as e:
        pytest.fail(f"GG assembly simulation failed: {e}")

    # Validate assembly simulation
    assert hasattr(designer, 'assembly_sim'), "Assembly simulation not created"
    assert designer.assembly_sim is not None, "Assembly simulation is None"

    # Test 1: Check construct count (should be exactly 1)
    construct_records = designer.assembly_sim.construct_records
    assert len(construct_records) == 1, \
        f"Expected 1 construct, got {len(construct_records)}"

    # Test 2: Validate assembled sequence matches expected
    actual_assembly = construct_records[0]
    expected_assembly = SeqIO.read(
        Path('tests/fixtures')/test_case_dir/'expected'/'golden_gg.gb',
        'genbank'
    )

    # Compare sequences case-insensitively (case may vary in FASTA files)
    assert str(actual_assembly.seq).upper() == str(expected_assembly.seq).upper(), \
        "Assembled DNA sequence doesn't match golden master"

    # Test 3: Validate assembly properties
    assert len(actual_assembly.seq) == 2278, \
        f"Expected 2278 bp assembly, got {len(actual_assembly.seq)}"

    # Check topology (circular for MoClo level 1)
    topology = actual_assembly.annotations.get('topology', 'linear')
    assert topology == 'circular', \
        f"Expected circular topology, got {topology}"
