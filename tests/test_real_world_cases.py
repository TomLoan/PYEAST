# tests/test_real_world_cases.py

import pytest
import os 
import pandas as pd
import pytest
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from unittest.mock import Mock, patch

from pyeast.core.tar import TARDesigner
from pyeast.utils.sequence_utils import load_sequences

#Helper functions to load test data 
def load_expected_primers(test_case_dir): 
    '''Load expected primers from a tsv file'''
    primers_file = Path('test/data')/test_case_dir/'expected_outputs'/'primers.tsv'
    if not primers_file.exists(): 
        return {}

    df = pd.read_csv(primers_file, sep = '\t')
    return dict(zip(df['Name'], df['Sequence']))    

def load_expected_instructions(test_case_dir):
    '''Load expected instructions from TSV file'''
    instructions_file = Path('tests/data')/test_case_dir/'expected_output'/'instructions.tsv'
    if not instructions_file.exists(): 
        return []
    df = pd.read_csv(instructions_file, sep = '\t')
    return df.values.tolist()

@pytest.fixture
def mock_console_interaction():
    """Mock prompt_toolkit components to avoid console interaction during tests"""
    with patch('pyeast.core.tar.PromptSession') as mock_session, \
         patch('pyeast.core.tar.Console') as mock_console, \
         patch('pyeast.core.tar.confirm') as mock_confirm:
        
        # Mock the session to avoid prompt_toolkit calls
        mock_session_instance = Mock()
        mock_session.return_value = mock_session_instance
        
        # Mock console to avoid rich output issues
        mock_console_instance = Mock()
        mock_console.return_value = mock_console_instance
        
        # Mock confirm to always return True
        mock_confirm.return_value = True
        
        yield {
            'session': mock_session_instance,
            'console': mock_console_instance,
            'confirm': mock_confirm
        }


def test_tar_designer_real_world_case(mock_console_interaction): 
    '''Test TARDesigner against a known real world case.'''
    
    test_case_dir = 'known_good_tar_case'

    #skip if test data doesn't exist yet 
    input_dir = Path('tests/data')/test_case_dir/'input'
    if not input_dir.exists() or not any(input_dir.glob('*.fasta')): 
        pytest.skip('Test data not yet set up')

    #Load test sequences 

    sequences = load_sequences(input_dir)
    expected_primers = load_expected_primers(test_case_dir)

    designer = TARDesigner(homology_length=25)
    designer.available_sequences = sequences

    #defined assembly order matches the order in the example case 
    assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']

    designer.set_assembly_order(assembly_order)
    designer.design_tar_primers()

    #test primer seqeunces 
    if expected_primers: 
        for name, expected_sequence in expected_primers.items(): 
            assert name in designer.primers, f'Expected primer {name} not found'
            actaul_sequence = str(designer.primers[name])
            assert actaul_sequence == expected_sequence, f'Primer {name}: {expected_sequence}. got {actaul_sequence}'

def test_tar_designer_full_workflow(mock_console_interaction): 
    '''Test the complete TAR workflow including instruction generation'''
    test_case_dir = 'known_good_tar_case'
    input_dir = Path('tests/data')/test_case_dir/'input'
    
    # Skip if test data doesn't exist
    if not input_dir.exists() or not any(input_dir.glob('*.fasta')): 
        pytest.skip('Test data not yet set up')
        
    assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']

    # Load sequences and run complete workflow
    sequences = load_sequences(input_dir)
    designer = TARDesigner(homology_length=25)
    designer.available_sequences = sequences
    designer.set_assembly_order(assembly_order)
    designer.design_tar_primers()


    # Mock template and primer locations for testing 
    designer.check_primer_locations(Path('tests/data')/test_case_dir/'input')
    designer.template_dict = {seq.name: ["TestTemplate"] for seq in sequences.values()}


    # Run rationalization and instruction generation 
    designer.rationalize_selections()
    instructions = designer.write_instructions()

    # Basic validation of instructions format
    assert len(instructions) == len(sequences), "Should have instructions for each sequence"

    # Check that instructions have the correct structure 
    for instruction in instructions:
        assert len(instruction) == 9, 'Each instruction should have 9 elements'
        assert instruction[0] in [seq.name for seq in sequences.values()], 'Part name should match sequence name'

def test_assembly_creation(mock_console_interaction): 
    '''Test that the final assembly is correct'''  
    test_case_dir = 'known_good_tar_case'
    input_dir = Path('tests/data')/test_case_dir/'input'
    assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']

    sequences = load_sequences(input_dir)

    #load_seqeunces and run complete work flow

    sequences = load_sequences(input_dir)
    designer = TARDesigner(homology_length=25)
    designer.available_sequences = sequences
    designer.set_assembly_order(assembly_order)
    designer.design_tar_primers()

    #create the final assembly 
    assembly = designer.create_assembly()

    #basic validation 
    assert assembly is not None, 'Assembly aught to be created'
    assert len(assembly.seq)>0, 'Assembly should have a sequence'
    assert assembly.annotations.get('topology') =='circular', ' TAR assembly aught to be circular'  

    #chack that all input sequences are represeented in the assmebly 
    assembly_str = str(assembly.seq)
    for seq in sequences.values(): 
        assert str(seq.seq) in assembly_str, f'Sequence {seq.name} should be in the assembly'

def test_tar_designer_without_interactive_methods(mock_console_interaction):
    '''Test TARDesigner using only non-interactive methods'''
    # This test doesn't need mocking since it doesn't call interactive methods
    
    test_case_dir = 'known_good_tar_case'
    input_dir = Path('tests/data')/test_case_dir/'input'
    
    # Skip if test data doesn't exist
    if not input_dir.exists() or not any(input_dir.glob('*.fasta')): 
        pytest.skip('Test data not yet set up')

    # Load test sequences 
    sequences = load_sequences(input_dir)
    
    # Create designer without any interactive calls
    designer = TARDesigner(homology_length=25)
    
    # Test basic functionality
    assert designer.homology_length == 25
    assert designer.annealing_temp == 50  # default value
    
    # Test sequence loading
    designer.available_sequences = sequences
    assert len(designer.available_sequences) > 0
    
    # Test assembly order setting
    assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']
    designer.set_assembly_order(assembly_order)
    assert len(designer.assembly_sequences) == len(assembly_order)
    
    # Test primer design
    designer.design_tar_primers()
    assert len(designer.primers) > 0
    
    # Verify primer naming convention
    for primer_name in designer.primers.keys():
        assert '_' in primer_name, f"Primer {primer_name} should contain underscore separator"

def test_primer_validation(mock_console_interaction):
    '''Test primer sequence validation and properties'''
    test_case_dir = 'known_good_tar_case'
    input_dir = Path('tests/data')/test_case_dir/'input'
    
    # Skip if test data doesn't exist
    if not input_dir.exists() or not any(input_dir.glob('*.fasta')): 
        pytest.skip('Test data not yet set up')

    sequences = load_sequences(input_dir)
    designer = TARDesigner(homology_length=25)
    designer.available_sequences = sequences
    
    assembly_order = ['pTEF1', 'YeRFP', 'tDIT1', 'Ura3', 'AmpR_ColE1', '2Micron']
    designer.set_assembly_order(assembly_order)
    designer.design_tar_primers()
    
    # Test primer properties
    for primer_name, primer_seq in designer.primers.items():
        # Check primer length is reasonable (typically 40-80 bp for TAR)
        assert 30 <= len(primer_seq) <= 100, f"Primer {primer_name} length {len(primer_seq)} is outside expected range"
        
        # Check primer contains only valid DNA bases
        valid_bases = set('ATCG')
        primer_bases = set(str(primer_seq).upper())
        assert primer_bases.issubset(valid_bases), f"Primer {primer_name} contains invalid bases: {primer_bases - valid_bases}"
        
        # Check primer name format (should contain sequence names and direction)
        assert any(seq_name in primer_name for seq_name in assembly_order), f"Primer {primer_name} should reference assembly sequences"