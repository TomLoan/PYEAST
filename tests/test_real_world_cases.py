# tests/test_real_world_cases.py

import pytest
import os 
import pandas as pd
import pytest
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO

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



def test_tar_designer_real_world_case(): 
    '''Test TARDesigner against a known real world case.'''
    
    test_case_dir = 'known_good_tar_case'

    #skip if test data doesn't exist yet 
    input_dir = Path('test/data')/test_case_dir/'input'
    if not input_dir.exisits() or not any(input_dir.glob('*.fasta')): 
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

def test_tar_designer_full_workflow(): 
    '''test the complete TAR workflow including instruction generation'''
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

    #Mock tmeplate and primer locations for testing 
    designer.primers_found = {name: [{'Location': 'TestPlate', 'Position': 'A1', 'sequence': 'seq'}] for name, seq in designer.primers.items()}    
    designer.missing_primers = {}
    designer.template_dict = {seq.name: ["TestTemplate"] for seq in sequences.values()}

    #Run rationalisation and instruction generation 
    designer.rationalize_selections()
    instructions = designer.write_instructions()

    #Basic validation of instructions format
    assert len(instructions) == len(sequences, "Should have instructions for each sequence")

    #Check that instructions have the correct structure 
    for instruction in instructions:
        assert len(instruction) == 9, 'Each instruction should have 9 elements'
        assert instruction[0] in [seq.name for seq in sequences.values()], 'Part name should match sequnece name'

def test_assembly_creation(): 
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
