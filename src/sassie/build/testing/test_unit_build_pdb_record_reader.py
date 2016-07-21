# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 11:12:58 2015

@author: dave
"""

from unittest import main 

import warnings; warnings.filterwarnings('ignore')

import pdbscan.pdb_record_reader as record_reader

def test_runon_line():
    """
    Test the combination of multiple lines from a PDB as parsed by 
    record_reader
    """
    
    test_field = 'text'
    test_text1 = 'So long'
    test_text2 = '& thanks for all the fish'
    
    # Note: Due to PDB format not guaranteeing correct spacing one is inserted
    expected_result = test_text1 + ' ' + '& thanks for all the fish'
    
    
    test_lines = [{'part':'1', test_field: test_text1}, 
                  {'part':'2', test_field: test_text2}]
    
    combined = record_reader.process_runon_line(test_lines, test_field)
        
    assert combined == expected_result
    
def test_parse_seqres_line():

    schema = record_reader.rec_schemas['SEQRES']

    line = 'SEQRES   1 A  444  SER ARG ALA PRO ALA PRO ALA THR PRO HIS ALA PRO ASP\n'

    expected_result = {
                        'record': 'SEQRES',
                        'cont': '1',
                        'chain': 'A',
                        'num_res': 444,
                        'resnames': ['SER', 'ARG', 'ALA', 'PRO', 'ALA', 'PRO', 
                                     'ALA', 'THR', 'PRO', 'HIS', 'ALA', 'PRO', 
                                     'ASP'],
                      }

    try:
        data, error = record_reader.parse_line(line, schema)
        
    finally:
        pass
        
    assert data == expected_result

def test_parse_ssbond_line():

    schema = record_reader.rec_schemas['SSBOND']

    line = 'SSBOND   1 CYS A  155    CYS A  359                          1555   1555  2.04  \n'

    expected_result = {
                        'record': 'SSBOND',
                        'serial': 1,
                        'resname1': 'CYS' ,
                        'chain1': 'A', 
                        'resid1': 155, 
                        'insert1': '', 
                        'resname2': 'CYS', 
                        'chain2': 'A', 
                        'resid2': 359, 
                        'insert2': '', 
                        'sym1': 1555, 
                        'sym2': 1555,
                      }

    try:
        data, error = record_reader.parse_line(line, schema)
        
    finally:
        pass

    assert data == expected_result

def test_parse_remark_line():
    
    schema = record_reader.rec_schemas['REMARK']
    
    line = 'REMARK 500  M RES CSSEQI        PSI       PHI\n'
    
    expected_result = {
                        'record': 'REMARK',
                        'num' : 500,
                        'text': ' M RES CSSEQI        PSI       PHI\n'
                      }
                      
    try:
        data, error = record_reader.parse_line(line, schema)
        
    finally:
        pass

    
    assert data == expected_result
    
if __name__ == '__main__': 
   main() 
    
    