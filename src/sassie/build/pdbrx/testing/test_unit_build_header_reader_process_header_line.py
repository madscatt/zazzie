# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 13:31:00 2015

@author: dave
"""

from unittest import main 

from mocker import Mocker, MockerTestCase
import warnings; warnings.filterwarnings('ignore')

import pdbscan.header_reader as header_reader

class Test_unit_build_PdbHeader_process_header_line(MockerTestCase):
    
    def setUp(self):
        self.h=header_reader.PdbHeader()

    def test_process_header_line_seqres_ok(self):
        
        line = 'SEQRES   1 A  444  SER ARG ALA PRO ALA PRO ALA THR PRO HIS ALA PRO ASP\n'

        expected_result = [{
                          'record': 'SEQRES',
                          'cont': '1',
                          'chain': 'A',
                          'num_res': 444,
                          'resnames': ['SER', 'ARG', 'ALA', 'PRO', 'ALA', 'PRO', 
                                       'ALA', 'THR', 'PRO', 'HIS', 'ALA', 'PRO', 
                                       'ASP'],
                          }]
                          
        self.h.process_header_line(line)

        self.assertEqual(self.h.pdb_recs['SEQRES'], expected_result)     

    def test_process_header_line_seqres_notok(self):
        
        line = 'SEQRES   1 A  BBB  SER ARG ALA PRO ALA PRO ALA THR PRO HIS ALA PRO ASP\n'
        
        with self.assertRaises(IOError):                   
            self.h.process_header_line(line)

    def test_process_header_line_nonstd_ok(self):
        
        line = 'NOTOK1   1 A  444  SER ARG ALA PRO ALA PRO ALA THR PRO HIS ALA PRO ASP\n'

        expected_result = [{
                          'record': 'NOTOK1',
                          'text': '1 A  444  SER ARG ALA PRO ALA PRO ALA THR PRO HIS ALA PRO ASP',
                          }]
                          
        self.h.process_header_line(line)

        self.assertEqual(self.h.pdb_recs['NONSTD'], expected_result)  

    def test_process_header_line_scale_ok(self):
        
        line = 'SCALE1      0.010601  0.006120  0.000000        0.00000\n'

        expected_result = [{
                           'record': 'SCALE1',
                           'sn1': 0.010601,
                           'sn2': 0.006120,
                           'sn3': 0.000000,
                           'un': 0.00000,
                          }]
        
        self.h.process_header_line(line)

        self.assertEqual(self.h.pdb_recs['SCALEn'], expected_result) 

        
    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 