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

    def test_is_obsolete_unset(self):
        
        self.assertEquals(self.h.is_obsolete(), False)

    def test_is_obsolete_set(self):

        self.h.pdb_recs['OBSLTE'] = [{'record':'dummy'}]
        self.assertEquals(self.h.is_obsolete(), True)

    def test_is_split_unset(self):
        
        self.assertEquals(self.h.is_split(), False)

    def test_is_split_set(self):

        self.h.pdb_recs['SPLIT'] = [{'record':'dummy'}]
        self.assertEquals(self.h.is_split(), True)
        
    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 