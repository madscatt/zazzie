'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import os
import sassie.tools.gui_mimic_align as gui_mimic_align
#import gui_mimic_align as gui_mimic_align

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'align') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Align_Filter(MockerTestCase):

    '''
    System integration test for align_filter.py / sassie 1.0

    Test to see whether align_filter catches improper input.

    Inputs tested:

    runname:      string      project name                          
    path:         string      input/output filepath                 
    pdbmol1:      string      reference pdb (mol 1)                 
    pdbmol2:      string      input pdb file (mol 2)                
    infile:       string      input (pdb or dcd) filename (mol 2)   
    basis1:       string      basis for molecule 1                 
    basis2:       string      basis for molecule 2                  
    lowres1:      integer     low residue for overlap molecule 1    
    highres1:     integer     high residue for overlap molecule 1   
    lowres2:      integer     low residue for overlap molecule 2    
    highres2:     integer     high residue for overlap molecule 2   

    Inputs not tested (deprecated functionality):

    ebasis1:      string      extra basis statement molecule 1 
    ebasis2:      string      extra basis statement molecule 2 

    Use cases tested:

    1.  check pdbmol1
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    2.  check pdbmol2
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    3.  check infile
        a.  infile is .pdb (value = 0)
        b.  infile is .dcd (value = 1)
        c.  infile doesn't have .pdb or .dcd suffix
    4.  check lowres1 and highres1
        q.  pdbmol1 contains lowres1 and/or highres1 
        b.  pdbmol1 doesn't contain lowres1 and/or highres1
    5.  check lowres2 and highres2
        a.  pdbmol2 contains lowres2 and/or highres2 
        b.  pdbmol2 doesn't contain lowres2 and/or highres2
    6.  check size of alignment basis
        a.  alignment basis is too small
    7.  check basis atoms
        a.  mol1 and mol2 have a different number of basis atoms
        b.  basis1 and basis2 don't match  
    8.  check input file path permissions 
        a.  no permission error
        b. permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed     
    '''

    def setUp(self):

        gui_mimic_align.test_variables(self, paths)


    def test_1(self):
        '''
        test if pdbmol1 exists
        '''
        self.pdbmol1 = os.path.join(
            module_data_path, 'does_not_exist!&@#X.pdb')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['mol 1 reference pdb file, ' +
                          self.pdbmol1 + ', does not exist']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if pdbmol1 is a valid pdb file
        '''
        self.pdbmol1 = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['mol 1 reference pdb file, ' +
                          self.pdbmol1 + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if pdbmol2 exists
        '''
        self.pdbmol2 = os.path.join(
            module_data_path, 'does_not_exist!&@#XXX.pdb')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['mol 2 reference pdb file, ' +
                          self.pdbmol2 + ', does not exist']
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test if pdbmol2 is a valid pdb file
        '''
        self.pdbmol2 = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['mol 2 reference pdb file, ' +
                          self.pdbmol2 + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if "pdb" infile exists
        '''
        self.infile = os.path.join(
            module_data_path, 'does_not_exist!&@#XXX.pdb')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input (pdb or dcd) file ' +
                          self.infile + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if "dcd" infile exists
        '''
        self.infile = os.path.join(
            module_data_path, 'does_not_exist!&@#XXX.dcd')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input (pdb or dcd) file ' +
                          self.infile + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if infile is a has a "pdb" or "dcd" suffix
        '''
        self.infile = os.path.join(module_data_path, 'wrong_suffix.txt')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file suffix error '''
        expected_error = ['infile needs to have a "pdb" or "dcd" suffix']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if "dcd" infile is a valid dcd file
        '''
        self.infile = os.path.join(module_data_path, 'not_valid.dcd')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['mol 2 input file, ' +
                          self.infile + ', is not a valid pdb or dcd file']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if "dcd" infile has the same number of atoms as pdbmol2
        '''
        self.infile = os.path.join(module_data_path, 'non_matching.dcd')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['mol 2 pdbfile ' + self.pdbmol2 + ' and dcdfile ' +
                          self.infile + ' are not compatible (different number of atoms)']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if "pdb" infile has the same number of atoms as pdbmol2
        '''
        self.infile = os.path.join(module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_align.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['mol 2 reference pdb file ' + self.pdbmol2 +
                          ' and input pdb file ' + self.infile + ' are not compatible']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if lowres1 is in pdbmol1
        '''
        self.lowres1 = '500'
        return_error = gui_mimic_align.run_module(
            self, test_filter=True)
        ''' check for file error '''
        expected_error = ['mol 1 pdb file does not have low residue amino acid, 500, range = 1 : 431']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if highres1 is in pdbmol1
        '''
        self.highres1 = '500'
        return_error = gui_mimic_align.run_module(
            self, test_filter=True)
        ''' check for file error '''
        expected_error = ['mol 1 pdb file does not have high residue amino acid, 500, range = 1 : 431']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if lowres2 is in pdbmol2
        '''
        self.lowres2 = '500'
        return_error = gui_mimic_align.run_module(
            self, test_filter=True)
        ''' check for file error '''
        expected_error = ['mol 2 pdb file does not have low residue amino acid, 500, range = 1 : 431']
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test if highres2 is in pdbmol2
        '''
        self.highres2 = '500'
        return_error = gui_mimic_align.run_module(
            self, test_filter=True)
        ''' check for file error '''
        expected_error = ['mol 2 pdb file does not have high residue amino acid, 500, range = 1 : 431']
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test if pdbmol1 alignment basis is too small
        '''
        self.lowres1 = '349'
        return_error = gui_mimic_align.run_module(self, test_filter=True)
        ''' check for file error '''
        expected_error = [
            'mol 1 alignment basis is too small (less than 3 points) or low residue > high residue']
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if pdbmol2 alignment basis is too small
        '''
        self.lowres2 = '349'
        return_error = gui_mimic_align.run_module(self, test_filter=True)
        ''' check for file error '''
        expected_error = [
            'mol 2 alignment basis is too small (less than 3 points) or low residue > high residue']
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test if mol1 and mol2 have a different number of basis atoms
        '''
        self.lowres1 = '144'
        return_error = gui_mimic_align.run_module(self, test_filter=True)
        ''' check for file error '''
        expected_error = [
            'mol 1 and mol2 have a different number of basis atoms, check low and high residue input values']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if basis1 and basis2 match
        '''
        self.basis1 = 'P'
        return_error = gui_mimic_align.run_module(self, test_filter=True)
        ''' check for file error '''
        expected_error = ['basis1 = ' + self.basis1 +
                          ' and basis2 = ' + self.basis2 + ' do not match']
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test if path exists
        '''
        self.inpath = os.path.join(module_data_path,'non_existent_path')
        return_error = gui_mimic_align.run_module(self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.inpath + '  [code = FalseFalseFalse]',
   'path does not exist']
        assert_equals(return_error, expected_error)    
        
    def test_20(self):
        '''
        test if directory has read permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder')
        ''' see if you can read the directory '''
        print os.access('empty_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_folder')
        ''' see if you can read the directory '''
        print os.access('empty_folder', os.R_OK)

        self.inpath= os.path.join('./','empty_folder')
        return_error = gui_mimic_align.run_module(self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.inpath + '  [code = TrueFalseTrue]', 'read permission not allowed']
        assert_equals(return_error, expected_error)  

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')

    def test_21(self):
        '''
        test if directory has write permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder1')
        ''' see if you can write to the directory '''
#        print os.access('empty_folder1', os.W_OK)
        ''' make the directory un-writeable'''
        os.system('chmod a-w empty_folder1')
        ''' see if you can write to the directory '''
#        print os.access('empty_folder', os.W_OK)

        self.inpath = os.path.join('./', 'empty_folder1')
        return_error = gui_mimic_align.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.inpath + '  [code = TrueTrueFalse]', 'write permission not allowed']
        assert_equals(return_error, expected_error)

        ''' make the directory writeable'''
        os.system('chmod a+w empty_folder1')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder1')        


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
