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
import string
import sassie.analyze.gui_mimic_apbs as gui_mimic_apbs
#import gui_mimic_apbs as gui_mimic_apbs

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'interface', 'apbs') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Apbs_Filter(MockerTestCase):

    '''
    System integration test for apbs_filter.py / sassie 1.0

    Test to see whether apbs_filter catches improper input.

    Inputs tested:

            runname:        project name
            pdbfile:        reference pdb name
            infile:         input trajectory filename (pdb or dcd)
            ph:             pH value
            temperature:    temperature value (K)
            ion_conc:       concentration of solute ion (M)
            ion_radius:     radius of solute ion (angstroms)            

    Inputs not tested:  (Not yet implemented)

            manual_flag:    user input file flag (0=default inputs; 1=user input file)
            manual_file:    user input file name


    Use cases tested:

    1.  check runname
        a.  check for invalid characters in runname    
    2.  check input PDB file
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    3.  check input trajectory file
        a.  file doesn't exist
        b.  file exists
            i.  file name doesn't start with "." and must end in ".pdb" or ".dcd"
            ii. file is .pdb (value = 0)
                1.  PDB file is valid
                2.  PDB file isn't valid
            iii.file is .dcd (value = 1)
                1.  DCD file is valid
                2.  DCD file isn't valid
    4.  check if input PDB and trajectory files are compatible
        a.  trajectory file is a DCD file
            i.  input PDB and DCD are compatible
            ii. input PDB and DCD are not compatible
        b.  trajectory file is a PDB file
            i.  input PDB and PDB are compatible
            ii. input PDB and PDB are not compatible
    5.  check if temperature >0
    6.  check if pH >=0
    7.  check if ion concentration >=0
    8.  check if ion radius > 0

    '''

    def setUp(self):

        gui_mimic_apbs.test_variables(self, paths)


    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if input PDB file exists
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['reference pdb file, ' +
                          self.pdbfile[3:] + ', does not exist']                 
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if input PDB file is a valid file
        '''
        self.pdbfile = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['reference pdb file, ' +
                          self.pdbfile[3:] + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test if input trajectory file exists
        '''
        self.infile = os.path.join(
            module_data_path, 'does_not_exist.dcd')
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' +
                          self.infile[3:] + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if input trajectory file is a valid dcd file
        '''
        self.infile = os.path.join(
            module_data_path, 'not_valid.dcd')
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file, ' +
                          self.infile[3:] + ', is not a valid pdb or dcd file']
        assert_equals(return_error, expected_error)


    def test_6(self):
        '''
        test if input trajectory file is a valid pdb file
        '''
        self.infile = os.path.join(
            module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file, ' +
                          self.infile[3:] + ', is not a valid pdb or dcd file']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if "dcd" trajectory file has the same number of atoms as input PDB
        '''
        self.infile = os.path.join(
            module_data_path, 'non_matching.dcd')
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb ' + self.pdbfile[3:] + ' and dcd file ' +
                          self.infile[3:] + ' are not compatible']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if "pdb" trajectory file has the same number of atoms as input PDB
        '''
        self.infile = os.path.join(
            module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['reference pdb file ' + self.pdbfile[3:] + ' and input pdb file ' +
                          self.infile[3:] + ' are not compatible']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if temperature is > 0
        '''
        self.temperature = '0'
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['temperature needs to be greater than zero : 0.0']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if pH >=0
        '''
        self.ph = '-1'
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ph needs to be greater than or equal to zero : -1.0']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if ion concentration >=0
        '''
        self.ion_conc = '-1'
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ion concentration needs to be greater than or equal to zero : -1.0']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if ion radius >0
        '''
        self.ion_radius = '0'
        return_error = gui_mimic_apbs.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ion radius needs to be greater than zero : 0.0']
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)
            
if __name__ == '__main__':
    main()
