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
import sassie.interface.input_filter as input_filter
#import input_filter as input_filter



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
    __file__)), '..', '..', 'data', 'interface', 'input_filter') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Input_Filter(MockerTestCase):

    '''
    System integration test for extract_utilities_filter.py / sassie 1.0

    Test to see whether extract_utilities_filter catches improper input.

    Inputs tested:

    runname:                string          project name
    pdbfile                 string          input pdb file
    second_pdbfile          string          second input pdb file (compatible with input pdb file)
    dcdfile                 string          input dcd file (compatible with input pdb file)                                          
    complex_pdbfile:        string          input pdb file with more than one segment
    psffile                 string          input psf file (compatible with input pdbfile)   
    exefile                 string          executable file
    integer_number          integer         single integer number
    integer_array           integer array   array of integer numbers
    float_number            float           single floating number
    float_array             float array     array of floating numbers
    test_flag1              boolean         'False' string conversion to Boolean False
    test_flag2              boolean         'True' string conversion to Boolean True

    Inputs not tested:

    Strings other than 'True' or 'False' for conversion to Boolean were not tested
    since 'true', 'T', 't', 'TRUE', 'false', 'F', 'f', and 'FALSE' are deprecated in SASSIE 2.0

    Use cases tested:

    1.  type_check_and_convert
        a.  convert all variables:  no errors
        b.  float_number is not a floating number
        c.  integer_number is not an integer
        d.  float array can not be read
        e.  float array is not an array of floating numbers
        f.  integer array can not be read
        g.  integer array is not an array of integer numbers
        h.  flag string can be converted to boolean True or False 
    2.  check_pdb_dcd
        a.  pdb and dcd files exist and are valid:  no errors
        b.  pdb file does not exist
        c.  pdb file is not a valid pdb file
        d   dcd file does not exist
        e.  dcd file is not a valid dcd file
    3.  check_binary
        a.  file is binary      
        b.  file is not binary
    4.  check_file_exists
        a.  path/file does not exist    
        c.  path/file exists 
    5.  check_name
        a.  name has no incorrect characters
        b.  name has incorrect characters
    6.  check_exe
        a.  exe file exists and is executable
        b.  exe file does not exist                          
        c.  exe file is not accessible
        d   exe file is a directory                          
    7.  check_permissions
        a.  no path permission error
        b.  path has permission error
            i.   path does not exist
            ii.  read permission not allowed   NOTE: made a directory with no read permission on the fly, tested and then deleted
            iii. write permission not allowed
    8.  certify_pdb_pdb
        a.  both pdb files exist, can be read and are compatible (2nd pdb file is a trajectory file)
        b.  pdb file 1 does not exist
        c.  pdb file 2 does not exist
        d.  pdb file 1 can not be read
        e.  pdb file 2 can not be read
        f.  pdb files are not compatible
    9.  read_psf_file
        a.  test if the psf file can be read:  no errors    NOTE:  didn't test for errors since they are captured in the next two tests below
    10. certify_pdb_psf 
        a.  psf and pdb files exist, can be read and are compatible
        b.  psf file does not exist
        c.  dcd file does not exist
        d.  psf file can not be read
        e.  pdb file can not be read
        f.  psf and pdb files are not compatible
    11. certify_dcd_psf 
        a.  psf and dcd files exist, can be read and are compatible
        b.  psf file does not exist
        c.  dcd file does not exist                         NOT TESTED:  method assumes dcd file exists; segmentation fault occurs for non-existent file
        d.  psf file can not be read
        e.  dcd file can not be read
        f.  psf and dcd files are not compatible
    12. certify_pdb_dcd
        a.  pdb and dcd files exist, can be read and are compatible
        b.  pdb file does not exist
        c.  dcd file does not exist                         NOT TESTED:  method assumes dcd file exists; segmentation fault occurs for non-existent file
        d.  pdb file can not be read
        e.  dcd file can not be read
        f.  pdb and dcd files are not compatible    
    13. get_pdb_stats
        a.  pdb stats can be read 
        b.  pdb stats can not be read
    14. get_pdb_complex_stats
        a.  complex pdb stats can be read for protein and RNA segments
        b.  complex pdb stats can not be read for protein segment
    15. check_and_convert_formula
        a.  formula string can be converted
        b.  formula string can be converted but there is an incorrect element symbol
        c.  formula string can not be converted, i.e, string has bad character   
        d.  formula array cannot be read                         

    '''

    def setUp(self):

        '''
        defines variables that will be used to test the input filter
        '''

        pdb_data_path = paths['pdb_data_path']
        dcd_data_path = paths['dcd_data_path']
        other_data_path = paths['other_data_path']
        module_data_path = paths['module_data_path']

        self.runname         = 'run_0'
        self.pdbfile		     = os.path.join(pdb_data_path, 'hiv1_gag.pdb')
        self.second_pdbfile  = os.path.join(pdb_data_path, 'hiv1_gag_20_frames.pdb')
        self.complex_pdbfile = os.path.join(pdb_data_path, 'rna_protein_complex.pdb')
        self.psffile		     = os.path.join(other_data_path, 'hiv1_gag.psf')
        self.dcdfile         = os.path.join(dcd_data_path, 'hiv1_gag_20_frames.dcd')
        self.exefile  	     = os.path.join(other_data_path, 'xtal2sas.exe')
        self.integer_number  = '2'
        self.float_number    = '5.5'
        self.integer_array   = '1,2,3,4,5'
        self.float_array     = '1.5,2.5,3.5,3.5,5.5'
        self.formula_array   = ['KCl','C4H11NO3','(C42H82NO8P)10']
        self.test_flag1      = False
        self.test_flag2      = True
        self.precision = 3

    def run_filter(self, **kwargs):
        '''
        method to "run" the input filter to type check and convert the variables
        '''

        svariables = {}

        svariables['runname'] = (self.runname, 'string')
        svariables['pdbfile'] = (self.pdbfile, 'string')
        svariables['second_pdbfile'] = (self.second_pdbfile, 'string')
        svariables['complex_pdbfile'] = (self.complex_pdbfile, 'string')
        svariables['dcdfile'] = (self.dcdfile, 'string')
        svariables['psffile'] = (self.psffile, 'string')
        svariables['exefile'] = (self.exefile, 'string')
        svariables['integer_number'] = (self.integer_number, 'int')
        svariables['float_number'] = (self.float_number, 'float')
        svariables['integer_array'] = (self.integer_array, 'int_array')
        svariables['float_array'] = (self.float_array, 'float_array')
        svariables['formula_array'] = (self.formula_array, 'string')
        svariables['test_flag1'] = (self.test_flag1, 'boolean')       
        svariables['test_flag2'] = (self.test_flag2, 'boolean')       
           

        error, self.variables = input_filter.type_check_and_convert(svariables)

        if(len(error) > 0):
            print 'error = ', error
#           sys.exit()
            return error


    def extract_important_path0(self, return_error):

        string_error = string.split(return_error[0])
#        print 'string_error: ', string_error
        path_list = string.split(string_error[-4], '..')
#        print 'path_list: ', path_list
        list1 = string_error[:2]
#        print 'list1: ', list1
        list2 = string_error[-3:]
#        print 'list2: ',list2
        important_path = string.split(path_list[-1], "/")[1:]
#        print 'important_path: ', important_path
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
#            print 'error: ', error
        error = ' '.join(list1)+' : '+error+' '+' '.join(list2)
#        print 'final error: ', error    
        return error

    def extract_important_path1(self, return_error):

        string_error = string.split(return_error[0])
#        print 'string_error: ', string_error
        path_list = string.split(string_error[-5], '..')
#        print 'path_list: ', path_list
        list1 = string_error[:2]
#        print 'list1: ', list1
        list2 = string_error[-4:]
#        print 'list2: ',list2
        important_path = string.split(path_list[-1], "/")[1:]
#        print 'important_path: ', important_path
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
#            print 'error: ', error
        error = ' '.join(list1)+' : '+error+' '+' '.join(list2)
#        print 'final error: ', error    
        return error

    def assert_list_almost_equal(self, a, b, places=5):
        if (len(a) != len(b)):
            raise TypeError
        else:
            for i in range(len(a)):
                if isinstance(a[i], (int, float, numpy.generic)):
                    if (numpy.isnan(a[i]) and numpy.isnan(b[i])):
                        continue
                    self.assertAlmostEqual(a[i], b[i], places)
                else:
                    self.assert_list_almost_equal(a[i], b[i], places)



    def test_27(self):
        '''
        check if path has write permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder1')
        ''' see if you can write to the directory '''
        print os.access('empty_folder1', os.W_OK)
        ''' make the directory un-writable'''
        os.system('chmod a-w empty_folder1')
        ''' see if you can write to the directory '''
        print os.access('empty_folder', os.W_OK)

        self.pdbfile = os.path.join('./','empty_folder1')
        self.run_filter()
            
        existvalue, readvalue, writevalue = input_filter.check_permissions(self.pdbfile)
        print 'existvalue, readvalue, writevalue: ',existvalue, readvalue, writevalue

        ''' check for file error '''
        expected_writevalue = False
        assert_equals(writevalue, expected_writevalue)

        ''' make the directory writable'''
        os.system('chmod a+w empty_folder1')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder1')
                                         



    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
