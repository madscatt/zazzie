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
        self.test_flag1      = 'False'
        self.test_flag2      = 'True'
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

    '''Testing type_check_and_convert'''

    def test_1(self):
        '''
        test type check and convert variables (string, integer, integer array, float, float array, boolean):  no errors
        '''

        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = None
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if float_number is a floating number
        '''

        self.float_number = 'xy'
        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = ['float_number is not a float : xy']
        assert_equals(return_error, expected_error)
        

    def test_3(self):
        '''
        test if integer_number is an integer
        '''

        self.integer_number = 'xy'
        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = ['integer_number is not a int : xy']
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test if float array is an array of floating numbers
        '''

        self.float_array = '1.5,x,6.0'
        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = ['float_array is not a float_array : 1.5,x,6.0']
        assert_equals(return_error, expected_error)
                
    def test_5(self):
        '''
        test if float array can be read
        '''

        self.float_array = 10
        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = ['float_array: could not read array of values']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if integer array is an array of integer numbers
        '''

        self.integer_array = '1,2.5,6'
        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = ['integer_array is not a int_array : 1,2.5,6']
        assert_equals(return_error, expected_error)
                
    def test_7(self):
        '''
        test if integer array can be read
        '''

        self.integer_array = 10
        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = ['integer_array: could not read array of values']
        assert_equals(return_error, expected_error) 

    def test_8(self):
        '''
        test if string can be converted to boolean True or False
        '''

        self.test_flag1 = 'Y'
        return_error = self.run_filter()

        ''' check for variable error '''
        expected_error = ['test_flag1: could not boolean input type']
        assert_equals(return_error, expected_error)

    ''' Testing check_pdb_dcd'''

    def test_9(self):
        '''
        test if PDB file exists and is valid:  no errors
        '''
        
        self.run_filter()
            
        fileexist,value = input_filter.check_pdb_dcd(self.pdbfile,'pdb')

        ''' check fileexist and value flags '''
        expected_fileexist = 1
        assert_equals(fileexist, expected_fileexist)
        expected_value = 1
        assert_equals(value, expected_value)

    def test_10(self):
        '''
        test if PDB file exists
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        self.run_filter()
            
        fileexist,value = input_filter.check_pdb_dcd(self.pdbfile,'pdb')

        ''' check fileexist flag '''
        expected_fileexist = False
        assert_equals(fileexist, expected_fileexist)

    def test_11(self):
        '''
        test if PDB file is a valid PDB file
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'not_valid.pdb')
        self.run_filter()
            
        fileexist,value = input_filter.check_pdb_dcd(self.pdbfile,'pdb')

        ''' check if value = 0 '''
        expected_value = 0
        assert_equals(value, expected_value)

    def test_12(self):
        '''
        test if DCD file exists and is valid:  no errors
        '''
        
        self.run_filter()
            
        fileexist,value = input_filter.check_pdb_dcd(self.dcdfile,'dcd')

        ''' check fileexist and value flags '''
        expected_fileexist = 1
        assert_equals(fileexist, expected_fileexist)
        expected_value = 1
        assert_equals(value, expected_value)

    def test_13(self):
        '''
        test if DCD file exists
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'does_not_exist.dcd')
        self.run_filter()
            
        fileexist,value = input_filter.check_pdb_dcd(self.dcdfile,'dcd')

        ''' check fileexist flag '''
        expected_fileexist = False
        assert_equals(fileexist, expected_fileexist)

    def test_14(self):
        '''
        test if DCD file is a valid DCD file
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'not_valid.dcd')
        self.run_filter()
            
        fileexist,value = input_filter.check_pdb_dcd(self.dcdfile,'dcd')

        ''' check if value = 0 '''
        expected_value = 0
        assert_equals(value, expected_value)

    '''Testing check_binary'''

    def test_15(self):
        '''
        check if DCD file is binary
        '''

        self.run_filter()
            
        flag = input_filter.check_binary(self.dcdfile)

        ''' check binary flag '''
        expected_flag = True
        assert_equals(flag, expected_flag)

    def test_16(self):
        '''
        check if PDB file not binary
        '''

        self.run_filter()
            
        flag = input_filter.check_binary(self.pdbfile)

        ''' check binary flag '''
        expected_flag = False
        assert_equals(flag, expected_flag)
    
    ''' Testing check_file_exists'''
    
    def test_17(self):
        '''
        check if path/file doesn't exist
        '''
            
        self.pdbfile = os.path.join('./','non_existent.pdb')
        self.run_filter()
            
        return_error = input_filter.check_file_exists(self.pdbfile)

        ''' check for file error '''
        expected_error = ['file : ./non_existent.pdb does not exist']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        check if path/file exists
        '''

        self.run_filter()
            
        return_error = input_filter.check_file_exists(self.pdbfile)

        ''' check for variable error '''
        expected_error = []
        assert_equals(return_error, expected_error)

    '''Testing check_name'''

    def test_19(self):
        '''
        check if name has no incorrect characters
        '''

        self.run_filter()
            
        return_error = input_filter.check_name(self.pdbfile)

        ''' check for file error '''
        expected_error = []
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        check if name has incorrect characters
        '''

        self.run_filter()
            
        return_error = input_filter.check_name('bad_&.pdb')

        ''' check for file error '''
        expected_error = ['file or path : bad_&.pdb has incorrect character : &']
        assert_equals(return_error, expected_error)

    '''Testing check_exe'''

    def test_21(self):
        '''
        check if exe file exists and is executable:  no errors
        '''

        self.run_filter()
            
        return_error = input_filter.check_exe(self.exefile)

        ''' check for file error '''
        expected_error = []
        assert_equals(return_error, expected_error)

    def test_22a(self):
        '''
        check if exe file is not accessible
        '''

        self.exefile = os.path.join(module_data_path,'not_executable.exe')
        self.run_filter()
            
        return_error = input_filter.check_exe(self.exefile)
#        print 'return_error: ',return_error
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path0(return_error)
#        print 'relative_path: ', relative_path
        new_error = [relative_path]
#        print 'new error: ', new_error

        ''' check for file error '''
        expected_error = ['Executable file : ../../data/interface/input_filter/not_executable.exe is not accessible']
        assert_equals(new_error, expected_error)


    def test_22b(self):
        '''
        check if exe file is not a file
        '''

        self.exefile = os.path.join(module_data_path,'non_existent.exe')
        self.run_filter()
            
        return_error = input_filter.check_exe(self.exefile)
#        print 'return_error: ',return_error
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path1(return_error)
#        print 'relative_path: ', relative_path
        new_error = [relative_path]
#        print 'new error: ', new_error

        ''' check for file error '''
        expected_error = ['Executable file : ../../data/interface/input_filter/non_existent.exe is not a file']
        assert_equals(new_error, expected_error)
                                         

    def test_22c(self):
        '''
        check if exe file is a directory
        '''

        self.exefile = os.path.join(module_data_path,'is_a_directory')
        self.run_filter()
            
        return_error = input_filter.check_exe(self.exefile)
#        print 'return_error: ',return_error
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path0(return_error)
#        print 'relative_path: ', relative_path
        new_error = [relative_path]
#        print 'new error: ', new_error

        ''' check for file error '''
        expected_error = ['Executable file : ../../data/interface/input_filter/is_a_directory is a directory!']
        assert_equals(new_error, expected_error)
                                         

                                                                   
    '''Testing check_permissions'''

    def test_23(self):
        '''
        check if path has exists and has read and write permissions:  no errors
        '''

        self.run_filter()
            
        existvalue, readvalue, writevalue = input_filter.check_permissions(self.pdbfile)
#        print 'existvalue, readvalue, writevalue: ',existvalue, readvalue, writevalue

        ''' check for file error '''
        expected_existvalue = True
        expected_readvalue = True
        expected_writevalue = True
        assert_equals(existvalue, expected_existvalue)
        assert_equals(readvalue, expected_readvalue)
        assert_equals(writevalue, expected_writevalue)

    def test_24(self):
        '''
        check if path has write permission
        '''
        self.pdbfile = os.path.join(module_data_path,'non_existent_path')
        self.run_filter()
            
        existvalue, readvalue, writevalue = input_filter.check_permissions(self.pdbfile)
#        print 'existvalue, readvalue, writevalue: ',existvalue, readvalue, writevalue

        ''' check for file error '''
        expected_existvalue = False
        assert_equals(existvalue, expected_existvalue)
                                         
    def test_25(self):
        '''
        check if path has read permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)


        self.pdbfile = os.path.join('./','empty_folder')
        self.run_filter()
            
        existvalue, readvalue, writevalue = input_filter.check_permissions(self.pdbfile)
#        print 'existvalue, readvalue, writevalue: ',existvalue, readvalue, writevalue

        ''' check for file error '''
        expected_readvalue = False
        assert_equals(readvalue, expected_readvalue)
                                         
        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')

    def test_26(self):
        '''
        check if path has write permission
        '''
        self.pdbfile = os.path.join(module_data_path,'no_write_permission')
        self.run_filter()
            
        existvalue, readvalue, writevalue = input_filter.check_permissions(self.pdbfile)
#        print 'existvalue, readvalue, writevalue: ',existvalue, readvalue, writevalue

        ''' check for file error '''
        expected_writevalue = False
        assert_equals(writevalue, expected_writevalue)
                                         
    '''Testing certify_pdb_pdb'''

    def test_27(self):
        '''
        check if pdb files exist, can be read and are compatible:  no errors
        '''

        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_pdb(self.pdbfile,self.second_pdbfile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 1
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_28(self):
        '''
        check if pdb file exists
        '''

        self.pdbfile = os.path.join(module_data_path,'non_existent.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_pdb(self.pdbfile,self.second_pdbfile)

        ''' check for file error '''
        expected_filexist = 0
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)
                                          
    def test_29(self):
        '''
        check if second pdb file exists
        '''

        self.second_pdbfile = os.path.join(module_data_path,'non_existent.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_pdb(self.pdbfile,self.second_pdbfile)

        ''' check for file error '''
        expected_filexist = 0
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_30(self):
        '''
        check if pdb file can be read
        '''

        self.pdbfile = os.path.join(module_data_path,'not_valid.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_pdb(self.pdbfile,self.second_pdbfile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)
                                          
    def test_31(self):
        '''
        check if second pdb file can be read
        '''

        self.second_pdbfile = os.path.join(module_data_path,'not_valid.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_pdb(self.pdbfile,self.second_pdbfile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_32(self):
        '''
        check if pdb files are compatible
        '''

        self.second_pdbfile = os.path.join(module_data_path,'non_matching.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_pdb(self.pdbfile,self.second_pdbfile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)
                                          
    '''Testing certify_pdb_psf'''

    def test_33(self):
        '''
        check if pdb and psf files exist, can be read and are compatible:  no errors
        '''

        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_psf(self.pdbfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 1
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_34(self):
        '''
        check if pdb file exists
        '''

        self.pdbfile = os.path.join(module_data_path,'non_existent.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_psf(self.pdbfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)
                                          
    def test_35(self):
        '''
        check if psf file exists
        '''

        self.psffile = os.path.join(module_data_path,'non_existent.psf')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_psf(self.pdbfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 0
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_36(self):
        '''
        check if pdb file can be read
        '''

        self.pdbfile = os.path.join(module_data_path,'not_valid.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_psf(self.pdbfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)
                                          
    def test_37(self):
        '''
        check if psf file can be read
        '''

        self.psffile = os.path.join(module_data_path,'not_valid.psf')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_psf(self.pdbfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_38(self):
        '''
        check if pdb and psf files are compatible
        '''

        self.pdbfile = os.path.join(module_data_path,'non_matching.pdb')
        self.run_filter()
            
        fileexist, value = input_filter.certify_pdb_psf(self.pdbfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    '''Testing certify_dcd_psf''' 

    def test_39(self):
        '''
        check if dcd and psf files exist, can be read and are compatible:  no errors
        '''

        self.run_filter()
            
        fileexist, value = input_filter.certify_dcd_psf(self.dcdfile,self.psffile)       

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 1
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)                                             

    def test_40(self):
        '''
        check if psf file exists
        '''

        self.psffile = os.path.join(module_data_path,'non_existent.psf')
        self.run_filter()
            
        fileexist, value = input_filter.certify_dcd_psf(self.dcdfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 0
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_41(self):
        '''
        check if dcd file can be read
        '''

        self.dcdfile = os.path.join(module_data_path,'not_valid.dcd')
        self.run_filter()
            
        fileexist, value = input_filter.certify_dcd_psf(self.dcdfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)
                                          
    def test_42(self):
        '''
        check if psf file can be read
        '''

        self.psffile = os.path.join(module_data_path,'not_valid.psf')
        self.run_filter()
            
        fileexist, value = input_filter.certify_dcd_psf(self.dcdfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)

    def test_43(self):
        '''
        check if dcd and psf files are compatible
        '''

        self.dcdfile = os.path.join(module_data_path,'non_matching.dcd')
        self.run_filter()
            
        fileexist, value = input_filter.certify_dcd_psf(self.dcdfile,self.psffile)

        ''' check for file error '''
        expected_filexist = 1
        expected_value = 0
        assert_equals(fileexist, expected_filexist)
        assert_equals(value, expected_value)
                                          
    '''Testing read_psf_file'''

    def test_44(self):
        '''
        check if psf file can be read:  no errors
        '''

        self.run_filter()
            
        natoms, names = input_filter.read_psf_file(self.psffile)

        ''' check for file error '''
        expected_natoms = 6730
        assert_equals(natoms, expected_natoms)

    '''Testing certify_pdb_dcd'''


    def test_45(self):
        '''
        check if pdb and dcd files can be read and are compatible:  no errors
        '''

        self.run_filter()
            
        value = input_filter.certify_pdb_dcd(self.pdbfile,self.dcdfile)       

        ''' check for file error '''
        expected_value = 1
        assert_equals(value, expected_value)                                             

    def test_46(self):
        '''
        check if pdb file exists
        '''

        self.pdbfile = os.path.join(module_data_path,'non_existent.pdb')
        self.run_filter()
            
        value = input_filter.certify_pdb_dcd(self.pdbfile,self.dcdfile)       

        ''' check for file error '''
        expected_value = 0
        assert_equals(value, expected_value)

    def test_47(self):
        '''
        check if dcd file can be read
        '''

        self.dcdfile = os.path.join(module_data_path,'not_valid.dcd')
        self.run_filter()
            
        value = input_filter.certify_pdb_dcd(self.pdbfile,self.dcdfile)       

        ''' check for file error '''
        expected_value = 0
        assert_equals(value, expected_value)
                                          
    def test_48(self):
        '''
        check if pdb file can be read
        '''

        self.pdbfile = os.path.join(module_data_path,'not_valid.pdb')
        self.run_filter()
            
        value = input_filter.certify_pdb_dcd(self.pdbfile,self.dcdfile)      

        ''' check for file error '''
        expected_value = 0
        assert_equals(value, expected_value)

    def test_49(self):
        '''
        check if pdb and dcd files are compatible
        '''

        self.dcdfile = os.path.join(module_data_path,'non_matching.dcd')
        self.run_filter()
            
        value = input_filter.certify_pdb_dcd(self.pdbfile,self.dcdfile)       

        ''' check for file error '''
        expected_value = 0
        assert_equals(value, expected_value)
                                          
    '''Testing get_pdb_stats'''

    def test_50(self):
        '''
        check if pdb file stats can be read:  no errors
        '''

        self.run_filter()

        locvariables = ['name','resid']
            
        value, result = input_filter.get_pdb_stats(self.pdbfile,locvariables)
        result00 = result[0][0]
        result10 = result[1][0]

        ''' check for file error '''
        expected_value = 1
        expected_result00 = 'N'
        expected_result10 = 1
        assert_equals(value, expected_value)
        assert_equals(result00, expected_result00)
        assert_equals(result10, expected_result10)

    def test_51(self):
        '''
        check if pdb file stats can be read
        '''

        self.pdbfile = os.path.join(module_data_path,'not_valid.pdb')
        self.run_filter()

        locvariables = ['name','resid']
            
        value, result = input_filter.get_pdb_stats(self.pdbfile,locvariables)

        ''' check for file error '''
        expected_value = 0
        expected_result = None
        assert_equals(value, expected_value)
        assert_equals(result, expected_result)

    '''Testing get_pdb_complex_stats'''

    def test_52(self):
        '''
        check if complex pdb file protein segment stats can be read:  no errors
        '''

        self.run_filter()

        segname = 'HFQ1'
        locvariables = ['resid','moltype']
            
        value, result = input_filter.get_pdb_complex_stats(self.complex_pdbfile,segname,locvariables)
        print 'result: ', result[0][0], result[1][0]
        result00 = result[0][0]
        result10 = result[1][0]

        ''' check for file error '''
        expected_value = 1
        expected_result00 = 1
        expected_result10 = 'protein'
        assert_equals(value, expected_value)
        assert_equals(result00, expected_result00)
        assert_equals(result10, expected_result10)

    def test_52(self):
        '''
        check if complex pdb file RNA segment stats can be read:  no errors
        '''

        self.run_filter()

        segname = 'RNA1'
        locvariables = ['resid','moltype']
            
        value, result = input_filter.get_pdb_complex_stats(self.complex_pdbfile,segname,locvariables)
        result00 = result[0][0]
        result10 = result[1][0]

        ''' check for file error '''
        expected_value = 1
        expected_result00 = 1
        expected_result10 = 'rna'
        assert_equals(value, expected_value)
        assert_equals(result00, expected_result00)
        assert_equals(result10, expected_result10)        

    def test_54(self):
        '''
        check if complex pdb file protein segment stats can be read
        '''

        self.complex_pdbfile = os.path.join(module_data_path,'not_valid.pdb')
        self.run_filter()

        segname = 'TEST'
        locvariables = ['resid','moltype']
            
        value, result = input_filter.get_pdb_complex_stats(self.complex_pdbfile,segname,locvariables)

        ''' check for file error '''
        expected_value = 0
        expected_result = None
        assert_equals(value, expected_value)
        assert_equals(result, expected_result)
                                          
    '''Testing check_and_convert_formula'''

    def test_55(self):
        '''
        check if formulas can be read and converted:  no errors
        '''

        self.run_filter()
                    
        error, formulas = input_filter.check_and_convert_formula(self.formula_array)

        ''' check for file error '''
        expected_error = []
        expected_formulas = [{'K': 1, 'Cl': 1}, {'H': 11, 'C': 4, 'O': 3, 'N': 1}, {'H': 820, 'C': 420, 'P': 10, 'O': 80, 'N': 10}]
        assert_equals(error, expected_error)
        assert_equals(formulas, expected_formulas)

    def test_56(self):
        '''
        check if formulas has non-alphanumeric character
        '''

        self.formula_array = ['N@Cl']
        self.run_filter()
            
        error, formulas = input_filter.check_and_convert_formula(self.formula_array)


        ''' check for file error '''

        expected_formulas = []
        assert_equals(formulas, expected_formulas)
        expected_error = [ValueError('unexpected character:\nN@Cl\n ^\n',)]
        print 'expected_error: ', expected_error
        print expected_error[0].message
        assert_equals(error[0].message, expected_error[0].message)

    def test_57(self):
        '''
        check if formulas has an incorrect character
        '''

        self.formula_array = ['Nacl']
        self.run_filter()
            
        error, formulas = input_filter.check_and_convert_formula(self.formula_array)


        ''' check for file error '''
        print 'error: ', error
        print error[0].message

        expected_formulas = []
        assert_equals(formulas, expected_formulas)
        expected_error = [ValueError("'Nacl' is not an element symbol:\nNacl\n^\n",)]
        print 'expected_error: ', expected_error
        print expected_error[0].message
        assert_equals(error[0].message, expected_error[0].message)

    def test_58(self):
        '''
        check if formula array length can be found
        '''

        self.formula_array=None
#        self.run_filter()

        error, formulas = input_filter.check_and_convert_formula(self.formula_array)


        ''' check for file error '''

        expected_formulas = []
        assert_equals(formulas, expected_formulas)
        expected_error = ['unable to read formula']
        assert_equals(error, expected_error)




    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
