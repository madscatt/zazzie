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
import sassie.tools.extract_utilities.gui_mimic_extract_utilities as gui_mimic_extract_utilities
#import gui_mimic_extract_utilities as gui_mimic_extract_utilities

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
    __file__)), '..', '..', 'data', 'interface', 'extract_utilities') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Extract_Utilities_Filter(MockerTestCase):

    '''
    System integration test for extract_utilities_filter.py / sassie 1.0

    Test to see whether extract_utilities_filter catches improper input.

    Inputs tested:

    runname:                string      project name
    pdb_filename            string      input pdb file
    trajectory_filename     string      input dcd file                          
    path:                   string      input/output filepath                 
    option:                 string      extract option (single_frame, range, all, weight_file, text_file)                 
    local_value:            string      value of option (frame value, range, all, weight_file name, text_file name)                
    output_filename:        string      output dcd file   
    extract_trajectory:     boolean     extract frames from trajectory (True or False)                
    extract_sas :           boolean     extract corresponding SAS files (True or False)                 
    sas_type:               integer     integer depends on SAS file type (SasCalc, Xtal2sas, Cryson, Crysol)    
    sas_paths:              string      path(s) to SAS files 


    Use cases tested:

    1.  check extract_sas and extract_trajectory
        a.  extract_trajectory and/or extract_sas is checked
        b.  neither is checked
    2.  check trajectory file
        a.  extract trajectory is True
            i.  output filename is greater than 4 characters 
            ii. output filename ends with pdb or dcd
        b.  extract trajectory is False
    3.  check input PDB file
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    4.  check input trajectory file
        a.  file is .pdb (value = 0)
            i.  PDB file is valid
            ii. PDB file isn't valid
        b.  file is .dcd (value = 1)
            i.  DCD file is valid
            ii. DCD file isn't valid
    5.  check if input PDB and trajectory files are compatible
        a.  trajectory file is a DCD file
            i.  input PDB and DCD are compatible
            ii. input PDB and DCD are not compatible
        b.  trajectory file is a PDB file
            i.  input PDB and PDB are compatible
            ii. input PDB and PDB are not compatible
    6.  check if SAS type is in range 0-3
    7.  check if SAS path permissions
        a.  no permission error
        b.  permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
    8.  check for files of the right SAS type in SAS path
    9.  check if number of SAS files matches the number of frames in the trajectory file               
   10.  check extract options
        a.  single frame
            i.  value must be an integer
            ii. value must not be greater than the number of frames in the trajectory file (or SAS curves in the SAS directory)
            iii.value must be positive
        b.  range
            i.  both number in range are integers
            ii. range is from low to higher integer
            iii.lower limit of range is greater than 0
            iv. lower and higher limits of range are different
            v.  lower limit of range must be equal to or smaller than the maximum number of frames in the trajectory file
            vi. higher limit of range must be equal to or smaller than the maximum number of frames in the trajectory file
            vii.range must be two integers separated by a hyphen 
        c.  text file
            i.  text file doesn't exist
            ii. text file must have only positive integers
            iii.text file must have numbers less than or equal to the number of frames in the trajectory file
            iv. text file can't contain duplicate numbers
            v.  unknown error encountered 
        d.  weight file
            i.  weight file doesn't exist
            ii. weight file column 3 can only have 0.0 or 1.0
            iii.weights can't all be 0
            iv. weight file column 1 can only have positive integers
            v.  weight file must have structure numbers less than or equal to the number of frames in the trajectory file
            vi. weight file can't contain duplicate structure numbers
            vii.unknown error encountered
        e.  sampling frequency
            i.  value must be positive
            ii. value must be smaller than the number of frames in the trajectory file
            iii.unknown error encountered
   11.  check input file path permissions 
        a.  no permission error
        b.  permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed
  12.  check runname
        a.  check for invalid characters in runname                     
    '''

    def setUp(self):

        gui_mimic_extract_utilities.test_variables(self, paths)

    def extract_important_path(self, return_error):

        string_error = string.split(return_error[0])
        path_list = string.split(string_error[-1], '..')
        important_path = string.split(path_list[-1], "/")[1:]
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
        return error[:-1]

    def test_1(self):
        '''
        test if extract_sas and/or extract_trajectory is checked
        '''
        self.extract_trajectory = False
        self.extract_sas = False
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            "at least one of the options ('extract trajectory' and 'extract SAS') needs to be checked"]
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test output trajectory filename
        '''
        self.output_filename = os.path.join(module_data_path, 'file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'output filename must be greater than four characters long and end with .pdb or .dcd : ' + self.output_filename]
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test output trajectory filename
        '''
        self.output_filename = '.text'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'output filename must be greater than four characters long and end with .pdb or .dcd : ' + self.output_filename]
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test output trajectory filename
        '''
        self.output_filename = 'ofil'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'output filename must be greater than four characters long and end with .pdb or .dcd : ' + self.output_filename]
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if input PDB file exists
        '''
        self.pdb_filename = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)
        print return_error

        ''' check for file error '''
        expected_error = ['file : ' +
                          self.pdb_filename + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if input PDB file is a valid file
        '''
        self.pdb_filename = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          self.pdb_filename + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if input trajectory file exists
        '''
        self.trajectory_filename = os.path.join(
            module_data_path, 'does_not_exist.dcd')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' +
                          self.trajectory_filename + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if input trajectory file is a valid pdb file
        '''
        self.trajectory_filename = os.path.join(
            module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input trajectory file, ' +
                          self.trajectory_filename + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if "dcd" trajectory file has the same number of atoms as input PDB
        '''
        self.trajectory_filename = os.path.join(
            module_data_path, 'non_matching.dcd')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file ' + self.pdb_filename + ' and dcd file ' +
                          self.trajectory_filename + ', are not compatiable']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if "pdb" trajectory file has the same number of atoms as input PDB
        '''
        self.trajectory_filename = os.path.join(
            module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file ' + self.pdb_filename + ' and dcd file ' +
                          self.trajectory_filename + ', are not compatiable']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if SAS type is in range 0-3
        '''
        self.sas_type = '4'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for SAS type error '''
        expected_error = ["sas_type %d not supported!" % int(self.sas_type)]
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if SAS path exists
        '''
        self.sas_paths = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.sas_paths + '  [code = FalseFalseFalse]',
                          'sas path "' + self.sas_paths + '" does not exist']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if SAS directory has read permission
        '''
        ''' make a directory '''
        os.system('mkdir empty_SAS_folder')
        ''' see if you can read the directory '''
        print os.access('empty__SAS_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_SAS_folder')
        ''' see if you can read the directory '''
        print os.access('empty_SAS_folder', os.R_OK)

        self.sas_paths= os.path.join('./','empty_SAS_folder')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.sas_paths +
                          '  [code = TrueFalseTrue]', 'read permission not allowed for sas path "' + self.sas_paths + '"']
        assert_equals(return_error, expected_error)  

        ''' make the directory readable'''
        os.system('chmod a+r empty_SAS_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_SAS_folder')        


    def test_14a(self):
        '''
        test if there are files of proper SAS type (sascalc) in SAS directories (no files in either neutron_D2Op_100 nor neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_100, as it won't move on to neutron_D2Op_0
        '''
        self.sas_paths = os.path.join(module_data_path, 'no_sas_files', 'sascalc', 'neutron_D2Op_100') + \
            ',' + os.path.join(module_data_path, 'no_sas_files',
                               'sascalc', 'neutron_D2Op_0')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path(return_error)
        new_error = [return_error[0][0:84] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: 'sascalc' in folder: ../../data/interface/extract_utilities/no_sas_files/sascalc/neutron_D2Op_100"]
        assert_equals(new_error, expected_error)

    def test_14b(self):
        '''
        test if there are files of proper SAS type (sascalc) in SAS directories (proper files in neutron_D2Op_100 but none in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_0 after successful check of neutron_D2Op_100
        '''
        self.sas_paths = os.path.join(module_data_path, 'no_sas_files1', 'sascalc', 'neutron_D2Op_100') + \
            ',' + os.path.join(module_data_path, 'no_sas_files1',
                               'sascalc', 'neutron_D2Op_0')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path(return_error)
        new_error = [return_error[0][0:84] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: 'sascalc' in folder: ../../data/interface/extract_utilities/no_sas_files1/sascalc/neutron_D2Op_0"]
        assert_equals(new_error, expected_error)

    def test_15(self):
        '''
        test if the number of SAS files matches the number of frames in the trajectory file
        '''
        self.sas_paths = os.path.join(
            module_data_path, 'wrong_number_of_sas_files', 'sascalc', 'neutron_D2Op_100')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            "number of SAS files does not match number of frames in the trajectory files"]
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if single frame option value is an integer
        '''
        self.option = 'single_frame'
        self.local_value = '2.5'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'you entered "2.5" : the number of frames and/or SAS files to be extracted needs to be an integer']
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test if single frame option value is larger than number of frames in dcd file
        '''
        self.option = 'single_frame'
        self.local_value = '300'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'there are 20 frames and/or SAS files in your data path : you requested frame and/or SAS file number 300']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if single frame option value is a positive integer
        '''
        self.option = 'single_frame'
        self.local_value = '-1'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'you entered: "-1": the frame and/or SAS file number to be extracted needs to be a positive integer']
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test if range option value is two integers separated by a hyphen
        '''
        self.option = 'range'
        self.local_value = '50'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'range needs to be two integers separated by a hyphen! : you entered "' + self.local_value + '"']
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        test if range option values are both integers
        '''
        self.option = 'range'
        self.local_value = '10-50.5'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'range needs to be two integers separated by a hyphen : you entered "' + self.local_value + '"']
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test if range option values are from low to higher integer
        '''
        self.option = 'range'
        self.local_value = '100-90'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'range needs to be from low to higher integer : you entered "' + self.local_value + '"']
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test if range option lower limit value is greater than 0
        '''
        self.option = 'range'
        self.local_value = '-1 - 100'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'lower limit of the range needs to be greater than 0 : you entered "-1"']
        assert_equals(return_error, expected_error)

    def test_23(self):
        '''
        test if range option values are different
        '''
        self.option = 'range'
        self.local_value = '50 - 50'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'lower and higher limits in the range should be different : you entered "' + self.local_value + '"']
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test if range option lower limit value is equal to or smaller than the maximum number of frames in the trajectory file
        '''
        self.option = 'range'
        self.local_value = '165-170'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'lower limit of the range needs to be equal or smaller than the maximum number of frames and/or SAS files : you entered "165"']
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test if range option upper limit value is equal to or smaller than the maximum number of frames in the trajectory file
        '''
        self.option = 'range'
        self.local_value = '15-25'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'higher limit of the range needs to be equal or smaller than the maximum number of frames and/or SAS files : you entered "25"']
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test if text file exists
        '''
        self.option = 'text_file'
        self.local_value = 'non_existent_file.txt'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' + self.local_value + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test for unknown error encountered when reading text file
        '''
        self.option = 'text_file'
        self.local_value = os.path.join(module_data_path, 'weird_file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'encountered an unknown error reading text_file: ' + self.local_value]
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test for non-positive integer in text file
        '''
        self.option = 'text_file'
        self.local_value = os.path.join(
            module_data_path, 'negative_number.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'text file can only have positive integers : "-1" was found in the text file']
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test if there are numbers greater than number of frames in trajectory file in text file
        '''
        self.option = 'text_file'
        self.local_value = os.path.join(
            module_data_path, 'number_too_large.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'there are 20 frames and/or SAS files in your data path : you requested frame and/or SAS file number 25 in the text file']
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test for redundant number in text file
        '''
        self.option = 'text_file'
        self.local_value = os.path.join(
            module_data_path, 'redundant_number.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'redundant frame and/or SAS file number "12" found in the text file']
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test if weight file exists
        '''
        self.option = 'weight_file'
        self.local_value = 'non_existent_file.txt'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' + self.local_value + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test for unknown error encountered when reading weight file
        '''
        self.option = 'weight_file'
        self.local_value = os.path.join(
            module_data_path, 'weird_weight_file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'encountered an unknown error reading weight_file: ' + self.local_value]
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test for all zeros in column 3 of weight file
        '''
        self.option = 'weight_file'
        self.local_value = os.path.join(
            module_data_path, 'all_zeros_weight_file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'all weights in weight file are zero which means you will not extract any structure or SAS profiles']
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test for numbers other than 0 or 1 in column 3 of weight file
        '''
        self.option = 'weight_file'
        self.local_value = os.path.join(
            module_data_path, 'bad_column_three_weight_file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'weight file column 3 can only have 0 or 1 : 5.000000 was found']
        assert_equals(return_error, expected_error)

    def test_35(self):
        '''
        test for non-positive integer in weight file
        '''
        self.option = 'weight_file'
        self.local_value = os.path.join(
            module_data_path, 'negative_number_weight_file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'weight file column 1 can only have positive integers : "-7" was found in the weight file']
        assert_equals(return_error, expected_error)

    def test_36(self):
        '''
        test if there are structure numbers greater than number of frames in trajectory file in weight file
        '''
        self.option = 'weight_file'
        self.local_value = os.path.join(
            module_data_path, 'number_too_large_weight_file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'there are 20 frames and/or SAS files in your data path : frame and/or SAS file number 25 was found in the weight file']
        assert_equals(return_error, expected_error)

    def test_37(self):
        '''
        test for redundant structure number in weight file
        '''
        self.option = 'weight_file'
        self.local_value = os.path.join(
            module_data_path, 'redundant_number_weight_file.txt')
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'redundant frame and/or SAS file number "19" found in the weight file']
        assert_equals(return_error, expected_error)

    def test_38(self):
        '''
        test for non-positive integer in sampling frequency
        '''
        self.option = 'sampling_frequency'
        self.local_value = '-1'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'the sampling frequency needs to be a positive number : you entered "' + self.local_value + '"']
        assert_equals(return_error, expected_error)

    def test_39(self):
        '''
        test if the sampling frequency is greater than number of frames in trajectory file 
        '''
        self.option = 'sampling_frequency'
        self.local_value = '165'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'the sampling frequency needs to be smaller than the number of frames and/or SAS files in your data path : you entered "165"']
        assert_equals(return_error, expected_error)

    def test_40(self):
        '''
        test for unknow error in sampling frequency
        '''
        self.option = 'sampling_frequency'
        self.local_value = 'not a number'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'encountered an unknown error reading the sampling frequency : ' + self.local_value]
        assert_equals(return_error, expected_error)

    def test_41(self):
        '''
        test if path exists
        '''
        self.path = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_extract_utilities.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path + '  [code = FalseFalseFalse]',
                          'path does not exist']
        assert_equals(return_error, expected_error)

    def test_42(self):
        '''
        test if directory has read permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)

        self.path = os.path.join('./', 'empty_folder')
        return_error = gui_mimic_extract_utilities.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.path + '  [code = TrueFalseTrue]', 'read permission not allowed']
        assert_equals(return_error, expected_error)

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')        

    def test_43(self):
        '''
        test if directory has write permission
        '''
        self.path = os.path.join(module_data_path, 'no_write_permission')
        return_error = gui_mimic_extract_utilities.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path +
                          '  [code = TrueTrueFalse]', 'write permission not allowed']
        assert_equals(return_error, expected_error)

    def test_44(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
