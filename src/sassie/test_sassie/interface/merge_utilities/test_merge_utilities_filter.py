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
import sassie.tools.merge_utilities.gui_mimic_merge_utilities as gui_mimic_merge_utilities
#import gui_mimic_merge_utilities as gui_mimic_merge_utilities

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
    __file__)), '..', '..', 'data', 'interface', 'merge_utilities') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Merge_Utilities_Filter(MockerTestCase):

    '''
    System integration test for merge_utilities_filter.py / sassie 1.0

    Test to see whether merge_utilities_filter catches improper input.

    Inputs tested:

    runname:                string      project name
    pdb_file                string      input pdb file
    trajectory_names        string      input dcd files to be merged (number of files = number of runs)
    output_filename:        string      output dcd file 
    number_of_runs:         integer     number of dcd files and/or SAS runs to be merged                                           
    local_value:            string      value of merge_type_option (not used, list of weight file names, periodic value to skip) 
    merge_option:           integer     merge option (0==merge dcd/pdb/sas profiles, 1==merge sas only, 2==merge dcd/pdb only )                 
    merge_type_option:      integer     merge type option (0==all, 1==weight files, 2==periodic)               
    sas_type:               integer     integer depends on SAS file type (0==SasCalc, 1==Xtal2sas, 2==Cryson, 3==Crysol)    
    sas_paths:              string      paths to SAS files 


    Use cases tested:

    1.  check options/types
        a.  merge_option is not 0, 1 or 2
        b.  merge_type_option is not 0, 1 or 2
        c.  sas_type is not 0, 1, 2 or 3
    2.  check for pdb/dcd files
        a.  check if output filename is greater than four characters long and ends with .pdb or .dcd
        b.  check if number of trajectory files matches number of runs
        c.  check input pdb file
            i.  check if input pdb file exists
            ii. check if input pdb file is a valid pdb file
        d.  check trajectory file (dcd)
            i.  check if dcd file exists
            ii. check if dcd file is valid
            iii.check if dcd file and input pdb file are compatible
        e.  check trajectory file (pdb)    
            i.  check if pdb file exists
            ii. check if pdb file is valid
            iii.check if pdb file and input pdb file are compatible
    3.  check input file path permissions                           
        a.  no permission error
        b.  permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
    4.  check if number of SAS folders matches number of runs               
    5.  check if chosen SAS type is compatible with files in SAS paths
    6.  check if SAS files of proper type are found in SAS directories      NOTE: tested for sas_type 0 and 1
    7.  check if number of SAS files matches number of frames in dcd file   NOTE: tested for sas_type 0 and 1
    8.  check if number of weight files matches number of runs
    9.  weight file checks
        a.  weight file doesn't exist
        b. weight file column 3 can only have 0.0 or 1.0
        c. weights can't all be 0
        d. weight file column 1 can only have positive integers
        e. weight file must have structure numbers less than or equal to the number of frames in the trajectory file
        f. weight file can't contain duplicate structure numbers
        g. unknown error encountered
    10.  sampling frequency checks
        a. value must be positive
        b. value must be smaller than the number of frames in the trajectory file
        c. unknown error encountered   
    11. check runname
        a.  check for invalid characters in runname 
                    
    '''

    def setUp(self):

        gui_mimic_merge_utilities.test_variables(self, paths)

    def extract_important_path0(self, return_error):

        string_error = string.split(return_error[0])
        path_list = string.split(string_error[-1], '..')
        important_path = string.split(path_list[-1], "/")[1:]
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
        return error

    def extract_important_path1(self, return_error):

        string_error = string.split(return_error[0])
        path_list = string.split(string_error[-13], '..')
        important_path1 = string.split(path_list[-1], "/")[1:]
        error = os.path.join('..', '..')
        for this_path in important_path1:
            error += os.sep + this_path
        return error        

    def extract_important_path2(self, return_error):

        string_error = string.split(return_error[0])
        path_list = string.split(string_error[-4], '..')
        important_path2 = string.split(path_list[-1], "/")[1:]
        error = os.path.join('..', '..')
        for this_path in important_path2:
            error += os.sep + this_path
        return error


    def test_1(self):
        '''
        test if merge option is 0, 1 or 2
        '''

        self.merge_option = '3'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'merge option needs to be 0, 1, or 2, you entered : 3']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if merge_type option is 0, 1 or 2
        '''

        self.merge_type_option = '3'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'merge_type option needs to be 0, 1, or 2, you entered : 3']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if sas_type is 0, 1, 2 or 3
        '''

        self.sas_type = '4'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'sas type needs to be 0, 1, 2, or 3, you entered : 4']
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test output trajectory filename
        '''
        self.output_filename = os.path.join(module_data_path, 'file.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'output filename must be greater than four characters long and end with .pdb or .dcd : ' + self.output_filename]
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if number of trajectory files matches number of runs
        '''
        self.number_of_runs = '3'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of trajectory files "2" does not match number of runs "3"']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if input PDB file exists
        '''
        self.pdb_file = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' +
                          self.pdb_file+' does not exist']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if input PDB file is a valid file
        '''
        self.pdb_file = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          self.pdb_file + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_8a(self):
        '''
        test if DCD trajectory files exist (first listed file doesn't exist) 
        '''
        self.trajectory_names = os.path.join(dcd_data_path,'does_not_exist.dcd')+','+os.path.join(dcd_data_path,'run_m2.dcd')       
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' +
                          self.trajectory_names.split(',')[0]+' does not exist']               
        assert_equals(return_error, expected_error)

    def test_8b(self):
        '''
        test if DCD trajectory files exist (second listed file doesn't exist) 
        '''
        self.trajectory_names = os.path.join(dcd_data_path,'run_m1.dcd')+','+os.path.join(dcd_data_path,'does_not_exist.dcd')        
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' +
                          self.trajectory_names.split(',')[1]+' does not exist']                 
        assert_equals(return_error, expected_error)
 
    def test_9(self):
        '''
        test if input trajectory file is a valid dcd file
        '''
        self.number_of_runs = '1'
        self.trajectory_names = os.path.join(
            module_data_path, 'not_valid.dcd')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input trajectory file, ' +
                          self.trajectory_names + ', is not a valid dcd file']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if dcd trajectory file has the same number of atoms as input PDB
        '''
        self.number_of_runs = '1'
        self.trajectory_names = os.path.join(
            module_data_path, 'non_matching.dcd')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file ' + self.pdb_file + ' and dcd file ' +
                          self.trajectory_names + ', are not compatible']                 
        assert_equals(return_error, expected_error)
 
    def test_11(self):
        '''
        test if input trajectory file is a valid pdb file
        '''
        self.number_of_runs = '1'
        self.trajectory_names = os.path.join(
            module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input trajectory file, ' +
                          self.trajectory_names + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if pdb trajectory file has the same number of atoms as input PDB
        '''
        self.number_of_runs = '1'
        self.trajectory_names = os.path.join(
            module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file ' + self.pdb_file + ' and pdb file ' +
                          self.trajectory_names + ', are not compatible']                 
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if number of SAS directories matches number of runs
        '''
        self.merge_option = '1'
        self.number_of_runs = '1'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of SAS folders "2" does not match number of runs "1"']                 
        assert_equals(return_error, expected_error)
 
    def test_14(self):
        '''
        check if chosen SAS type is compatible with files in sas paths
        '''
        self.sas_type = '3'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the SAS type "crysol" you entered is not compatiable with the SAS type in the SAS data path you selected']             
        assert_equals(return_error, expected_error)
        
    def test_15(self):
        '''
        test if SAS path exists
        '''
        self.number_of_runs = '1'
        self.merge_option = '1'
        self.sas_paths = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_merge_utilities.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.sas_paths + '  [code = FalseFalseFalse]',
                          'sas path "' + self.sas_paths + '" does not exist']
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if SAS path has read permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)

        self.number_of_runs = '1'
        self.merge_option = '1'
        self.sas_paths = os.path.join('./', 'empty_folder')
        return_error = gui_mimic_merge_utilities.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.sas_paths + '  [code = TrueFalseTrue]', 'read permission not allowed for sas path "'+self.sas_paths+'"']
        assert_equals(return_error, expected_error)

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder') 

    def test_17a(self):
        '''
        test if there are files of proper SAS type (xtal2sas) in SAS directories (no files in run_0 but proper files in run_1)
        test should fail upon checking run_0, as it won't move on to run_1
        '''
        self.sas_paths = os.path.join(module_data_path, 'no_sas_files_0', 'run_0', 'xtal2sas') + \
            ',' + os.path.join(module_data_path, 'no_sas_files_0',
                               'run_1', 'xtal2sas')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path0(return_error)
        print 'relative path: ', relative_path
        new_error = [return_error[0][0:83] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: xtal2sas in folder: ../../data/interface/merge_utilities/no_sas_files_0/run_0/xtal2sas"]
        assert_equals(new_error, expected_error)

    def test_17b(self):
        '''
        test if there are files of proper SAS type (xtal2sas) in SAS directories (proper files in run_0 but none in run_1)
        test should fail upon checking run_1 after successful check of run_0
        '''
        self.sas_paths = os.path.join(module_data_path, 'no_sas_files_1', 'run_0', 'xtal2sas') + \
            ',' + os.path.join(module_data_path, 'no_sas_files_1',
                               'run_1', 'xtal2sas')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path0(return_error)
        new_error = [return_error[0][0:83] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: xtal2sas in folder: ../../data/interface/merge_utilities/no_sas_files_1/run_1/xtal2sas"]
        assert_equals(new_error, expected_error)

    def test_18a(self):
        '''
        test if there are files of proper SAS type (sascalc) in SAS directories (no files in neutron_D2Op_80/run_0 but proper files in run_1)
        test should fail upon checking run_0, as it won't move on to run_1
        '''
        self.sas_type = '0'
        self.sas_paths = os.path.join(module_data_path, 'no_sas_files_0', 'run_0', 'sascalc') + \
            ',' + os.path.join(module_data_path, 'no_sas_files_0',
                               'run_1', 'sascalc')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path0(return_error)
        print 'relative path: ', relative_path
        new_error = [return_error[0][0:82] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: sascalc in folder: ../../data/interface/merge_utilities/no_sas_files_0/run_0/sascalc/neutron_D2Op_80"]
        assert_equals(new_error, expected_error)

    def test_18b(self):
        '''
        test if there are files of proper SAS type (sascalc) in SAS directories (proper files in neutron_D2Op_80/run_0 but none in run_1)
        test should fail upon checking run_1 after successful check of run_0
        '''
        self.sas_type = '0'
        self.sas_paths = os.path.join(module_data_path, 'no_sas_files_1', 'run_0', 'sascalc') + \
            ',' + os.path.join(module_data_path, 'no_sas_files_1',
                               'run_1', 'sascalc')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path0(return_error)
        new_error = [return_error[0][0:82] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: sascalc in folder: ../../data/interface/merge_utilities/no_sas_files_1/run_1/sascalc/neutron_D2Op_80"]
        assert_equals(new_error, expected_error)

    def test_19a(self):
        '''
        test if there are the right number of SAS type (xtal2sas) files in SAS directories (wrong number in run_0 but correct number in run_1)
        test should fail upon checking run_0, as it won't move on to run_1
        '''
        self.sas_paths = os.path.join(module_data_path, 'wrong_number_of_sas_files_0', 'run_0', 'xtal2sas') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_sas_files_0',
                               'run_1', 'xtal2sas')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path1 = self.extract_important_path1(return_error)
        relative_path2 = self.extract_important_path0(return_error)
        new_error = [return_error[0][0:31] + relative_path1 + ' does not match number of frames of the PDB/DCD file in "' + relative_path2]

        ''' check for file error '''
        expected_error = [
            'number of SAS files in folder "'+'../../data/interface/merge_utilities/wrong_number_of_sas_files_0/run_0/xtal2sas" does not match number of frames of the PDB/DCD file in "../../data/dcd_common/run_m1.dcd"']
        assert_equals(new_error, expected_error)

    def test_19b(self):
        '''
        test if there the right number of SAS type (xtal2sas) files in SAS directories (proper files in run_0 but wrong number in run_1)
        test should fail upon checking run_1 after successful check of run_0
        '''
        self.sas_paths = os.path.join(module_data_path, 'wrong_number_of_sas_files_1', 'run_0', 'xtal2sas') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_sas_files_1',
                               'run_1', 'xtal2sas')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path1 = self.extract_important_path1(return_error)
        relative_path2 = self.extract_important_path0(return_error)
        new_error = [return_error[0][0:31] + relative_path1 + ' does not match number of frames of the PDB/DCD file in "' + relative_path2]

        expected_error = [
            'number of SAS files in folder "'+'../../data/interface/merge_utilities/wrong_number_of_sas_files_1/run_1/xtal2sas" does not match number of frames of the PDB/DCD file in "../../data/dcd_common/run_m2.dcd"']
        assert_equals(new_error, expected_error)

    def test_20a(self):
        '''
        test if there are the right number of SAS type (sascalc) files in SAS directories (wrong number in neutron_D2Op_80/run_0 but correct number in run_1)
        test should fail upon checking run_0, as it won't move on to run_1
        '''
        self.sas_type = '0'
        self.sas_paths = os.path.join(module_data_path, 'wrong_number_of_sas_files_0', 'run_0', 'sascalc') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_sas_files_0',
                               'run_1', 'sascalc')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path1 = self.extract_important_path1(return_error)
        relative_path2 = self.extract_important_path0(return_error)
        new_error = [return_error[0][0:31] + relative_path1 + ' does not match number of frames of the PDB/DCD file in "' + relative_path2]

        ''' check for file error '''
        expected_error = [
            'number of SAS files in folder "'+'../../data/interface/merge_utilities/wrong_number_of_sas_files_0/run_0/sascalc/neutron_D2Op_80" does not match number of frames of the PDB/DCD file in "../../data/dcd_common/run_m1.dcd"']
        assert_equals(new_error, expected_error)

    def test_20b(self):
        '''
        test if there the right number of SAS type (sascalc) files in SAS directories (proper files in neutron_D2Op_80/run_0 but wrong number in run_1)
        test should fail upon checking run_1 after successful check of run_0
        '''
        self.sas_type = '0'
        self.sas_paths = os.path.join(module_data_path, 'wrong_number_of_sas_files_1', 'run_0', 'sascalc') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_sas_files_1',
                               'run_1', 'sascalc')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path1 = self.extract_important_path1(return_error)
        relative_path2 = self.extract_important_path0(return_error)
        new_error = [return_error[0][0:31] + relative_path1 + ' does not match number of frames of the PDB/DCD file in "' + relative_path2]

        expected_error = [
            'number of SAS files in folder "'+'../../data/interface/merge_utilities/wrong_number_of_sas_files_1/run_1/sascalc/neutron_D2Op_80" does not match number of frames of the PDB/DCD file in "../../data/dcd_common/run_m2.dcd"']
        assert_equals(new_error, expected_error)


    def test_21(self):
        '''
        test if number of weight files matches number of runs
        '''

        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path, 'weights_file_m1.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of weight files "1" does not match number of runs "2"']                 
        assert_equals(return_error, expected_error)        

    def test_22a(self):
        '''
        test if weight file exists (weight file 1 doesn't exist but weight file 2 exists)
        test should fail upon checking weight file 1, as it won't move on to weight file 2
        '''

        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path, 'non_existent.txt')+','+os.path.join(other_data_path, 'weights_file_m2.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        print 'return_error: ', return_error
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path2(return_error)
        new_error = [return_error[0][0:7] + relative_path + ' does not exist']

        expected_error = [
            'file : ../../data/other_common/non_existent.txt does not exist']
        assert_equals(new_error, expected_error)

    def test_22b(self):
        '''
        test if weight file exists (weight file 1 exisst but weight file 2 doesn't exist)
        test should fail upon checking weight file 2 after successfully checking weight file 1
        '''

        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path, 'weights_file_m1.txt')+','+os.path.join(other_data_path, 'non_existent.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)
        print 'return_error: ', return_error
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path2(return_error)
        new_error = [return_error[0][0:7] + relative_path + ' does not exist']
        print 'new_error: ', new_error

        expected_error = [
            'file : ../../data/other_common/non_existent.txt does not exist']
        assert_equals(new_error, expected_error)

    def test_23(self):
        '''
        test for numbers other than 0 or 1 in column 3 of weight file
        '''

        self.merge_option = '1'
        self.merge_type_option = '1'
        self.number_of_runs = '1'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = os.path.join(
            module_data_path, 'bad_column_three_weight_file.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'weight file column 3 can only have 0 or 1 : 5.000000 was found']
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test for all zeros in column 3 of weight file
        '''

        self.merge_option = '1'
        self.merge_type_option = '1'
        self.number_of_runs = '1'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = os.path.join(
            module_data_path, 'all_zeros_weight_file.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'all weights in weight file are zero which means you will not extract any structure or SAS profiles']
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test for non-positive integer in weight file
        '''

        self.merge_option = '1'
        self.merge_type_option = '1'
        self.number_of_runs = '1'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = os.path.join(
            module_data_path, 'negative_number_weight_file.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'weight file column 1 can only have positive integers : "-7" was found in the weight file']
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test if there are structure numbers greater than number of frames in trajectory file in weight file
        '''

        self.merge_type_option = '1'
        self.number_of_runs = '1'
        self.trajectory_names = os.path.join(dcd_data_path,'run_m1.dcd')
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = os.path.join(
            module_data_path, 'number_too_large_weight_file.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'there are 39 frames and/or SAS files in your data path : frame and/or SAS file number 45 was found in the weight file']
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test for redundant structure number in weight file
        '''

        self.merge_option = '1'
        self.merge_type_option = '1'
        self.number_of_runs = '1'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = os.path.join(
            module_data_path, 'redundant_number_weight_file.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'redundant SAS file number "19" found in the weight file']
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test for unknown error encountered when reading weight file
        '''

        self.merge_option = '1'
        self.merge_type_option = '1'
        self.number_of_runs = '1'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = os.path.join(
            module_data_path, 'weird_weight_file.txt')
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'encountered an unknown error reading weight_file: ' + self.local_value]
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test for non-positive integer in sampling frequency
        '''

        self.merge_option = '1'
        self.merge_type_option = '2'
        self.number_of_runs = '1'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = '-1'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'the sampling frequency needs to be a positive number : you entered "'+self.local_value+'"']
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test if the sampling frequency is larger than the number of frames in trajectory file 
        '''

        self.merge_type_option = '2'
        self.number_of_runs = '1'
        self.trajectory_names = os.path.join(dcd_data_path,'run_m1.dcd')
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = '40'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'the sampling frequency needs to be smaller than the number of frames and/or SAS files in your data path : you entered "'+self.local_value+'"']
        assert_equals(return_error, expected_error)


    def test_31(self):
        '''
        test for unknown error encountered when reading sampling frequency
        '''

        self.merge_option = '1'
        self.merge_type_option = '2'
        self.number_of_runs = '1'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')
        self.local_value = 'not a number'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [
            'encountered an unknown error reading the sampling frequency : '+ self.local_value]
        assert_equals(return_error, expected_error)


    def test_32(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_merge_utilities.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
