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
import sassie.analyze.chi_square_filter.gui_mimic_chi_square_filter as gui_mimic_chi_square_filter
#import gui_mimic_chi_square_filter as gui_mimic_chi_square_filter

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
    __file__)), '..', '..', 'data', 'interface', 'chi_square_filter') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Chi_Square_Filter_Filter(MockerTestCase):

    '''
    System integration test for chi_square_filter_filter.py / sassie 1.0

    Test to see whether chi_square_filter_filter catches improper input.

    Inputs tested:

    runname:                string  project name
    saspaths:               string  paths to synthetic SAS files (number of saspaths = number of sasintfiles = number of I(0) values)
    sasintfiles:            string  names of interpolated data sets
    io:                     float   I(0) values
    sastype:                integer scattering calculator used to generate synthetic SAS files (0=sascalc, 1=xtal2sas, 2=cryson, 3=crysol)
    reduced_x2:		        integer type of X2 comparison (0=chi-squared, 1=reduced chi-squared, 2=Pearson's chi-squared, 3=R-factor)
    number_of_weight_files: integer number of weight files to be generated
    weight_file_names:      string  names of weight files (number_of_weight_files = number of names = number of basis_strings)
    basis_string:           string  basis strings for weight files, i.e., 'x2<1', '(x2<5) and (Rg<30)'
    plot_flat               integer flag for plotting results

    Use cases tested:

    1.  check runname
        a.  check for invalid characters in runname              
    2.  check interpolated data file(s)
        a.  number of interpolated data files matches number of SAS data paths
        b.  interpolated data file exists
        c.  interpolated data file can be opened and read
        d.  interpolated data file isn't empty or doesn't have all-zero intensity entries
        e.  interpolated data file has 3 columns
            i.   file has less than 3 columns
            ii.  file has greater than 3 columns
            iii. file has 3 columns
        f.  interpolated data file has positive Q values
    3.  check plotflag
        a.  plotflag is 0 or 1
        b.  plotflag is not 0 or 1
    4.  check basis string
        a.  length of basis string matches number of weight files (check only if number of weight files > 0)
        b.  basis string doesn't contain upper case characters, i.e., Rg, RG, X2
        c.  basis string doesn't contain incorrect characters, i.e., 'rg<30/'
        d.  basis string doesn't contain incorrect expression, i.e., 'rg-30' 
        d.  basis string doesn't contain mismatched parenthesis, i.e., '(rg<20) and (x2<5))'
        e.  basis string doesn't contain wrong parenthesis, i.e., '(rg<20( and (x2<5)'
        f.  basis string doesn't contain incorrect syntax, i.e., 'rg lt 10'  
    5.  check weight file names
        a.  number of weight files matches number of weight file names      
    6.  check sastype
        a.  sastype is 0, 1, 2 or 3
        b.  sastype is not 0, 1, 2 or 3
    7.  check SAS data paths            
        a.  no permission error
        b.  permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
    8.  check if there are files of proper SAS type in SAS data paths    NOTE:  two SAS data paths tested
    9.  check if number of q values in SAS files matches number of Q values in interpolated data file
    10. check if q value or spacing in SAS files matches q value or spacing in interpolated data file
    11. check if number of q values is correct in SAS files     NOTE:  two SAS data paths tested
    12. check if q values are correct in SAS files              NOTE:  two SAS data paths tested
    13. check if number of log files is equal to number of SAS files    NOTE:  two SAS data paths tested
    14. check if I(0) is greater than 0                         NOTE:  two SAS data paths tested
        a.  I(0) > 0
        b.  I(0) </= 0


    '''

    def setUp(self):

        gui_mimic_chi_square_filter.test_variables(self, paths)

    def extract_important_path(self, return_error):

        string_error = string.split(return_error[0])
        path_list = string.split(string_error[-1], '..')
        important_path = string.split(path_list[-1], "/")[1:]
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
        return error[:-1]

    def extract_important_path1(self, return_error):

        string_error = string.split(return_error[0])
#        print 'string_error: ', string_error
        path_list = string.split(string_error[-4], '..')
#        print 'path_list: ', path_list
        important_path = string.split(path_list[-1], "/")[1:]
#        print 'important_path: ', important_path
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
        return error

    def extract_important_path2(self, return_error):

        string_error = string.split(return_error[0])
#        print 'string_error: ', string_error
        path_list = string.split(string_error[-15], '..')
#        print 'path_list: ', path_list
        important_path = string.split(path_list[-1], "/")[1:]
#        print 'important_path: ', important_path
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
        return error

    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if number of interpolated data files matches number of SAS data paths
        '''
        self.saspaths = os.path.join(other_data_path, 'pai-vn', 'run_2', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(other_data_path, 'pai-vn',
                               'run_2', 'sascalc', 'neutron_D2Op_0')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'number of interpolated data files does not match number of SAS data paths']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if interpolated data file exists
        '''
        self.sasintfiles = os.path.join(other_data_path, 'non-existent.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
#        print 'return error: ', return_error
        relative_path = self.extract_important_path1(return_error)
 #       print 'relative_path: ', relative_path
        new_error = [return_error[0][:7] + relative_path + ' does not exist']
        new_error.append(
            'check interpolated data file path + filename : non-existent.dat')
#        print 'new_error: ', new_error

        ''' check for value error '''
        expected_error = ['file : ../../data/other_common/non-existent.dat does not exist',
                          'check interpolated data file path + filename : non-existent.dat']
#        print 'expected error: ', expected_error
        assert_equals(new_error, expected_error)

    def test_4(self):
        '''
        test if interpolated data file can be opened and read
        '''
        ''' make a directory '''
        os.system('mkdir empty_SAS_folder')
        ''' see if you can read the directory '''
        print os.access('empty__SAS_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_SAS_folder')
        ''' see if you can read the directory '''
        print os.access('empty_SAS_folder', os.R_OK)

        self.sasintfiles = os.path.join('./', 'empty_SAS_folder')
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for path error '''
        expected_error = [
            'unable to open and read your interpolated data file : empty_SAS_folder']
        assert_equals(return_error, expected_error)

        ''' make the directory readable'''
        os.system('chmod a+r empty_SAS_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_SAS_folder')

    def test_5(self):
        '''
        test if interpolated data file is empty 
        '''
        self.sasintfiles = os.path.join(module_data_path, 'empty_file.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'no lines or all-zero intensity entries found in your SAS data file : empty_file.dat']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if interpolated data file has all zero-intensity values
        '''
        self.sasintfiles = os.path.join(module_data_path, 'zero_intensity.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'no lines or all-zero intensity entries found in your SAS data file : zero_intensity.dat']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if interpolated data file has more than 3 columns
        '''
        self.sasintfiles = os.path.join(module_data_path, 'extra_columns.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'the number of columns should be 3 in your interpolated data file : extra_columns.dat  \n0.000          0.8500              0.1040E-01     0.11\n']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if interpolated data file has less than 3 columns
        '''
        self.sasintfiles = os.path.join(module_data_path, 'two_columns.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'the number of columns should be 3 in your interpolated data file : two_columns.dat  \n0.000          0.8500 \n']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if interpolated data file has 3 columns with no # character in the 3rd column
        '''
        self.sasintfiles = os.path.join(
            module_data_path, 'comment_in_column_3.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'the number of columns should be 3 in your interpolated data file : comment_in_column_3.dat  \n0.000          0.8500   #testing\n']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if interpolated data file has negative Q values
        '''
        self.sasintfiles = os.path.join(
            module_data_path, 'negative_q_value.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'negative Q value found in your SAS data file : negative_q_value.dat  \n-0.030          0.7027              0.1712E-02\n']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if plotflag is 0 or 1
        '''
        self.plotflag = '2'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['plot flag needs to be 1 or 0 ... you entered: 2']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if length of basis strings matches number of weight files
        '''
        self.number_of_weight_files = '1'
        self.basis_string = 'x2<50' + ',' + '(rg<40) and (rg>20) and (x2<8)'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'number of weight files does not match length of basis string']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test for upper case characters in basis string
        '''
        self.number_of_weight_files = '2'
        self.basis_string = 'x2<50' + ',' + '(Rg<40) and (rg>20) and (x2<8)'
        self.weight_file_names = 'x2_lt_50.txt' + ',' + '0_x2_8_20_rg_40.txt'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['expression "(Rg<40) and (rg>20) and (x2<8)" not understood!',
                          'you may need to change "Rg|RG|X2" to lower cases']
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test for grammar error '-' instead of '<', '>', '='
        '''
        self.number_of_weight_files = '1'
        self.basis_string = 'x2-10'
        self.weight_file_names = 'x2_lt_50.txt'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['grammar wrong in expression "x2-10"']
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test for grammar error '/' added to expression
        '''
        self.number_of_weight_files = '2'
        self.basis_string = 'x2<50/' + ',' + '(rg<40) and (rg>20) and (x2<8)'
        self.weight_file_names = 'x2_lt_50.txt' + ',' + '0_x2_8_20_rg_40.txt'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['grammar wrong in expression "x2<50/"']
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test for grammar error '()' in expression
        '''
        self.number_of_weight_files = '1'
        self.basis_string = '()x2<50'
        self.weight_file_names = 'x2_lt_50.txt'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['grammar wrong in expression "()x2<50"']
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test for grammar error wrong parenthesis
        '''
        self.number_of_weight_files = '1'
        self.basis_string = '(rg<40( and (rg>20) and (x2<8)'
        self.weight_file_names = '0_x2_8_20_rg_40.txt'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'grammar is wrong in expression "(rg<40( and (rg>20) and (x2<8)"']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test for wrong syntax '(rg lt 40)'
        '''
        self.number_of_weight_files = '2'
        self.basis_string = 'x2<50' + ',' + '(rg lt 40) and (rg>20) and (x2<8)'
        self.weight_file_names = 'x2_lt_50.txt' + ',' + '0_x2_8_20_rg_40.txt'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'expression "(rg lt 40) and (rg>20) and (x2<8)" not understood!']
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test for mismatched parenthesis
        '''
        self.number_of_weight_files = '2'
        self.basis_string = 'x2<50' + ',' + '(rg<40) and (rg>20) and (x2<8))'
        self.weight_file_names = 'x2_lt_50.txt' + ',' + '0_x2_8_20_rg_40.txt'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'parenthesis in expression "(rg<40) and (rg>20) and (x2<8))" does not match!']
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        test if number of weight files matches number of weight file names
        '''
        self.number_of_weight_files = '2'
        self.basis_string = 'x2<50' + ',' + '(rg<40) and (rg>20) and (x2<8)'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'number of weight files does not match number of weight file names']
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test for proper sastype value (0, 1, 2 or 3)
        '''
        self.sastype = '4'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ["sastype 4 not supported!"]
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test if SAS path exists
        '''
        self.saspaths = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.saspaths + '  [code = FalseFalseFalse]',
                          'sas path "' + self.saspaths + '" does not exist']
        assert_equals(return_error, expected_error)

    def test_23(self):
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

        self.saspaths = os.path.join('./', 'empty_SAS_folder')
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.saspaths +
                          '  [code = TrueFalseTrue]', 'read permission not allowed for sas path "' + self.saspaths + '"']
        assert_equals(return_error, expected_error)

        ''' make the directory readable'''
        os.system('chmod a+r empty_SAS_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_SAS_folder')

    def test_24(self):
        '''
        test if there are scattering files of proper SAS type (sascalc) in SAS directories (no files in either neutron_D2Op_85 nor neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_85, as it won't move on to neutron_D2Op_0
        '''
        self.saspaths = os.path.join(module_data_path, 'no_sas_files', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'no_sas_files',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path(return_error)
        new_error = [return_error[0][0:84] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: 'sascalc' in folder: ../../data/interface/chi_square_filter/no_sas_files/sascalc/neutron_D2Op_85"]
        assert_equals(new_error, expected_error)

    def test_25(self):
        '''
        test if there are scattering files of proper SAS type (sascalc) in SAS directories (proper files in neutron_D2Op_85 but none in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_0 after successful check of neutron_D2Op_85
        '''
        self.saspaths = os.path.join(module_data_path, 'no_sas_files1', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'no_sas_files1',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path(return_error)
        new_error = [return_error[0][:84] + relative_path]

        ''' check for file error '''
        expected_error = [
            "there are no scattering files found for the selected sas-type: 'sascalc' in folder: ../../data/interface/chi_square_filter/no_sas_files1/sascalc/neutron_D2Op_0"]
        assert_equals(new_error, expected_error)

    def test_26(self):
        '''
        test if number of Q values in theoretical SAS profiles matches that in interpolated data file
        '''
        self.sasintfiles = os.path.join(module_data_path, 'wrong_q_values.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'The number of Q values in supplied interpolated data file does not match Q value in theoretical SAS profiles\n']
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test if Q value or spacing in theoretical SAS profiles matches that in interpolated data file
        '''
        self.sasintfiles = os.path.join(
            module_data_path, 'wrong_q_spacing.dat')

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['Q value or spacing in supplied interpolated data file does not match Q value or spacing in theoretical SAS profiles\n',
                          'Q value in your SAS file : run_0_00001.iq:\n...[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]...\n  does not match that in your experimental file:\n...[0.015, 0.02, 0.03, 0.04, 0.05, 0.065, 0.07, 0.08, 0.09, 0.1]...']
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test if the number of q values are correct in SAS files (file with incorrect values in neutron_D2Op_85 but correct values in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_85, as it won't move on to neutron_D2Op_0
        '''
        self.saspaths = os.path.join(module_data_path, 'wrong_number_of_q_values_in_sas_files', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_q_values_in_sas_files',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path2(return_error)
        new_error = [return_error[0][:24] +
                     relative_path + return_error[0][-89:]]

        ''' check for file error '''
        expected_error = [
            'number of q-values of "/../../data/interface/chi_square_filter/wrong_number_of_q_values_in_sas_files/sascalc/neutron_D2Op_85/run_2_00136.iq" does not match "run_2_00001.iq".\nPlease verify the files in the SAS path are consistent!']
        assert_equals(new_error, expected_error)

    def test_29(self):
        '''
        test if the number of q values are correct in SAS files (files with correct values in neutron_D2Op_85 but incorrect values in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_0 after checking neutron_D2Op_85
        '''
        self.saspaths = os.path.join(module_data_path, 'wrong_number_of_q_values_in_sas_files1', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_q_values_in_sas_files1',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path2(return_error)
        new_error = [return_error[0][:24] +
                     relative_path + return_error[0][-89:]]

        ''' check for file error '''
        expected_error = [
            'number of q-values of "/../../data/interface/chi_square_filter/wrong_number_of_q_values_in_sas_files1/sascalc/neutron_D2Op_0/run_2_00042.iq" does not match "run_2_00001.iq".\nPlease verify the files in the SAS path are consistent!']
        assert_equals(new_error, expected_error)

    def test_30(self):
        '''
        test if the q values are correct in SAS files (file with incorrect values in neutron_D2Op_85 but correct values in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_85, as it won't move on to neutron_D2Op_0
        '''
        self.saspaths = os.path.join(module_data_path, 'wrong_q_values_in_sas_files', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'wrong_q_values_in_sas_files',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path2(return_error)
        new_error = [return_error[0][:13] +
                     relative_path + return_error[0][-87:]]

        ''' check for file error '''
        expected_error = [
            'q-values of "../../data/interface/chi_square_filter/wrong_q_values_in_sas_files/sascalc/neutron_D2Op_85/run_2_00153.iq" do not match "run_2_00001.iq".\nPlease verify the files in the SAS path are consistent!']
        assert_equals(new_error, expected_error)

    def test_31(self):
        '''
        test if the q values are correct in SAS files (files with correct values in neutron_D2Op_85 but incorrect values in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_0 after checking neutron_D2OP_85
        '''
        self.saspaths = os.path.join(module_data_path, 'wrong_q_values_in_sas_files1', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'wrong_q_values_in_sas_files1',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path2(return_error)
        new_error = [return_error[0][:13] +
                     relative_path + return_error[0][-87:]]

        ''' check for file error '''
        expected_error = [
            'q-values of "../../data/interface/chi_square_filter/wrong_q_values_in_sas_files1/sascalc/neutron_D2Op_0/run_2_00072.iq" do not match "run_2_00001.iq".\nPlease verify the files in the SAS path are consistent!']
        assert_equals(new_error, expected_error)

    def test_32(self):
        '''
        test if the number of log files is the same as the number of SAS files (incorrect number of log files 
        in neutron_D2Op_85 but correct number in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_85, as it won't move on to neutron_D2Op_0
        '''
        self.saspaths = os.path.join(module_data_path, 'wrong_number_of_log_files', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_log_files',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path(return_error)
        new_error = [return_error[0][:91] + relative_path]

        ''' check for file error '''
        expected_error = [
            'the number of log files does not match the number of theoretical SAS profiles in sas path: ../../data/interface/chi_square_filter/wrong_number_of_log_files/sascalc/neutron_D2Op_85']
        assert_equals(new_error, expected_error)

    def test_33(self):
        '''
        test if the number of log files is the same as the number of SAS files (correct number of log files 
        in neutron_D2Op_85 but incorrect number in neutron_D2Op_0)
        test should fail upon checking neutron_D2Op_0 after checking neutron_D2Op_85
        '''
        self.saspaths = os.path.join(module_data_path, 'wrong_number_of_log_files1', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(module_data_path, 'wrong_number_of_log_files1',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)
        '''extract the relative path of the files for new error message'''
        relative_path = self.extract_important_path(return_error)
        new_error = [return_error[0][:91] + relative_path]

        ''' check for file error '''
        expected_error = [
            'the number of log files does not match the number of theoretical SAS profiles in sas path: ../../data/interface/chi_square_filter/wrong_number_of_log_files1/sascalc/neutron_D2Op_0']
        assert_equals(new_error, expected_error)

    def test_34(self):
        '''
        test if number of I(0) values matches that of SAS data paths
        '''
        self.io = '0.013,0.85'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = [
            'number I(0) values does not match number of of SAS data paths']
        assert_equals(return_error, expected_error)

    def test_35(self):
        '''
        test if I(0) greater than 0 (first path not OK, second path OK)
        '''
        self.saspaths = os.path.join(other_data_path, 'pai-vn', 'run_2', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(other_data_path, 'pai-vn', 'run_2',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0, 0.85'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['I(0): "0.0" should be greater than 0']
        assert_equals(return_error, expected_error)

    def test_36(self):
        '''
        test if I(0) greater than 0 (first path OK, second path not OK)
        '''
        self.saspaths = os.path.join(other_data_path, 'pai-vn', 'run_2', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(other_data_path, 'pai-vn', 'run_2',
                               'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,-0.85'

        return_error = gui_mimic_chi_square_filter.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['I(0): "-0.85" should be greater than 0']
        assert_equals(return_error, expected_error)

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
