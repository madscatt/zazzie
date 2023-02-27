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
import sys
import shutil
import numpy
import multiprocessing

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
    __file__)), '..', '..', 'data', 'analyze', 'chi_square_filter') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Merge_Utilities(MockerTestCase):

    '''
    System integration test for chi_square_filter.py / sassie 1.0

    CHI_SQUARE_FILTER is the main module that contains the functions
    that are used to compare interpolated experimental data
    to the synthetic data sets generated from structures files
    using the modules in CALCULATE.

    Some filtering options based on Rg and X2 are provided.


    INPUT:	

            Name of output directories
            Paths to synthetic SAS files
            Names of interpolated data sets
            I(0) values corresponding to interpolated data sets
            Scattering calculator used to obtain synthetic SAS files
            Type of X2 file comparison
            Number of weight files to be generated
            Names of weight files
            Basis strings for weight files

    OUTPUT:

            Each output directory contains:

            Best and worst sas profiles from single structure
            Average profile for the entire ensemble
            List of X2 and Rg with each filename
            Rg and X2 without filename (for plotting)
            Text file of sas profile plot from Gnuplot
            Directory of synthetic SAS files scaled to I(0)


    Use cases:

    1.  Input paths
        a.  more than one input path is designated (SasCalc: multiple contrasts, neutron and x-ray)
        b.  one input path is designated

    2.  Synthetic data option
        a.  sastype 0 (SasCalc) is chosen
        b.  sastype 1 (xtal2sas) is chosen
        c.  sastype 2 (Cryson) is chosen
        d.  sastype 3 (Crysol) is chosen

    3.  X2 option
        a.  option 0 (X2) is chosen
        b.  option 1 (reduced X2) is chosen
        c.  option 2 (Pearson's X2) is chosen
        d.  option 3 (R-factor) is chosen

    4.  Weight files
        a.  no weight files are written
        b.  one or more weight files are written 

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


    Test tree:

    project name

******************************
    single input file path
******************************
    saspath             saspath             saspath             saspath        
    sasinfile           sasintfile          sasintfile          sasintfile            
    io                  io                  io                  io           
    no weight files     no weight files     no weight files     no weight files
    SasCalc             Xtal2sas            Cryson              Crysol                         
                        X2, reduced X2, Pearson's X2, R-factor


    saspath             saspath             saspath             saspath        
    sasinfile           sasintfile          sasintfile          sasintfile            
    io                  io                  io                  io           
    two weight files    two weight files    two weight files    two weight files
    SasCalc             Xtal2sas            Cryson              Crysol                         
                        X2, reduced X2, Pearson's X2, R-factor


**********************************
    multiple input file paths
**********************************

    NOTE:  Since all of the synthetic data and reduced X2 options are tested individually
    above, the following two options were chosen to test multiple input file paths:

    input path 1, input path 2          input path 1, input path 2
    sasintfile 1, sasintfile 2          sasintfile 1, sasinfile 2
    io 1, io 2                          io 1, io 2
    SasCalc                             SasCalc
    reduced X2                          reduced X2
    no weight files                     weight file 1, weight file 2**

    **The same two basis_strings and weight_file_names are used for both input paths.      

    '''

    module = 'chi_square_filter'

    def setUp(self):

        gui_mimic_chi_square_filter.test_variables(self, paths)

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

    def check_dir_trees_equal(self, dir1, dir2):
        '''
        compares directories recursively as well as files within them
        ignoring '.DS_Store','averagefile.txt', 'bestworstfile.txt', 'x2file.txt
        '''
        dirs_cmp = filecmp.dircmp(
            dir1, dir2, ['.DS_Store', 'averagefile.txt', 'bestworstfile.txt', 'x2file.txt'])
        if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or \
                len(dirs_cmp.funny_files) > 0:
            return False
        (_, mismatch, errors) = filecmp.cmpfiles(
            dir1, dir2, dirs_cmp.common_files, shallow=False)
        if len(mismatch) > 0 or len(errors) > 0:
            return False
        for common_dir in dirs_cmp.common_dirs:
            new_dir1 = os.path.join(dir1, common_dir)
            new_dir2 = os.path.join(dir2, common_dir)
            if not self.check_dir_trees_equal(new_dir1, new_dir2):
                return False
        return True

    def check_dir_trees_equal1(self, dir1, dir2):
        '''
        compares directories recursively as well as files within them
        ignoring '.DS_Store','averagefile.txt', 'bestworstfile.txt', 'x2file.txt
        as well as *.sassie_json and *.sassie_log files
        applies to the xtal2sas, cryson and crysol sas data types since
        *.sassie_json and *.sassie_log files are in the chi_square_filter directory 
        along with the outputs from chi_square_filter
        these files must be ignored since they have date stamps in their names
        '''

        dirs_cmp = filecmp.dircmp(dir1, dir2,['.DS_Store','averagefile.txt', 'bestworstfile.txt', 'x2file.txt'])
#        if len(dirs_cmp.left_only)>0 or len(dirs_cmp.right_only)>0 or \
#            len(dirs_cmp.funny_files)>0:
        if len(dirs_cmp.right_only)>0 or len(dirs_cmp.funny_files)>0:
            return False
        (_, mismatch, errors) =  filecmp.cmpfiles(
            dir1, dir2, dirs_cmp.common_files)
        if len(mismatch)>0 or len(errors)>0:
            return False
        for common_dir in dirs_cmp.common_dirs:
            new_dir1 = os.path.join(dir1, common_dir)
            new_dir2 = os.path.join(dir2, common_dir)
            if not self.check_dir_trees_equal1(new_dir1, new_dir2):
                return False
        return True


    def test_1(self):
        '''
        single input path:  test sascalc, no weight files, reduced X2
        '''

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:
            #            print 'base_path: ', base_path

            outdirectory = os.path.join(self.runname, self.module, base_path)
#            print 'outdirectory: ', outdirectory
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_red_x2_nowt', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_2(self):
        '''
        single input path:  test sascalc, no weight files, X2
        '''

        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))
        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_x2_nowt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_3(self):
        '''
        single input path:  test sascalc, no weight files, Pearson's X2
        '''

        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))
        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_p_x2_nowt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_4(self):
        '''
        single input path:  test sascalc, no weight files, R-factor
        '''

        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))
        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_r_nowt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_5(self):
        '''
        single input path:  test sascalc, two weight files, reduced X2
        '''

        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_25.txt' + ',' + 'x2_19_27.txt'
        self.basis_string = 'x2<25' + ',' + '(rg<19.4) and (x2<27)'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))
        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_red_x2_wt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_6(self):
        '''
        single input path:  test sascalc, two weight files, X2
        '''

        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_900.txt' + ',' + 'x2_19_850.txt'
        self.basis_string = 'x2<900' + ',' + '(rg<19.4) and (x2<850)'
        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))
        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_x2_wt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_7(self):
        '''
        single input path:  test sascalc, two weight files, Pearson's X2
        '''

        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_0046.txt' + ',' + 'x2_19_0045.txt'
        self.basis_string = 'x2<0.0046' + ',' + '(rg<19.4) and (x2<0.0045)'
        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))
        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_p_x2_wt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_8(self):
        '''
        single input path:  test sascalc, two weight files, R-factor
        '''

        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_039.txt' + ',' + 'x2_19_0375.txt'
        self.basis_string = 'x2<0.039' + ',' + '(rg<19.4) and (x2<0.0375)'
        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))
        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_r_wt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_9(self):
        '''
        single input path:  test xtal2sas, no weight files, reduced X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
#        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_red_x2_nowt')
#        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_10(self):
        '''
        single input path:  test xtal2sas, no weight files, X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'
        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_x2_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_11(self):
        '''
        single input path:  test xtal2sas, no weight files, Pearson's X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'
        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_p_x2_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_12(self):
        '''
        single input path:  test xtal2sas, no weight files, R-factor
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'
        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_r_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_13(self):
        '''
        single input path:  test xtal2sas, two weight files, reduced X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_10.txt' + ',' + 'x2_19_12.txt'
        self.basis_string = 'x2<10' + ',' + '(rg<19.4) and (x2<12.0)'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_red_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_14(self):
        '''
        single input path:  test xtal2sas, two weight files, X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_400.txt' + ',' + 'x2_19_300.txt'
        self.basis_string = 'x2<400' + ',' + '(rg<19.4) and (x2<300)'
        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_15(self):
        '''
        single input path:  test xtal2sas, two weight files, Pearson's X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_003.txt' + ',' + 'x2_19_003.txt'
        self.basis_string = 'x2<0.003' + ',' + '(rg<19.4) and (x2<0.003)'
        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_p_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_16(self):
        '''
        single input path:  test xtal2sas, two weight files, R-factor
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'xtal2sas')
        self.sastype = '1'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_03.txt' + ',' + 'x2_19_03.txt'
        self.basis_string = 'x2<0.03' + ',' + '(rg<19.4) and (x2<0.03)'
        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_r_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_17(self):
        '''
        single input path:  test cryson, no weight files, reduced X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
#        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_red_x2_nowt')
#        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_18(self):
        '''
        single input path:  test cryson, no weight files, X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'
        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_x2_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_19(self):
        '''
        single input path:  test cryson, no weight files, Pearson's X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'
        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_p_x2_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_20(self):
        '''
        single input path:  test cryson, no weight files, R-factor
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'
        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_r_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_21(self):
        '''
        single input path:  test cryson, two weight files, reduced X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_82.txt' + ',' + 'x2_19_84.txt'
        self.basis_string = 'x2<82' + ',' + '(rg<19.64) and (x2<84)'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_red_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_22(self):
        '''
        single input path:  test cryson, two weight files, X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_2570.txt' + ',' + 'x2_19_2580.txt'
        self.basis_string = 'x2<2570' + ',' + '(rg<19.64) and (x2<2580)'
        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_23(self):
        '''
        single input path:  test cryson, two weight files, Pearson's X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_0123.txt' + ',' + 'x2_19_0125.txt'
        self.basis_string = 'x2<0.0123' + ',' + '(rg<19.64) and (x2<0.0125)'
        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_p_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_24(self):
        '''
        single input path:  test cryson, two weight files, R-factor
        '''

        self.saspaths = os.path.join(
            other_data_path, 'diUb', 'run_00', 'cryson')
        self.sastype = '2'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_062.txt' + ',' + 'x2_19_0625.txt'
        self.basis_string = 'x2<0.062' + ',' + '(rg<19.64) and (x2<0.0625)'
        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_r_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_25(self):
        '''
        single input path:  test crysol, no weight files, reduced X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
#        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_red_x2_nowt')
#        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_26(self):
        '''
        single input path:  test crysol, no weight files, X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'
        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_x2_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_27(self):
        '''
        single input path:  test crysol, no weight files, Pearson's X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'
        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_p_x2_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_28(self):
        '''
        single input path:  test crysol, no weight files, R-factor
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'
        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_r_nowt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_29(self):
        '''
        single input path:  test crysol, two weight files, reduced X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_25.txt' + ',' + 'x2_26_28.txt'
        self.basis_string = 'x2<25' + ',' + '(rg<26.6) and (x2<28)'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_red_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_30(self):
        '''
        single input path:  test crysol, two weight files, X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_1950.txt' + ',' + 'x2_26_2000.txt'
        self.basis_string = 'x2<1950' + ',' + '(rg<26.6) and (x2<2000)'
        self.reduced_x2 = '0'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_31(self):
        '''
        single input path:  test crysol, two weight files, Pearson's X2
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_008.txt' + ',' + 'x2_26_01.txt'
        self.basis_string = 'x2<0.008' + ',' + '(rg<26.6) and (x2<0.01)'
        self.reduced_x2 = '2'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_p_x2_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_32(self):
        '''
        single input path:  test crysol, two weight files, R-factor
        '''

        self.saspaths = os.path.join(
            other_data_path, 'RNA', 'run_00', 'crysol')
        self.sasintfiles = os.path.join(other_data_path, 'trunc.dat')
        self.io = '0.031'
        self.sastype = '3'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_012.txt' + ',' + 'x2_26_013.txt'
        self.basis_string = 'x2<0.12' + ',' + '(rg<26.6) and (x2<0.13)'
        self.reduced_x2 = '3'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_r_wt')
        assert_equals(self.check_dir_trees_equal1(
            outdirectory, correct_outdirectory), True)

    def test_33(self):
        '''
        two input paths:  test sascalc, no weight files, reduced X2
        '''

        self.saspaths = os.path.join(other_data_path, 'pai-vn', 'run_2', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(other_data_path, 'pai-vn',
                               'run_2', 'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:
            #            print 'base_path: ', base_path

            outdirectory = os.path.join(self.runname, self.module, base_path)
#            print 'outdirectory: ', outdirectory
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_mult_red_x2_nowt', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def test_34(self):
        '''
        two input paths:  test sascalc, two weight files, reduced X2
        '''

        self.saspaths = os.path.join(other_data_path, 'pai-vn', 'run_2', 'sascalc', 'neutron_D2Op_85') + \
            ',' + os.path.join(other_data_path, 'pai-vn',
                               'run_2', 'sascalc', 'neutron_D2Op_0')
        self.sasintfiles = os.path.join(
            other_data_path, '85p1.dat') + ',' + os.path.join(other_data_path, '0p.dat')
        self.io = '0.013,0.85'
        self.number_of_weight_files = '2'
        self.weight_file_names = 'x2_lt_50.txt' + ',' + '0_x2_8_20_rg_40.txt'
        self.basis_string = 'x2<50' + ',' + '(rg<40) and (rg>20) and (x2<8)'

        gui_mimic_chi_square_filter.run_module(self)

        ''' confirm correct file are in output directory '''

        sas_paths = [x.strip() for x in self.saspaths.split(',')]
        base_paths = []

        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:

            outdirectory = os.path.join(self.runname, self.module, base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_mult_red_x2_wt', base_path)
            assert_equals(self.check_dir_trees_equal(
                outdirectory, correct_outdirectory), True)

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__ == '__main__':
    main()
