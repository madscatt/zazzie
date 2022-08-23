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
import sassie.tools.data_interpolation.gui_mimic_data_interpolation as gui_mimic_data_interpolation
#import gui_mimic_data_interpolation as gui_mimic_data_interpolation

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'data_interpolation') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Data_Interpolation_Filter(MockerTestCase):

    '''
    System integration test for data_interpolation.py / sassie 1.0

    data_interpolation is the module that calculates an approximate data set
	with a defined number of points and grid spacing.  The output of
	this file is used to compare to synthetic profiles generated in
	subsequent modules.

    Inputs tested:

    run_name:      string      project name                                          
    expdata:      string      experimental data file                : 2 files (SANS, SAXS)
    ofile:        string      output interpolated data file         : 1 for each above
    io:           float       I(0) value of the data, determined from Guinier, P(r), or other suitable method
    ioe:          float       error on I(0) value
    dq:           float       delta q value (grid spacing between q values)
    maxpoints:    integer     maximum number of q values
    plotflag:     integer     flag for plotting interpolated data   : 2 values (0: no plot; 1: plot)


    Use cases tested:

    1.  check if expdata exists
    2.  check if expdata is a directory
    3.  check if run_name has incorrect character
    4.  check for negative dq
    5.  check for maxpoints less than 2
    6.  check dq and maxpoints              NOT TESTED: non-float error occurs in input_filter before getting to this point
        a.  dq can not be read
        b.  maxpoints can not be read
    7.  check if dq*maxpoints is greater than maximum q
    8.  check if q values are increasing


     
    '''

    def setUp(self):

        gui_mimic_data_interpolation.test_variables(self, paths)


    def test_1(self):
        '''
        test if expdata file exists
        '''
        self.expdata = os.path.join(
            module_data_path, 'does_not_exist.sub')
        return_error = gui_mimic_data_interpolation.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input data file, ' +
                          self.expdata + ', does not exist']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if expdata file is a directory
        '''
        self.expdata = os.path.join(
            module_data_path, 'is_a_directory')
        return_error = gui_mimic_data_interpolation.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input data file, ' +
                          self.expdata + ', is a directory!']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if run_name has incorrect character
        '''
        self.run_name = 'run_&'
        return_error = gui_mimic_data_interpolation.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test for dq less than 0
        '''
        self.dq = '-0.05'
        return_error = gui_mimic_data_interpolation.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['dq needs to be > 0 : -0.05']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if maxpoints is less than 2
        '''
        self.maxpoints = '1'
        return_error = gui_mimic_data_interpolation.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['maxpoints needs to be greater than 2 : 1']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if dq*maxpoints is greater than maximum q
        '''
        self.dq = '0.05'
        return_error = gui_mimic_data_interpolation.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['input parameters compared to data in file, '+self.expdata+', are not valid : dq*maxpoints > q-range in data file [q-range = last q-value minus first q-value] : dq = '+self.dq+', maxpoints = '+self.maxpoints]
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test for non-increasing q values
        '''
        self.expdata = os.path.join(
            module_data_path, 'sans_data_q_not_increasing.sub')
        return_error = gui_mimic_data_interpolation.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['in your input file, '+self.expdata[3:]+', found that q-values that are not increasing: duplicate points or successive q-values where q[0] < q[1]']
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.run_name):
            shutil.rmtree(self.run_name)

if __name__ == '__main__':
    main()
