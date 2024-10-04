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

import sasmol.system as system
import sassie.tools.data_interpolation.gui_mimic_data_interpolation as gui_mimic_data_interpolation

import filecmp
from unittest import main, TestCase
from unittest.mock import Mock, patch
from nose.tools import assert_equals

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'tools', 'data_interpolation') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Data_Interpolation(TestCase):

    '''
    System integration test for data_interpolation.py / sassie 1.0

    data_interpolation is the module that calculates an approximate data set
        with a defined number of points and grid spacing.  The output of
        this file is used to compare to synthetic profiles generated in
        subsequent modules.

    Use cases:

    1.  input file is a SANS data file
    2.  input file is a SAXS data file

    Inputs tested:

    run_name:      string      project name
    expdata:      string      experimental data file                : 2 files (SANS, SAXS)
    ofile:        string      output interpolated data file         : 1 for each above
    io:           float       I(0) value of the data, determined from Guinier, P(r), or other suitable method
    ioe:          float       error on I(0) value
    dq:           float       delta q value (grid spacing between q values)
    maxpoints:    integer     maximum number of q values
    plotflag:     integer     flag for plotting interpolated data   : 2 values (0: no plot; 1: plot)



    Test tree:


                                                *********************
                                                *   project name    *
                                                * input/output path *
                                                *********************
                                                *                   *
                                              *                       *
                                            *                           *
                            ********************                   *********************
                            *    SANS data     *                   *     SAXS Data     *
                            ********************                   *********************
                                    *                                       *
                                    *                                       *
                           *********************                   *********************
                           *      I(0)         *                   *      I(0)         *
                           *   I(0) error      *                   *   I(0) error      *
                           *    delta q        *                   *    delta q        *
                           *number of q values *                   *number of q values *
                           *********************                   *********************

    '''

    module = 'data_interpolation'

    def setUp(self):

        gui_mimic_data_interpolation.test_variables(self, paths)

    @patch('sassie.tools.data_interpolation.gui_mimic_data_interpolation.run_module')
    def test_1(self, MockClass):
        '''
        test SANS data file
        '''
        instance = MockClass.return_value
        instance.some_method.return_value = 'expected_value'

        result = gui_mimic_data_interpolation.run_module(self)

        ''' confirm output data file is correct '''
        outfile = os.path.join(self.run_name, self.module, self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'sans_data.dat')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm that output stn_data file is correct '''
        outfile = os.path.join(self.run_name, self.module,
                               'stn_'+self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'stn_sans_data.dat')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)


    '''
    def test_2(self):
        """
        test SAXS data file
        """
        self.expdata = os.path.join(other_data_path, 'trunc2a_saxs.sub')
        self.ofile = 'trunc2a.dat'
        self.io = '0.031'
        self.dq = '0.007'
        self.maxpoints = '72'

        gui_mimic_data_interpolation.run_module(self)

        # confirm output data file is correct
        outfile = os.path.join(self.run_name, self.module, self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'trunc2a.dat')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        # confirm that output stn_data file is correct
        outfile = os.path.join(self.run_name, self.module,
                               'stn_'+self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'stn_trunc2a.dat')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)
    '''

    def tearDown(self):
        if os.path.exists(self.run_name):
            shutil.rmtree(self.run_name)


if __name__ == '__main__':
    main()
