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
import sassie.analyze.density_plot.gui_mimic_density_plot as gui_mimic_density_plot
#import gui_mimic_density_plot as gui_mimic_density_plot

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
    __file__)), '..', '..', 'data', 'interface', 'density_plot') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Density_Plot_Filter(MockerTestCase):

    '''
    System integration test for density_plot_filter.py / sassie 1.0

    Test to see whether density_plot_filter catches improper input.

    Inputs tested:

            runname:        project name
            path:           path name for input files
            dcdfile:        input trajectory filename (pdb or dcd)
            pdbfile:        reference pdb name
            xlength:        x boxlength
            ylength:        y boxlength
            zlength:        z boxlength
            gridsp:         grid spacing (angstroms)            
            nsegments:	    number of segments
            segvariables:   number of flexible regions, high and low residues, basis string, segment name
            save_occupancy: save the unweighted raw cube data ('Y' or 'N')
            equalweights:   use equalweights (1=yes) or weights from file (0=no)
            weightsfile:    filename containing weights per structure


    Use cases tested:

 
    '''

    def setUp(self):

        gui_mimic_density_plot.test_variables(self, paths)

    def extract_important_path(self, return_error):

        string_error = string.split(return_error[0])
        path_list = string.split(string_error[-1], '..')
        important_path = string.split(path_list[-1], "/")[1:]
        error = os.path.join('..', '..')
        for this_path in important_path:
            error += os.sep + this_path
        return error[:-1]

    def test_21(self):
        '''
        test for unknown error encountered when reading weight file
        '''

        self.equalweights = "0"
        self.weightsfile = os.path.join(
            module_data_path, 'weird_weight_file.txt')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ["unable to open and read your weight file : " + os.path.basename(self.weightsfile)]
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)
            
if __name__ == '__main__':
    main()
