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
        self.extract_trajectory = 'False'
        self.extract_sas = 'False'
        return_error = gui_mimic_extract_utilities.run_module(
            self, test_filter=True)
        print 'return error: ', return_error

        ''' check for file error '''
        expected_error = [
            "at least one of the options ('extract trajectory' and 'extract SAS') needs to be checked"]
        print 'expected error: ', expected_error
        assert_equals(return_error, expected_error)



    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
