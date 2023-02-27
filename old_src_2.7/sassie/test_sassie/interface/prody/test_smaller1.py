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
import sassie.simulate.prody.gui_mimic_prody_anm as gui_mimic_prody
#import gui_mimic_prody_anm as gui_mimic_prody

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'prody') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Prody_Filter(MockerTestCase):

    '''
    System integration test for prody_filter.py / sassie 1.0

    Test to see whether prody_filter catches improper input.

    Inputs tested:

    runname:                    string      project name                                           
    pdbfile:                    string      input pdb file             
    number_modes                string      number of normal modes to compute
    number_conformations_samp   string      number of conformation to generate by random sampling of modes
    number_steps_traverse       string      number of steps to traverse each mode in both directions
    rmsd_comformations_samp     string      average RMSD of randomly sampled conformations with respect to initial conformation
    rmsd_traverse               string      maximum RMSD of conformations for trajectory from traversed mode with respect to initial conformation
  
    Inputs not tested (not yet implemented):

    advanced_usage              string      advanced usage option: 0=no; 1=yes 
    advanced_usage_cmd          string      user-supplied ProDy command if advanced_usage = 1    

    Use cases tested:

    1.  check if runname has invalid character
    2.  check pdbfile
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    3.  check if number of normal modes > 0
    4.  check if number of conformations per sample > 0
    5.  check if number of frames to traverse per mode > 0
    6.  check if averaged RMSD > 0
    7.  check if maximum RMSD for traverse mode > 0

    '''

    def setUp(self):

        gui_mimic_prody.test_variables(self, paths)


    def test_2(self):
        '''
        test if pdbfile exists
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_prody.run_module(self, test_filter=True)
        print 'return error: ', return_error

        ''' check for file error ''' 
        expected_error = ['file : ' + self.pdbfile + ' does not exist']
        print 'expected error: ', expected_error        
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
