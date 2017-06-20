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
import sassie.simulate.energy_minimization.gui_mimic_energy_minimization as gui_mimic_energy_minimization
#import gui_mimic_energy_minimization as gui_mimic_energy_minimization

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'energy_minimization') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Energy_Minimization_Filter(MockerTestCase):

    '''
    System integration test for minimize_filter.py / sassie 1.0

    Test to see whether minimize_filter catches improper input.

    Inputs tested:

	runname:                        string      run_name 
    infile:                         string      input pdb or dcd file name
    pdbfile:                        string      input (reference) pdb file name
    outfile:                        string      output dcd file name
    nsteps:                         integer     number of minimization steps
    parmfile:                       string      path and name of topology file
    psffile: 		                string      psf file name
    ncpu	:                           integer     number of cpus to use
    keepout:                        integer     keep output files (0==no, 1==yes)
    dcdfreq:                        integer     save individual dcd frequency
    md                              integer     md flag (0=min, 1=min+md, 2=min+md+min)
    mdsteps                         integer     number of md steps
    dielect                         float       solvent dielectric constant
    temperature                     float       temperature

    Inputs not tested:

    path                            string      input path (not in variables)
    infiletype                      string      input file type (pdb or dcd) determined by file name suffix
    use_external_input_file         boolean     flag to use external input file (True or False)
    external_input_file             string      external input file name
    velocity_restart_file           string      velocity restart file name
    extended_system_restart_file    string      extended system restart file name
    charmm_parameter_file           string      name of user-supplied CHARMM parameter file

    Use cases tested:

    1.  check if runname has incorrect character
    2.  check if ncpu >= 1
    3.  check if dcdfreq > 0
    4.  check if nsteps >= dcdfreq
    5.  check if keepout = 0 or 1
    6.  check if nsteps >= 20
    7.  check if parmfile exists
    8.  check if psffile exists
    9.  check if input pdb or dcd file exists
        a)  file doesn't exist
        b)  file exists
            i.  file is not a valid pdb or dcd file
            ii. file is a valid pdb or dcd file
                A. file is not compatible with psf file
                B. file is compatible with psf file
    10.  check if reference pdb file exists
        a)  file doesn't exist
        b)  file exists
            i.  file is not a valid pdb or dcd file
            ii. file is a valid pdb or dcd file
                A. file is not compatible with input pdb or dcd file
                B. file is compatible with input pdb or dcd file
    11. check md parameters
        a)  md != 0, 1 or 2
        b)  md = 0
        c)  md = 1 or 2
            i.  mdsteps > 0                 
            ii. mdsteps is a multiple of 20
            iv. dielect >= 0                
            v.  temperature > 0                              
                        
    '''

    def setUp(self):

        gui_mimic_energy_minimization.test_variables(self, paths)


    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_energy_minimization.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if ncpu >= 1
        '''
        self.ncpu = '0'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ncpu needs to be >= 1, ncpu = 0']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if DCD write frequency > 0
        '''
        self.dcdfreq = '0'
        return_error = gui_mimic_energy_minimization.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['dcd write frequency must be greater than zero, dcd write frequency = 0']
        assert_equals(return_error, expected_error)


    def test_4(self):
        '''
        test if nsteps > dcdfreq 
        '''
        self.nsteps = '25'
        self.dcdfreq = '50'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['max steps < dcd write frequency, nsteps = 25']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if keepout = 0 or 1
        '''
        self.keepout = '2'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['keepout == 0 for "no" and 1 for "yes", keepout = 2']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if nsteps >= 20
        '''
        self.nsteps = '10'
        self.dcdfreq = '10'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['max steps < 20 : max steps = 10']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if parmfile exists
        '''
        self.resparmfile = './does_not_exist.pdb'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ./does_not_exist.pdb does not exist',
        'check parameter path and filename: ./does_not_exist.pdb JEC',
        'check parameter path and filename: ./does_not_exist.pdb JEC']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if psffile exists
        '''
        self.psffile = 'does_not_exist.psf'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : does_not_exist.psf does not exist', '\ncan not find psf file in "./"\n']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if input pdb or dcd file (infile) exists
        '''
        self.infile = 'does_not_exist.dcd'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : does_not_exist.dcd does not exist', 'check input filename (dcd)  path + filename : does_not_exist.dcd']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if reference pdb file exists
        '''
        self.pdbfile = 'does_not_exist.pdb'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : does_not_exist.pdb does not exist', 'check input filename (pdb)  path + filename : does_not_exist.pdb']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if input pdb file is a valid file
        '''
        self.infile = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file '+self.infile+' is not a valid pdb or dcd file']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if input dcd file is a valid file
        '''
        self.infile = os.path.join(module_data_path, 'not_valid.dcd')
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file '+self.infile+' is not a valid pdb or dcd file']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if input pdb file is compatible with psf file
        '''
        self.psffile = os.path.join(module_data_path, 'non_matching.psf')
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file '+self.infile+' and psf file '+self.psffile+' are not compatible']
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test if input dcd file is compatible with psf file
        '''
        self.infile = os.path.join(dcd_data_path, 'ten_mer.dcd')
        self.infiletype = 'dcd'
        self.psffile = os.path.join(module_data_path, 'non_matching.psf')
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file '+self.infile+' and psf file '+self.psffile+' are not compatible']
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test if reference pdb file is a valid pdb file
        '''
        self.pdbfile = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file' +
                          self.pdbfile + ' is not valid']
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if input dcd file is compatible with reference pdb file
        '''
        self.infile = os.path.join(dcd_data_path, 'ten_mer.dcd')
        self.infiletype = 'dcd'
        self.pdbfile = os.path.join(module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file '+self.pdbfile+' and input file '+self.infile+' are not compatible']
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test if input pdb file is compatible with reference pdb file
        '''
        self.pdbfile = os.path.join(module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['reference pdb file '+self.pdbfile+' and input pdb file '+self.infile+' are not compatible']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if md = 0, 1 or 2
        '''
        self.md = '3'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['md == 0 for "min only", 1 for "min + md" and 2 for "min + md + min", md = : 3']
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test if mdsteps > 0
        '''
        self.md = '1'
        self.mdsteps = '0'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of md steps needs to be greater than zero : 0']
        assert_equals(return_error, expected_error)


    def test_20(self):
        '''
        test if mdsteps is a multiple of 20
        '''
        self.md = '1'
        self.mdsteps = '18'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of md steps needs to be a multiple of 20 : 18']
        assert_equals(return_error, expected_error)



    def test_21(self):
        '''
        test if dielectric >= 0
        '''
        self.md = '1'
        self.dielect = '-1.0'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['dielectric needs to be >= zero : -1.0']
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test if temperature > 0
        '''
        self.md = '1'
        self.temperature = '0.0'
        return_error = gui_mimic_energy_minimization.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['temperature needs to be > zero : 0.0']
        assert_equals(return_error, expected_error)

                    

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
