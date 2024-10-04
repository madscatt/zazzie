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
import sassie.calculate.gui_mimic_em_to_sans as gui_mimic_em_to_sans
#import gui_mimic_em_to_sans as gui_mimic_em_to_sans

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'em_to_sans') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Em_To_Sans_Filter(MockerTestCase):

    '''
    System integration test for em_to_sans.py / sassie 1.0

    EM_TO_SANS is the module that reads a three-dimensional file
    (either Gaussian Cube or a binary MRC electron density file) and 
    then represent each voxel with an occupation > threshold
    as a scattering center represented as beads (C-atoms), 
    that is then written as PDB/xyz files
    and the SANS profile is calculated using the "scat" program

    Inputs tested:

    runname:        string      project name                                          
    emfiletype:     integer     EM map file type :  2 values (0: cube; 1: mrc)
    inputpath:      string      input file path
    emdensityfile   string      path and name of input EM density file
    threshold:      float       EM map threshold value
    pdbfile:        string      output PDB file name
    sansfile:       string      output SANS intensity file name
    npoints:        integer     number of SANS data points
    qmax:           float       maximum q value for SANS data
    plotflag:       integer     flag for plotting interpolated data   : 2 values (0: no plot; 1: plot)


    Use cases tested:

    1.  check if runname has incorrect character
    2.  check input file path permissions 
        a.  no permission error
        b. permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed     
    3.  check if emfiletype is 0 or 1
    4.  check if npoints > 1
    5.  check if qmax > 0
    6.  check if threshold > 0
    7.  check if plotflag is 0 or 1
    7.  check if emdensity file exists

     
    '''

    def setUp(self):

        gui_mimic_em_to_sans.test_variables(self, paths)

    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_em_to_sans.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)


    def test_2(self):
        '''
        test if path exists
        '''
        self.inputpath = os.path.join(module_data_path,'non_existent_path')
        return_error = gui_mimic_em_to_sans.run_module(self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.inputpath + '  [code = FalseFalseFalse]',
   'path does not exist']
        assert_equals(return_error, expected_error)    
        
    def test_3(self):
        '''
        test if directory has read permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder')
        ''' see if you can read the directory '''
        print os.access('empty_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_folder')
        ''' see if you can read the directory '''
        print os.access('empty_folder', os.R_OK)

        self.inputpath= os.path.join('./','empty_folder')
        return_error = gui_mimic_em_to_sans.run_module(self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.inputpath + '  [code = TrueFalseTrue]', 'read permission not allowed']
        assert_equals(return_error, expected_error)  

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')

    def test_4(self):
        '''
        test if directory has write permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder1')
        ''' see if you can write to the directory '''
#        print os.access('empty_folder1', os.W_OK)
        ''' make the directory un-writeable'''
        os.system('chmod a-w empty_folder1')
        ''' see if you can write to the directory '''
#        print os.access('empty_folder', os.W_OK)

        self.inputpath = os.path.join('./', 'empty_folder1')
        return_error = gui_mimic_em_to_sans.run_module(
            self, file_check=True)
#        print 'return_error: ', return_error

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.inputpath + '  [code = TrueTrueFalse]', 'write permission not allowed']
#        print 'expected_error: ', expected_error
        assert_equals(return_error, expected_error)

        ''' make the directory writeable'''
        os.system('chmod a+w empty_folder1')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder1')        


    def test_5(self):
        '''
        test if emfiletype is 0 or 1
        '''
        self.emfiletype = '2'
        return_error = gui_mimic_em_to_sans.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['emfiletype == 0 for "cube file" and 1 for "mrc file", emfiletype = 2']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if plotflag is 0 or 1
        '''
        self.plotflag = '2'
        return_error = gui_mimic_em_to_sans.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['plotflag == 0 for no plotting and 1 for plotting, plotflag = 2']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if npoints > 1
        '''
        self.npoints = '1'
        return_error = gui_mimic_em_to_sans.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['npoints needs to be > 1, npoints = 1']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test for qmax > 0
        '''
        self.qmax = '0'
        return_error = gui_mimic_em_to_sans.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['qmax needs to be > 0, qmax = 0.0']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test for threshold > 0
        '''
        self.threshold = '0'
        return_error = gui_mimic_em_to_sans.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['threshold need to be > 0, threshold = 0.0']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if emdensityfile exists
        '''
        self.emdensityfile = os.path.join(other_data_path,'does_not_exist.mrc')
        return_error = gui_mimic_em_to_sans.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file : '+self.emdensityfile+' does not exist', 'check input emdensity file path + filename']
        assert_equals(return_error, expected_error)        

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
