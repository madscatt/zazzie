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

import sasmol.sasmol as sasmol
import sassie.calculate.em_to_sas.gui_mimic_em_to_sas as gui_mimic_em_to_sas
#import gui_mimic_em_to_sas as gui_mimic_em_to_sas

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'calculate', 'em_to_sas') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Em_To_Sas(MockerTestCase):

    '''
    System integration test for em_to_sas.py / sassie 1.0

        EM is the function to read in three-dimensional voxel data and
        calculate the SAS profile through the binary program scat.

        INPUT:  variable descriptions:

                runname:        project name
                emfiletype:     em file type (0=cube; 1=mrc)
                inputpath:      input file path
                emdensityfile:  em filename
                threshold:      treshold cutoff
                npoints:        number of points in sas calculation
                qmax:           q-max in sas calculation

        OUTPUT:

                pdbfile:        output filename (pdb)
                
                files stored in ./"runname"/em_to_sas/ directory:

                sasfile*.sub:  output sas profile
                sasfile*.pr:   output p(r)
                dum.inp:        input file written for scat
                sasfile*.pdb:  pdb file of coordinates used for scat
                sasfile*.xyz:  xyz file of coordinates used for scat


    Use cases:

    1.  Input EM density file
        a.  input file is a mrc file
        b.  input file is a cube file


    Inputs tested:

    Inputs tested:

    runname:        string      project name                                          
    emfiletype:     integer     EM map file type :  2 values (0: cube; 1: mrc)
    inputpath:      string      input file path
    emdensityfile   string      path and name of input EM density file
    threshold:      float       EM map threshold value
    pdbfile:        string      output PDB file name
    sasfile:        string      output SAS intensity file name
    npoints:        integer     number of SAS data points
    qmax:           float       maximum q value for SAS data
    plotflag:       integer     flag for plotting interpolated data   : 2 values (0: no plot; 1: plot)      NO plotting for tests


    Test tree:

    project name
    input/output path

**************************
    input file type 
**************************
    input mrc file          input cube file
 
    '''

    module = 'em_to_sas'

    def setUp(self):

       gui_mimic_em_to_sas.test_variables(self, paths)

    
    def check_dir_trees_equal(self,dir1, dir2):
        '''
        compares directories recursively as well as files within them
        ignoring files with date stamps: *.sassie_json and *.sassie_log
        '''
        dirs_cmp = filecmp.dircmp(dir1, dir2)
#        if len(dirs_cmp.left_only)>0 or len(dirs_cmp.right_only)>0 or \
#            len(dirs_cmp.funny_files)>0:
        if len(dirs_cmp.right_only)>0 or len(dirs_cmp.funny_files)>0:
            return False
        (_, mismatch, errors) =  filecmp.cmpfiles(
            dir1, dir2, dirs_cmp.common_files, shallow=False)
        if len(mismatch)>0 or len(errors)>0:
            return False
        for common_dir in dirs_cmp.common_dirs:
            new_dir1 = os.path.join(dir1, common_dir)
            new_dir2 = os.path.join(dir2, common_dir)
            if not self.check_dir_trees_equal(new_dir1, new_dir2):
                return False
        return True

    def test_1(self):
        '''
        test mrc input
        '''

        gui_mimic_em_to_sas.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'mrc_file')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_2(self):
        '''
        test cube input
        '''

        self.emfiletype = '0'
        self.threshold = '0.01'
        self.emdensityfile = os.path.join(other_data_path,'gag_complete.cube')
        gui_mimic_em_to_sas.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cube_file')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)                      
                       

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

