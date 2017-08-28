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
import string
import numpy
import multiprocessing

import sasmol.sasmol as sasmol
import sassie.simulate.prody.gui_mimic_prody_anm as gui_mimic_prody
#import gui_mimic_prody_anm as gui_mimic_prody

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'simulate', 'prody') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Data_Interpolation(MockerTestCase):

    '''
    System integration test for prody_anm.py / sassie 2.0

    PRODY is the module that contains functions to setup and run normal mode
	analysis, using ProDy, based on the structure in a supplied
	pdb file.

    Use cases:

    1.  ANM (normal mode) analysis using ProDy 
    2.  user-supplied input command to run a wide variety of "prody" or "evol" commands  (NOT YET IMPLEMENTED)

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


    Test tree:

    project name

***************************
    ANM analysis (default)
***************************

    input PDB
    no advanced input           
                                    
    '''

    module = 'prody'

    def setUp(self):

       gui_mimic_prody.test_variables(self, paths)

    ''' compare the files that should be the same '''

    def test_1(self):
        '''
        test input PDB file, no advanced options
        '''

        #prody expects the input pdb file to be in the ./ directory for further processing, so copy the test pdb file to the test directory 
        cmd = 'cp '+self.pdbfile+' .'
#        print 'cmd: ', cmd
        os.system(cmd)
        #rename self.pdbfile to the file just copied
        filestring = string.split(self.pdbfile,'/')        
        self.pdbfile = filestring[-1]
#        print 'pdbfile: ', self.pdbfile
        gui_mimic_prody.run_module(self)

        ''' confirm _anm_beta.txt file is correct '''
        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_beta.txt')
        correct_outfile = os.path.join(
            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_beta.txt')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

#       file is too big for git -- not checked
#        ''' confirm _anm_covariance.txt file is correct '''
#        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_covariance.txt')
#        correct_outfile = os.path.join(
#            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_covariance.txt')
#        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm _anm_cross-correlations.txt file is correct '''
        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_cross-correlations.txt')
        correct_outfile = os.path.join(
            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_cross-correlations.txt')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm _anm_evalues.txt file is correct '''
        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_evalues.txt')
        correct_outfile = os.path.join(
            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_evalues.txt')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm _anm_evectors.txt file is correct '''
        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_evectors.txt')
        correct_outfile = os.path.join(
            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_evectors.txt')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

#       file is too big for git -- not checked
#        ''' confirm _anm_hessian.txt file is correct '''
#        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_hessian.txt')
#        correct_outfile = os.path.join(
#            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_hessian.txt')
#        assert_equals(filecmp.cmp(outfile, correct_outfile), True)        

        ''' confirm _anm_kirchhoff.txt file is correct '''
        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_kirchhoff.txt')
        correct_outfile = os.path.join(
            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_kirchhoff.txt')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)        

        ''' confirm _anm_sqflucts.txt file is correct '''
        outfile = os.path.join(self.runname, self.module, self.pdbfile[:-4]+'_anm_sqflucts.txt')
        correct_outfile = os.path.join(
            module_data_path,self.runname,self.module,'pdb_none', self.pdbfile[:-4]+'_anm_sqflucts.txt')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)        


        #use this test to see which files are different if any of the above tests fail
        #dir1 = os.path.join(self.runname, self.module)
        #print 'dir1: ', dir1
        #dir2 = os.path.join(
        #    module_data_path,self.runname,self.module,'pdb_none')
        #print 'dir2: ', dir2
        #comparison = filecmp.dircmp(dir2, dir1)
        #comparison.report_full_closure()

        try:
            os.system('rm -f *.pdb ')
        except:
            pass

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

