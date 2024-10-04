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
import tarfile

import sassie.sasmol.sasmol as sasmol
import sassie.tools.gui_mimic_merge_utilities as gui_mimic_merge_utilities
#import gui_mimic_merge_utilities as gui_mimic_merge_utilities

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'tools', 'merge_utilities') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Merge_Utilities(MockerTestCase):

    '''
    System integration test for merge_utilities.py / sassie 1.0

	MERGE_UTILITIES is a module that allows one to merge coordinates and
	scattering profiles from a series of runs.  Note that each coordinate
	scattering profile set is for the same molecule.

	INPUT: 
		Name of output directory		
		List of run_name directories and dcd filenames	
		List of scattering path names and sas profile type
		Extraction option

	OUTPUT:
	
		New single DCD file concatenated from multiple runs
		New single SAS profile directory with re-numbered filenames	

		
    Use cases:

    1.  Merge trajectory is chosen
        a.  input files are PDB files
        b.  input files are DCD files
        c.  output file is a PDB file
        d.  output file is a DCD file
        
    2.  Merge SAS is chosen
        a.  SAS option 0 (SasCalc) is chosen
        b.  SAS option 1 (xtal2sas) is chosen
        c.  SAS option 2 (cryson) is chosen
        d.  SAS option 3 (crysol) is chosen

    3.  Merge trajectory and merge SAS are chosen
        options a-d for cases 1 and 2 apply

    Selection options (apply to all cases above):
        a.  all
        b.  weight file
        c.  periodic
    

    Inputs tested:

    runname:                string      project name
    pdb_file                string      input pdb file
    trajectory_names        string      input dcd files to be merged (number of files = number of runs)
    output_filename:        string      output dcd file 
    number_of_runs:         integer     number of dcd files and/or SAS runs to be merged                                           
    local_value:            string      value of merge_type_option (not used, list of weight file names, periodic value to skip) 
    merge_option:           integer     merge option (0==merge dcd/pdb/sas profiles, 1==merge sas only, 2==merge dcd/pdb only )                 
    merge_type_option:      integer     merge type option (0==all, 1==weight files, 2==periodic)               
    sas_type:               integer     integer depends on SAS file type (0==SasCalc, 1==Xtal2sas, 2==Cryson, 3==Crysol)    
    sas_paths:              string      paths to SAS files 




    Test tree:

    project name


*******************************************
    merge trajectory and SAS, >100k files
********************************************
 
    reference PDB           reference PDB
    input DCD               input DCD
    output PDB              output PDB
    SasCalc                 SasCalc
    weight files            periodic

    NOTE:  all option not tested
    
    '''

    module = 'merge_utilities'

    @classmethod
    def setUpClass(cls):
        print 'extracting tarfile'
        filename = os.path.join(other_data_path,'extract_large_file_number_0.tar.gz')
        tar = tarfile.open(filename, "r:gz")
        tar.extractall()
        print 'finished extracting'
        tar.close

 
    @classmethod
    def tearDownClass(cls):
        print 'deleting extracted directory'
        def robust_rmtree(path, max_retries=6):
            dt = 1
            for i in range(max_retries):
                try:
                    shutil.rmtree(path)
                    return
                except OSError:
                    time.sleep(dt)
                    dt *= 2
            shutil.rmtree(path, ignore_error=True)

        tarpath = os.path.join('./','extract_large_file_number_0')
        if os.path.exists(tarpath):
            robust_rmtree(tarpath)

    def setUp(self):

       gui_mimic_merge_utilities.test_variables(self, paths)


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

    def check_dir_trees_equal(self,dir1, dir2):
        '''
        compares directories recursively as well as files within them
        '''
        dirs_cmp = filecmp.dircmp(dir1, dir2)
        if len(dirs_cmp.left_only)>0 or len(dirs_cmp.right_only)>0 or \
            len(dirs_cmp.funny_files)>0:
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
        merge trajectory and SAS:  test dcd input, pdb output, weight files, Sascalc, >100k files
        '''

        self.merge_option = '0' 
        self.merge_type_option = '1'
        self.pdb_file = os.path.join(pdb_data_path, 'final_small_peptide.pdb')
        self.trajectory_names = os.path.join(dcd_data_path,'vmd_test_small_peptide.dcd')+','+os.path.join(dcd_data_path,'vmd_test_small_peptide.dcd')
        self.output_filename = 'weight_files_small_peptide_sascalc.pdb'
        self.sas_type = '0'
        self.local_value = os.path.join(other_data_path,'small_peptide_weight_file_extract.txt')+','+os.path.join(other_data_path,'small_peptide_weight_file_merge.txt')         
        self.sas_paths = os.path.join('./', 'extract_large_file_number_0','run_1','sascalc')+','+os.path.join('./', 'extract_large_file_number_0','run_1', 'sascalc')  
  
        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'weight_files_small_peptide_sascalc.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)
        ''' confirm correct SAS files are in merged directories '''

        self.sas_paths = [x.strip() for x in self.sas_paths.split(',')]
        for sas_path in self.sas_paths:
            base_paths = []
#            print 'sas_path: ',sas_path

            for root, dirs, files in os.walk(sas_path, topdown=False):
                for name in dirs:
                    base_paths.append(os.path.basename(os.path.normpath(os.path.join(root, name))))

#            print 'base_paths: ', base_paths

        for base_path in base_paths:
            print 'base_path: ', base_path       

            outdirectory = os.path.join(self.runname, self.module, 'sascalc', base_path)
#            print 'outdirectory: ',outdirectory
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_small_peptide_weight_files', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                                                                                            
    def test_2(self):
        '''
        merge trajectory and SAS:  test dcd input, pdb output, periodic, Sascalc, >100k files
        '''

        self.merge_option = '0' 
        self.merge_type_option = '2'
        self.pdb_file = os.path.join(pdb_data_path, 'final_small_peptide.pdb')
        self.trajectory_names = os.path.join(dcd_data_path,'vmd_test_small_peptide.dcd')+','+os.path.join(dcd_data_path,'vmd_test_small_peptide.dcd')
        self.output_filename = 'periodic_small_peptide_sascalc.pdb'
        self.sas_type = '0'
        self.local_value = '10000'         
        self.sas_paths = os.path.join('./', 'extract_large_file_number_0','run_1','sascalc')+','+os.path.join('./', 'extract_large_file_number_0','run_1', 'sascalc')  
  
        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'periodic_small_peptide_sascalc.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)
        ''' confirm correct SAS files are in merged directories '''

        self.sas_paths = [x.strip() for x in self.sas_paths.split(',')]
        for sas_path in self.sas_paths:
            base_paths = []
#            print 'sas_path: ',sas_path

            for root, dirs, files in os.walk(sas_path, topdown=False):
                for name in dirs:
                    base_paths.append(os.path.basename(os.path.normpath(os.path.join(root, name))))

#            print 'base_paths: ', base_paths

        for base_path in base_paths:
            print 'base_path: ', base_path       

            outdirectory = os.path.join(self.runname, self.module, 'sascalc', base_path)
#            print 'outdirectory: ',outdirectory
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_small_peptide_periodic', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

        
    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

