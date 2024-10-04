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
import sassie.tools.gui_mimic_extract_utilities as gui_mimic_extract_utilities
#import gui_mimic_extract_utilities as gui_mimic_extract_utilities

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'tools', 'extract_utilities') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Extract_Utilities(MockerTestCase):

    '''
    System integration test for extract_utilities.py for LARGE files/directories / sassie 1.0

	EXTRACT_UTILITIES is a module that allows one to extract coordinates and/or
    SAS profiles from PDB/DCD files and/or a directory containing SAS profiles.
	The multi-frame files trajectory is saved into new a PDB/DCD file.  SAS profiles
    are saved to a new directory.   Options are provided to extract
	single structures and/or SAS, structures and/or SAS over a range, by reading a list of frame 
	numbers from a simple text file, via a 'weights' file, or by a user-supplied
    sampling frequency.

	INPUT: 
	
		PDB or DCD file with molecular structure data (multi-frames)
		Extraction option
		Option value
		Name of output file		
        SAS profile type
        Folder containing input SAS profiles

	OUTPUT:
        	
		PDB or DCD file with requested coordinates
        and/or
        Folder containing SAS profiles

    Use cases:

    1.  Extract trajectory is chosen
        a.  input file is a PDB file
        b.  input file is a DCD file
        c.  output file is a PDB file
        d.  output file is a DCD file
        
    2.  Extract SAS is chosen
        a.  SAS option 0 (SasCalc) is chosen
        b.  SAS option 1 (xtal2sas) is chosen
        c.  SAS option 2 (cryson) is chosen
        d.  SAS option 3 (crysol) is chosen

    3.  Extract trajectory and extract SAS are chosen
        options a-d for cases 1 and 2 apply

    Selection options (apply to all cases above):
        a.  single frame
        b.  range
        c.  text file
        d.  weight file
        e.  periodic
        f.  all
    

    Inputs tested:

    runname:                string      project name
    path:                   string      input/output filepath      
    pdb_filename            string      input pdb file
    trajectory_filename     string      input pdb or dcd file                          
    option:                 string      extract option (single_frame, range, all, weight_file, text_file)                 
    local_value:            string      value of option (frame value, range, all, weight_file name, text_file name)                
    output_filename:        string      output pdb or dcd file   
    extract_trajectory:     boolean     extract frames from trajectory (True or False)                
    extract_sas :           boolean     extract corresponding SAS files (True or False)                 
    sas_type:               integer     integer depends on SAS file type (SasCalc, Xtal2sas, Cryson, Crysol)    
    sas_path:               string      path to SAS files 



    Test tree:

    project name
    input/output path
    

*******************************************
    extract trajectory and SAS, >100k files
********************************************
 
    reference PDB       reference PDB   reference PDB   reference PDB   reference PDB
    input DCD           input DCD       input DCD       input DCD       input DCD
    output PDB          output PDB      output PDB      output PDB      output PDB 
    SasCalc             SasCalc         SasCalc         SasCalc         SasCalc
    weight file         text file       range           single frame    periodic

    NOTE:  all option not tested


    '''

    module = 'extract_utilities'

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

       gui_mimic_extract_utilities.test_variables(self, paths)


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
        extract_trajectory and SAS: test dcd input, pdb output, weight file, SasCalc, >100k files
        '''

        self.pdb_filename = os.path.join(pdb_data_path, 'final_small_peptide.pdb')
        self.trajectory_filename = os.path.join(dcd_data_path, 'vmd_test_small_peptide.dcd')
        self.output_filename = 'weight_file_small_peptide_sascalc.pdb'
        self.sas_paths = os.path.join('./', 'extract_large_file_number_0', 'run_1', 'sascalc','neutron_D2Op_100')
        self.option = 'weight_file'
        self.local_value = os.path.join(other_data_path,'small_peptide_weight_file_extract.txt')
  
        gui_mimic_extract_utilities.run_module(self)

        ''' confirm correct SAS files are chosen '''

        self.sas_paths = [x.strip() for x in self.sas_paths.split(',')]
        base_paths = []
        for sas_path in self.sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:
            outdirectory = os.path.join(self.runname, self.module, 'sascalc', base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_small_peptide_weight_file', base_path)
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_2(self):
        '''
        extract_trajectory and SAS: test dcd input, pdb output, text file, SasCalc, >100k files
        '''

        self.pdb_filename = os.path.join(pdb_data_path, 'final_small_peptide.pdb')
        self.trajectory_filename = os.path.join(dcd_data_path, 'vmd_test_small_peptide.dcd')
        self.output_filename = 'text_file_small_peptide_sascalc.pdb'
        self.sas_paths = os.path.join('./', 'extract_large_file_number_0', 'run_1', 'sascalc','neutron_D2Op_100')
        self.option = 'text_file'
        self.local_value = os.path.join(other_data_path,'small_peptide_text_file_extract.txt')
  
        gui_mimic_extract_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'text_file_small_peptide_sascalc.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm correct SAS files are chosen '''

        self.sas_paths = [x.strip() for x in self.sas_paths.split(',')]
        base_paths = []
        for sas_path in self.sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:
            outdirectory = os.path.join(self.runname, self.module, 'sascalc', base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_small_peptide_text_file', base_path)
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_3(self):
        '''
        extract_trajectory and SAS: test dcd input, pdb output, range, SasCalc, >100k files
        '''

        self.pdb_filename = os.path.join(pdb_data_path, 'final_small_peptide.pdb')
        self.trajectory_filename = os.path.join(dcd_data_path, 'vmd_test_small_peptide.dcd')
        self.output_filename = 'chosen_range_small_peptide_sascalc.pdb'
        self.sas_paths = os.path.join('./', 'extract_large_file_number_0', 'run_1', 'sascalc','neutron_D2Op_100')
        self.option = 'range'
        self.local_value = '10000-10005'
  
        gui_mimic_extract_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'chosen_range_small_peptide_sascalc.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm correct SAS files are chosen '''

        self.sas_paths = [x.strip() for x in self.sas_paths.split(',')]
        base_paths = []
        for sas_path in self.sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:
            outdirectory = os.path.join(self.runname, self.module, 'sascalc', base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_small_peptide_range', base_path)
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_4(self):
        '''
        extract_trajectory and SAS: test dcd input, pdb output, single frame, SasCalc, >100k files
        '''

        self.pdb_filename = os.path.join(pdb_data_path, 'final_small_peptide.pdb')
        self.trajectory_filename = os.path.join(dcd_data_path, 'vmd_test_small_peptide.dcd')
        self.output_filename = 'single_frame_small_peptide_sascalc.pdb'
        self.sas_paths = os.path.join('./', 'extract_large_file_number_0', 'run_1', 'sascalc','neutron_D2Op_100')
        self.option = 'single_frame'
        self.local_value = '100000'
  
        gui_mimic_extract_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'single_frame_small_peptide_sascalc.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm correct SAS files are chosen '''

        self.sas_paths = [x.strip() for x in self.sas_paths.split(',')]
        base_paths = []
        for sas_path in self.sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:
            outdirectory = os.path.join(self.runname, self.module, 'sascalc', base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_small_peptide_single_frame', base_path)
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_5(self):
        '''
        extract_trajectory and SAS: test dcd input, pdb output, periodic, SasCalc, >100k files
        '''

        self.pdb_filename = os.path.join(pdb_data_path, 'final_small_peptide.pdb')
        self.trajectory_filename = os.path.join(dcd_data_path, 'vmd_test_small_peptide.dcd')
        self.output_filename = 'periodic_small_peptide_sascalc.pdb'
        self.sas_paths = os.path.join('./', 'extract_large_file_number_0', 'run_1', 'sascalc','neutron_D2Op_100')
        self.option = 'sampling_frequency'
        self.local_value = '10000'
  
        gui_mimic_extract_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'periodic_small_peptide_sascalc.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm correct SAS files are chosen '''

        self.sas_paths = [x.strip() for x in self.sas_paths.split(',')]
        base_paths = []
        for sas_path in self.sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        for base_path in base_paths:
            outdirectory = os.path.join(self.runname, self.module, 'sascalc', base_path)
            correct_outdirectory = os.path.join(
                module_data_path, self.runname, self.module, 'sascalc_small_peptide_periodic', base_path)
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

                       

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

