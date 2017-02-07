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
import sassie.tools.merge_utilities.gui_mimic_merge_utilities as gui_mimic_merge_utilities
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

**************************
    merge trajectory 
**************************
    reference PDB           reference PDB           reference PDB           reference PDB        *reference PDB file
    input PDB               input DCD               input PDB               input DCD            *trajectory file
    output PDB              output PDB              output DCD              output DCD           *output file
                                        all

    reference PDB           reference PDB           reference PDB           reference PDB
    input PDB               input DCD               input PDB               input DCD
    output PDB              output PDB              output DCD              output DCD
                                        weight files

    reference PDB           reference PDB           reference PDB           reference PDB
    input PDB               input DCD               input PDB               input DCD
    output PDB              output PDB              output DCD              output DCD
                                         periodic

  
**************************
    extract SAS
**************************
   +SasCalc                 xtal2sas                cryson              crysol
                                       all

   +SasCalc                 xtal2sas                cryson              crysol
                                       weight files                                       

   +SasCalc                 xtal2sas                cryson              crysol
                                        periodic
                                    

    +SasCalc option can take sas_paths with subdirectories.  Two paths were tested, one with no subdirectories and one with 3 subdirectories.


**********************************
    extract trajectory and SAS
**********************************
    reference PDB           reference PDB           reference PDB           reference PDB
    input PDB               input DCD               input PDB               input DCD
    output PDB              output PDB              output DCD              output DCD
    sas path*               sas path*               sas path*               sas path*

    *sas path depends on SAS option:  SasCalc, xtal2sas, cryson, crysol

    selection options:  all, weight files, periodic

    NOTE:  Since all of these conditions are tested individually above, cryson and
           xtal2sas are retired, and crysol will be retired, the following two options were chosen:

    reference PDB           reference PDB
    input DCD               input DCD
    output PDB              output DCD
    SasCalc                 SasCalc
    weight files            all
    
    '''

    module = 'merge_utilities'

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
            if not are_dir_trees_equal(new_dir1, new_dir2):
                return False
        return True

    def test_1(self):
        '''
        merge_trajectory:  test dcd input, dcd output, all
        '''

        self.merge_option = '2' 
        self.merge_type_option = '0'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdb_file)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'all.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdb_file)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)
                        
    def test_2(self):
        '''
        merge_trajectory:  test dcd input, dcd output, weight files
        '''

        self.merge_option = '2'
        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt')
        self.output_filename = 'chosen_weights.dcd'  

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdb_file)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'chosen_weights.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdb_file)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

    def test_3(self):
        '''
        merge_trajectory:  test dcd input, dcd output, periodic
        '''

        self.merge_option = '2'
        self.merge_type_option = '2'
        self.local_value = '2'
        self.output_filename = 'periodic.dcd'  

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdb_file)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'periodic.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdb_file)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

    def test_4(self):
        '''
        merge_trajectory:  test dcd input, pdb output, all
        '''

        self.merge_option = '2' 
        self.merge_type_option = '0'
        self.output_filename = 'all.pdb'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'all.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_5(self):
        '''
        merge_trajectory:  test dcd input, pdb output, weight files
        '''

        self.merge_option = '2' 
        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt')        
        self.output_filename = 'chosen_weights.pdb'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'chosen_weights.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_6(self):
        '''
        merge_trajectory:  test dcd input, pdb output, periodic
        '''

        self.merge_option = '2' 
        self.merge_type_option = '2'
        self.local_value = '10'     
        self.output_filename = 'periodic.pdb'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'periodic.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_7(self):
        '''
        merge_trajectory:  test pdb input, dcd output, all
        '''

        self.merge_option = '2' 
        self.merge_type_option = '0'
        self.trajectory_names = os.path.join(pdb_data_path,'run_m1.pdb')+','+os.path.join(pdb_data_path,'run_m2.pdb')          

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdb_file)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'all_pdb_input.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdb_file)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

    def test_8(self):
        '''
        merge_trajectory:  test pdb input, dcd output, weight_files
        '''

        self.merge_option = '2' 
        self.merge_type_option = '1'
        self.trajectory_names = os.path.join(pdb_data_path,'run_m1.pdb')+','+os.path.join(pdb_data_path,'run_m2.pdb')
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt')
        self.output_filename = 'chosen_weights_pdb_input.dcd'                  

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdb_file)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'chosen_weights_pdb_input.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdb_file)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)
                        
    def test_9(self):
        '''
        merge_trajectory:  test pdb input, dcd output, periodic
        '''

        self.merge_option = '2' 
        self.merge_type_option = '2'
        self.trajectory_names = os.path.join(pdb_data_path,'run_m1.pdb')+','+os.path.join(pdb_data_path,'run_m2.pdb')
        self.local_value = '2'
        self.output_filename = 'periodic_pdb_input.dcd'                  

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdb_file)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'periodic_pdb_input.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdb_file)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

    def test_10(self):
        '''
        merge_trajectory:  test pdb input, pdb output, all
        '''

        self.merge_option = '2' 
        self.merge_type_option = '0'
        self.trajectory_names = os.path.join(pdb_data_path,'run_m1.pdb')+','+os.path.join(pdb_data_path,'run_m2.pdb')        
        self.output_filename = 'all_pdb_input.pdb'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'all_pdb_input.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_11(self):
        '''
        merge_trajectory:  test pdb input, pdb output, weight_files
        '''

        self.merge_option = '2' 
        self.merge_type_option = '1'
        self.trajectory_names = os.path.join(pdb_data_path,'run_m1.pdb')+','+os.path.join(pdb_data_path,'run_m2.pdb')
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt')                  
        self.output_filename = 'chosen_weights_pdb_input.pdb'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'chosen_weights_pdb_input.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_12(self):
        '''
        merge_trajectory:  test pdb input, pdb output, periodic
        '''

        self.merge_option = '2' 
        self.merge_type_option = '2'
        self.trajectory_names = os.path.join(pdb_data_path,'run_m1.pdb')+','+os.path.join(pdb_data_path,'run_m2.pdb')
        self.local_value = '10'                 
        self.output_filename = 'periodic_pdb_input.pdb'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'periodic_pdb_input.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_13(self):
        '''
        merge_sas  test SasCalc, all
        '''

        self.merge_option = '1'
        self.merge_type_option = '0'
        self.sas_type = '0'
        self.sas_paths = os.path.join(other_data_path,'merge_files_0','run_0','sascalc')+','+os.path.join(other_data_path, 'merge_files_0','run_1','sascalc')  
  
        gui_mimic_merge_utilities.run_module(self)

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
                module_data_path, self.runname, self.module, 'sascalc_all', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_14(self):
        '''
        merge_sas  test SasCalc, weight files
        '''

        self.merge_option = '1'
        self.merge_type_option = '1'
        self.sas_type = '0'
        self.sas_paths = os.path.join(other_data_path,'merge_files_0','run_0','sascalc')+','+os.path.join(other_data_path, 'merge_files_0','run_1','sascalc') 
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt') 
  
        gui_mimic_merge_utilities.run_module(self)

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
                module_data_path, self.runname, self.module, 'sascalc_weight_files', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_15(self):
        '''
        merge_sas  test SasCalc, periodic
        '''

        self.merge_option = '1'
        self.merge_type_option = '2'
        self.sas_type = '0'
        self.sas_paths = os.path.join(other_data_path,'merge_files_0','run_0','sascalc')+','+os.path.join(other_data_path, 'merge_files_0','run_1','sascalc') 
        self.local_value = '2'
  
        gui_mimic_merge_utilities.run_module(self)

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
                module_data_path, self.runname, self.module, 'sascalc_periodic', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_16(self):
        '''
        merge_sas  test xtal2sas, all
        '''

        self.merge_option = '1'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'xtal2sas')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_all')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_17(self):
        '''
        merge_sas  test xtal2sas, weight_files
        '''

        self.merge_option = '1'
        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt') 

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'xtal2sas')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_weight_files')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                                               
    def test_18(self):
        '''
        merge_sas  test xtal2sas, periodic
        '''

        self.merge_option = '1'
        self.merge_type_option = '2'
        self.local_value = '2'

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'xtal2sas')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'xtal2sas_periodic')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                                                                                              
    def test_19(self):
        '''
        merge_sas  test cryson, all
        '''

        self.merge_option = '1'
        self.sas_type = '2'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','cryson')+','+os.path.join(other_data_path, 'merge_files_0','run_1','cryson')        

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'cryson')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_all')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_20(self):
        '''
        merge_sas  test cryson, weight_files
        '''

        self.merge_option = '1'
        self.sas_type = '2'
        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt')         
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','cryson')+','+os.path.join(other_data_path, 'merge_files_0','run_1','cryson')        

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'cryson')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_weight_files')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_21(self):
        '''
        merge_sas  test cryson, periodic
        '''

        self.merge_option = '1'
        self.sas_type = '2'
        self.merge_type_option = '2'
        self.local_value = '2'   
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','cryson')+','+os.path.join(other_data_path, 'merge_files_0','run_1','cryson')        

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'cryson')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'cryson_periodic')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_22(self):
        '''
        merge_sas  test crysol, all
        '''

        self.merge_option = '1'
        self.sas_type = '3'
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','crysol')+','+os.path.join(other_data_path, 'merge_files_0','run_1','crysol')        

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'crysol')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_all')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_23(self):
        '''
        merge_sas  test crysol, weight_files
        '''

        self.merge_option = '1'
        self.sas_type = '3'
        self.merge_type_option = '1'
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt')         
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','crysol')+','+os.path.join(other_data_path, 'merge_files_0','run_1','crysol')        

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'crysol')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_weight_files')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_24(self):
        '''
        merge_sas  test crysol, periodic
        '''

        self.merge_option = '1'
        self.sas_type = '3'
        self.merge_type_option = '2'
        self.local_value = '2'   
        self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','crysol')+','+os.path.join(other_data_path, 'merge_files_0','run_1','crysol')        

        gui_mimic_merge_utilities.run_module(self)

        ''' confirm correct SAS files are in merged directories '''

        outdirectory = os.path.join(self.runname, self.module, 'crysol')
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'crysol_periodic')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_25(self):
        '''
        merge trajectory and SAS:  test dcd input, dcd output, all, Sascalc
        '''

        self.merge_option = '0' 
        self.merge_type_option = '0'
        self.output_filename = 'all_dcd_input_sascalc.dcd'
        self.sas_type = '0'
        self.sas_paths = os.path.join(other_data_path,'merge_files_0','run_0','sascalc')+','+os.path.join(other_data_path, 'merge_files_0','run_1','sascalc')  
  
        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdb_file)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'all_dcd_input_sascalc.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdb_file)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)                  

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
                module_data_path, self.runname, self.module, 'sascalc_all_dcd_input', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_26(self):
        '''
        merge trajectory and SAS:  test dcd input, pdb output, weight_files, Sascalc
        '''

        self.merge_option = '0' 
        self.merge_type_option = '1'
        self.output_filename = 'weight_files_dcd_input_sascalc.pdb'
        self.sas_type = '0'
        self.local_value = os.path.join(other_data_path,'weights_file_m1.txt')+','+os.path.join(other_data_path,'weights_file_m2.txt')         
        self.sas_paths = os.path.join(other_data_path,'merge_files_0','run_0','sascalc')+','+os.path.join(other_data_path, 'merge_files_0','run_1','sascalc')  
  
        gui_mimic_merge_utilities.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.output_filename)
        correct_outfile = os.path.join(
            module_data_path, 'weight_files_dcd_input_sascalc.pdb')
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
                module_data_path, self.runname, self.module, 'sascalc_weight_files_dcd_input', base_path)
#            print 'correct outdirectory: ', correct_outdirectory
            assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                                                                                            
        
    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

