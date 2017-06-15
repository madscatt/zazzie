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
import sassie.simulate.energy_minimization.gui_mimic_energy_minimization as gui_mimic_energy_minimization
#import gui_mimic_energy_minimization as gui_mimic_energy_minimization

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'simulate', 'energy_minimization') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Energy_Minimization(MockerTestCase):

    '''
        MINIMIZE is the module that contains the functions
        that are used to run a series of energy minimization calculations
	    on a set of structures in a supplied pdb/dcd file.

        This module is called from Structure Minization in the main 
        GUI through the graphical_minimize.py script.

        REFERENCE:

        J. C. Phillips et al.
        Journal of Computational Chemistry  26  1781-1802  (2005)

        INPUT:  variable descriptions:
              
        reference PDB file
        input PDB or DCD file
        number of minimization steps
        path and name of topology file
        input PSF file
        output (DCD) file
        number of CPUs
        keep output files (0==no, 1==yes)
        MD flag (0=min, 1=min+md, 2=min+md+min)
        number of MD steps (if md = 1 or 2)
        solvent dielectric constant
        temperature (K)
        frequency to save individual DCD files
        flag to use external input file (True or False)
        external input file name
        velocity restart file name
        extended system restart file name


        Advanced input:

        name of user-supplied CHARMM parameter file

        
        OUTPUT:
                
        files stored in ~/run_name/energy_minimization directory:
                
        original PDB file with molecular structure data
        original PSF file
        DCD file containing energy minimized (and md, if applicable) trajectories (size depends on DCD write frequency)
        PDB file containing final energy minimized (and md, if applicable) structure
        output file (if keep output files is chosen)
        input file for energy minimization 


     Use cases:

    1.  input:
        a.  input file is a PDB file
        b.  input file is a DCD file

    2.  MD option:
        a.  minimization only
        b.  minimization + MD
        c.  minimization + MD + minimization

    3.  other options:
        a.  keep output files
        b.  use external input file (new simulation)
        c.  use external input file (continued simulation) 
        

    Inputs tested:

	runname:                        string      run_name 
    infile:                         string      input pdb or dcd file name
    pdbfile:                        string      input (reference) pdb file name
    outfile:                        string      output dcd file name
    nsteps:                         integer     number of minimization steps
    resparmfile:                    string      path and name of topology file
    psffile: 		                string      psf file name
    ncpu	:                           integer     number of cpus to use
    keepout:                        integer     keep output files (0==no, 1==yes)
    dcdfreq:                        integer     save individual dcd frequency
    md                              integer     md flag (0=min, 1=min+md, 2=min+md+min)
    mdsteps                         integer     number of md steps
    dielect                         float       solvent dielectric constant
    temperature                     float       temperature
    use_external_input_file         boolean     flag to use external input file (True or False)
    external_input_file             string      external input file name
    velocity_restart_file           string      velocity restart file name
    extended_system_restart_file    string      extended system restart file name
    
    Inputs not tested:

    path                            string      input path (not in variables)
    infiletype                      string      input file type (pdb or dcd)    #determined by file name suffix
    charmm_parameter_file           string      name of user-supplied CHARMM parameter file

    Test tree:

    project name

*********************************************
    default input file (generated by program)
*********************************************

    reference PDB           reference PDB           reference PDB       reference PDB           *reference PDB file
    input PDB               input PDB               input DCD           input DCD               *trajectory file
    output DCD              output DCD              output DCD          output DCD              *output file
    minimization only       minimization only       minimization only   minimization only       
    1 processor             more than 1 processor   1 processor         more than 1 processor

                                           keep output files
                                           do not keep output files

**********************************************
    external input file (provided by user)
**********************************************
    reference PDB           reference PDB           reference PDB       reference PDB           *reference PDB file
    input PDB               input DCD               input DCD           input DCD               *trajectory file
    output DCD              output DCD              outpt DCD           output DCD              *output file
    minimization only       minimization only       minimization only   minimization only
    1 processor             more than 1 processor   1 processor         more than 1 processor

                                            keep output files
                                            do not keep output files

**********************************************
    MD options 1 & 2 
**********************************************
 
    NOTE:  Since the input options were tested above, the following options were chosen:

    reference PDB               reference PDB                   reference PDB
    input PDB                   input PDB                       input DCD (multi-frame)
    output DCD                  output DCD                      output DCD
    minimization + MD           minimization + MD               minimization + MD + minimization
    1 processor                 1 processor                     more than 1 processor
    keep output files           keep output files               keep output files
    external input file*        external input file*            default input file**
                                velocity restart file
                                extended system restart file

*includes a seed in the input file so the output dcd and pdb files will be the same each time it is run (1 processor must be used)
**output dcd and pdb files will not be the same; testing for completion and correct number of frames in output DCD file only

************************************************
    detect error in simulation
************************************************
    NOTE:  Determine if error in simulation can be detected by checking for junk.out file.

    reference PDB
    input PDB
    output DCD
    minimization + MD
    1 processor
    keep output files
    external input file*** 
    velocity restart file
    extended system restart file

*** external input file contains a temperature, which isn't compatible with the velocity restart file, so an error message appears in junk.out 
    '''


    module = 'energy_minimization'

    def setUp(self):

       gui_mimic_energy_minimization.test_variables(self, paths)


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
        ignoring self.outfile and min_00001.out, which have date stamps,
        temp.inp, which has a specific path that won't match the test path,
        and *.sassie_json and *.sassie_.log
        '''
        dirs_cmp = filecmp.dircmp(dir1, dir2,[self.outfile, 'min_00001.out', 'temp.inp'])
#        dirs_cmp = filecmp.dircmp(dir1, dir2)
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

    def check_for_completion(self):
        '''
        tests for completion of run without an error
        looks for junk.out file in current directory, which won't exist unless there was an error
        this test will be used when md=1 or 2 and the default input file is used (since the output files
        will not be the same and can't be compared)
        '''

        test = os.path.isfile('junk.out')
        if(test):
            return False
        return True


    def test_1(self):
        '''
        PDB input, minimization only, don't keep output files
        '''

        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

                                                                                  
    def test_2(self):
        '''
        DCD input, minimization only, don't keep output files
        '''

        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

  
    def test_3(self):
        '''
        PDB input, minimization only, keep output files
        '''

        self.keepout = '1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_4(self):
        '''
        DCD input, minimization only, keep output files
        '''

        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        self.keepout='1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_5(self):
        '''
        PDB input, minimization only, 1 processor, don't keep output files
        '''

        self.ncpu = '1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_1', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_1')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

                                                                                  
    def test_6(self):
        '''
        DCD input, minimization only, 1 processor, don't keep output files
        '''

        self.ncpu = '1'
        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_1', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_1')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

  
    def test_7(self):
        '''
        PDB input, minimization only, 1 processor, keep output files
        '''

        self.ncpu = '1'
        self.keepout = '1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep_1', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep_1')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_8(self):
        '''
        DCD input, minimization only, 1 processor, keep output files
        '''

        self.ncpu = '1'
        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        self.keepout='1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep_1', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep_1')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_9(self):
        '''
        PDB input, minimization only, don't keep output files, external input file
        '''

        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

                                                                                  
    def test_10(self):
        '''
        DCD input, minimization only, don't keep output files, external input file
        '''

        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

  
    def test_11(self):
        '''
        PDB input, minimization only, keep output files, external input file
        '''

        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        self.keepout = '1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_12(self):
        '''
        DCD input, minimization only, keep output files, external input file
        '''

        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        self.keepout='1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_13(self):
        '''
        PDB input, minimization only, 1 processor, don't keep output files, external input file
        '''

        self.ncpu = '1'
        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_1_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_1_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

                                                                                  
    def test_14(self):
        '''
        DCD input, minimization only, 1 processor, don't keep output files, external input file
        '''

        self.ncpu = '1'
        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_1_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_1_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

  
    def test_15(self):
        '''
        PDB input, minimization only, 1 processor, keep output files, external input file
        '''

        self.ncpu = '1'
        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        self.keepout = '1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep_1_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_keep_1_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_16(self):
        '''
        DCD input, minimization only, 1 processor, keep output files, external input file
        '''

        self.ncpu = '1'
        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input.inp')
        self.infile = os.path.join(dcd_data_path,'ten_mer.dcd')
        self.keepout='1'
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep_1_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_min_keep_1_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_17(self):
        '''
        PDB input, minimization + MD, 1 processor, keep output files, external input
        '''

        self.md='1'
        self.ncpu='1'
        self.keepout='1'
        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input_1.inp')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_md_keep_1_ext', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_md_keep_1_ext')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_18(self):
        '''
        PDB input, minimization + MD, 1 processor, keep output files, external input, velocity restart and extended system restart files
        '''

        self.md='1'
        self.ncpu='1'
        self.keepout='1'
        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input_2.inp')
        self.velocity_restart_file = os.path.join(other_data_path,'velocity_restart.vel')
        self.extended_system_restart_file = os.path.join(other_data_path,'extended_system_restart.xsc')
        gui_mimic_energy_minimization.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_md_keep_1_ext_vel_xsc', self.outfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_min_md_keep_1_ext_vel_xsc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_19(self):
        '''
        2-frame DCD input, minimization + MD + minimization, keep output files
        Test for completion and proper number of frames in output dcd file
        '''

        self.md='2'
        self.keepout='1'
        self.infile = os.path.join(dcd_data_path,'ten_mer2.dcd')
        gui_mimic_energy_minimization.run_module(self)

        '''check for completion'''
        assert_equals(self.check_for_completion(),True)

        '''check for proper number of frames in output DCD file'''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        dcdfile = molecule.open_dcd_read(outfile)
        nf = dcdfile[2]
#        print 'number of frames: ', nf
        assert_equals(nf,2)       

    def test_20(self):
        '''
        PDB input, minimization + MD, 1 processor, keep output files, external input, velocity restart and extended system restart files
        Determine if error in simulation can be detected
        If a velocity restart file is used, a temperature cannot also be given.  Thus, an error will occur during execution and an error
        message will appear in junk.out.
        '''

        self.md='1'
        self.ncpu='1'
        self.keepout='1'
        self.use_external_input_file = 'True'
        self.external_input_file = os.path.join(other_data_path,'external_input_1.inp')
        self.velocity_restart_file = os.path.join(other_data_path,'velocity_restart.vel')
        self.extended_system_restart_file = os.path.join(other_data_path,'extended_system_restart.xsc')
        gui_mimic_energy_minimization.run_module(self)

        '''check for completion error'''
        assert_equals(self.check_for_completion(),False)

        try:
            os.system('rm -f junk.out ')
        except:
            pass
        try:
            os.system('rm -f temp.inp ')
        except:
            pass              
  
    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

