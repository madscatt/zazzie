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
import sassie.simulate.monomer_monte_carlo.gui_mimic_monomer_monte_carlo as gui_mimic_monomer_monte_carlo
#import gui_mimic_monomer_monte_carlo as gui_mimic_monomer_monte_carlo

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'simulate', 'monomer_monte_carlo') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Monomer_Monte_Carlo(MockerTestCase):

    '''
    DIHEDRAL is the module that contains the functions
    that are used to generate ensembles of structures by varying
    protein dihedral angles.


    This module calls to C / Python extension modules to speed up
    calculations (see: overlap.c).

    REFERENCE:

    J. A. D. MacKerell et al.
    Journal of Physical Chemistry B,  102  3586-3616  (1998)

    B. R. Brooks et al.
    Journal of Computational Chemistry  4  187--217  (1983)


    INPUT:	
		PDB file with molecular structure data 
        number of trial attempts
        number of failed trial attempts before returning to previously-accepted structure
        temperature (K)
        molecule type (protein or rna)
        number of flexible regions to vary
        residue range for each flexible region (pairs of hyphenated integers separated by commas; converted to original reslow and numcont format)
        maximum angle that each torsion can sample in a single move for each flexible region
        structure alignment range (pairs of hyphenated integers separated by commas; converted to original lowres1 and highres1 format)
        overlap basis (heavy, backbone, all or user's choice of specific atom)

        Advanced input:

        low Rg cutoff
        high Rg cutoff
        directed Monte Carlo (0 or Rg value)
        Z cutoff value
        atomic constraints file
        

    OUTPUT:
        original PDB file with molecular structure data
        DCD file containing accepted structures
        file containing Rg values for all trial structures
        file containing Rg value for accepted structures
        file containing run statistics

 
    Use cases:

    1.  molecule type
        a.  protein
        b.  rna

    2.  overlap basis
        a.  heavy
        b.  backbone
        c.  all
        d.  individual atom (protein only; rna needs more than one basis atom to avoid overlap)

    3.  flexible regions
        a.  one flexible region
        b.  more than one flexible region

    4.  advanced input options
        a.  low Rg cutoff
        b.  high Rg cutoff
        c.  directed Monte Carlo
        d.  Z-cutoff
        e.  contraints

    Inputs tested:

    runname:        string      project name                          
    path:           string      input file path                 
    dcdfile:        string      name of output dcd file containing accepted structures       
    moltype:        string      molecule type ('protein' or 'rna')
    pdbfile:        string      name of input pdb file containing intial structure
    trials:         integer     number of Monte Carlo move attempts
    goback:         integer     number of failed Monte Carlo attempts before returning to previously accepted structure
    temp            float       run temperature (K)
    numranges       integer     number of flexible regions
    dtheta          float_array maximum angle that torsion can sample (in each flexible region)
    reslow          int_array   low residue number for each flexible region
    numcont         int_array   number of contiguous residues per flexible region (not enetered directly; parsed from entered residue range in GenApp)
    lowres1         integer     low residue for (non-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
    highres1        integer     high residue for (no-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
    basis           string      type of basis for overlap check ("all", "heavy", "backbone" or specific atom name, i.e., "CA")
    cutoff          float       overlap cutoff distance (all=0.8, heavy=0.8, backbone=1.0, specific atom=user's choice)
    lowrg           float       low Rg cutoff value if Advanced Input is chosen
    highrg          float       high Rg cutoff value if Advanced Input is chosen
    zflag           integer     enable zcutoff flag (0=no, 1=yes)
    zcutoff         float       zcutoff value (discard structures with any z-axis coordinates less than this value)
    cflag           integer     enable atomic constraint flag (0=no, 1=yes)
    confile         string      name of file describing additional constraints to check before accepting a structure
    directedmc      float       non-zero Rg value to guide Monte Carlo run; 0=no directed Monte Carlo (used if Advanced Input is chosen)

    Inputs not tested (options not yet implemented):

    nonbondflag     integer     flag for nonbonded option
    nonbondedscale  float       nonbondedscale value
    psffilepath     string      path to psf file
    psffilename     string      psf file name
    parmfilepath    string      path to CHARMM parameter file
    parmfilename    string      name of CHARMM parameter file
    plotflag        integer     option to plot structure number vs Rg


    Test tree:

    project name
    input PDB file
    output DCD file

******************************
    single flexible region
******************************
    protein             protein             protein             protein             rna             rna             rna                      
    heavy               backbone            all*                single atom         heavy           backbone       all*    

                                                            no advanced options 
                                                            low Rg cutoff
                                                            high Rg cutoff
                                                            directed MC
                                                            constraints
                                                            Z-cutoff
 
    *only run with no advanced options

**********************************
    multiple flexible regions
**********************************
 
    NOTE:  Since the basis atom options and advanced inputs were tested individually above, the following three options were chosen:

    protein             protein             rna
    5 regions           5 regions           2 regions
    heavy               heavy               heavy
    low Rg cutoff       Z-cutoff            no advanced options
    high Rg cutoff
    directed MC
    constraints

    '''

    module = 'monomer_monte_carlo'

    def setUp(self):

       gui_mimic_monomer_monte_carlo.test_variables(self, paths)


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
        ignoring files that have date stamps: self.dcdfile, *.sassie_json
        and *.sassie_log
        '''
        dirs_cmp = filecmp.dircmp(dir1, dir2,[self.dcdfile])    
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


    def test_1(self):
        '''
        protein, one flexible region, heavy atom basis, no advanced options
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'    
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_2(self):
        '''
        protein, one flexible region, heavy atom basis, low Rg cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1' 
        self.lowrg = '23.0'   
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_3(self):
        '''
        protein, one flexible region, heavy atom basis, high Rg cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1' 
        self.highrg = '24.0'   
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_4(self):
        '''
        protein, one flexible region, heavy atom basis, directed MC
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1' 
        self.directedmc = '25.0' 
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_5(self):
        '''
        protein, one flexible region, heavy atom basis, constraints
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1' 
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'ma_constraints.txt')
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_6(self):
        '''
        protein, one flexible region, heavy atom basis, Z-cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1' 
        self.zflag = '1'
        self.zcutoff = '25'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_heavy_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_7(self):
        '''
        protein, one flexible region, backbone basis, no advanced options
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'backbone'
        self.cutoff = '1.0'    
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_8(self):
        '''
        protein, one flexible region, backbone basis, low Rg cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'backbone'
        self.cutoff = '1.0' 
        self.lowrg = '23.0'    
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_9(self):
        '''
        protein, one flexible region, backbone basis, high Rg cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'backbone'
        self.cutoff = '1.0' 
        self.highrg = '24.0'    
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_10(self):
        '''
        protein, one flexible region, backbone basis, directed MC
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'backbone'
        self.cutoff = '1.0' 
        self.directedmc = '25.0'    
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_11(self):
        '''
        protein, one flexible region, backbone basis, constraints
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'backbone'
        self.cutoff = '1.0' 
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'ma_constraints.txt')    
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_12(self):
        '''
        protein, one flexible region, backbone basis, Z-cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'backbone'
        self.cutoff = '1.0' 
        self.zflag = '1'
        self.zcutoff = '25'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_backbone_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_13(self):
        '''
        protein, one flexible region, all basis, no advanced options
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'all'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_all_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_all_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_14(self):
        '''
        protein, one flexible region, atom basis, no advanced options
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'CA'
        self.cutoff = '3.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_15(self):
        '''
        protein, one flexible region, atom basis, low Rg cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'CA'
        self.cutoff = '3.0'
        self.lowrg = '23.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_16(self):
        '''
        protein, one flexible region, atom basis, high Rg cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'CA'
        self.cutoff = '3.0'
        self.highrg = '24.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_17(self):
        '''
        protein, one flexible region, atom basis, directed MC
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'CA'
        self.cutoff = '3.0'
        self.directedmc = '25.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_18(self):
        '''
        protein, one flexible region, atom basis, constraints
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'CA'
        self.cutoff = '3.0'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'ma_constraints.txt')         
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_19(self):
        '''
        protein, one flexible region, atom basis, Z-cutoff
        '''

        self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
        self.reslow = '115'
        self.numcont = '24'
        self.lowres1 = '1'
        self.highres1 = '110'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.basis = 'CA'
        self.cutoff = '3.0'
        self.zflag = '1'
        self.zcutoff = '25'      
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_atom_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_20(self):
        '''
        rna, one flexible region, heavy basis, no advanced options
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_21(self):
        '''
        rna, one flexible region, heavy basis, low Rg cutoff
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.lowrg = '26.5'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_22(self):
        '''
        rna, one flexible region, heavy basis, high Rg cutoff
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.highrg = '27.2'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_23(self):
        '''
        rna, one flexible region, heavy basis, directed MC
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.directedmc = '20.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_24(self):
        '''
        rna, one flexible region, heavy basis, constraints
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'rna_constraints.txt')
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_25(self):
        '''
        rna, one flexible region, heavy basis, Z-cutoff
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.zflag = '1'
        self.zcutoff = '-25.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_heavy_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_26(self):
        '''
        rna, one flexible region, backbone basis, no advanced options
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.basis = 'backbone'
        self.cutoff = '1.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_27(self):
        '''
        rna, one flexible region, backbone basis, low Rg cutoff
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.basis = 'backbone'
        self.cutoff = '1.0'
        self.lowrg = '26.5'        
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                                      
    def test_28(self):
        '''
        rna, one flexible region, backbone basis, high Rg cutoff
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.basis = 'backbone'
        self.cutoff = '1.0'
        self.highrg = '27.2'        
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_29(self):
        '''
        rna, one flexible region, backbone basis, directed MC
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.basis = 'backbone'
        self.cutoff = '1.0'
        self.directedmc = '20.0'        
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_30(self):
        '''
        rna, one flexible region, backbone basis, constraints
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.basis = 'backbone'
        self.cutoff = '1.0'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'rna_constraints.txt')   
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_31(self):
        '''
        rna, one flexible region, backbone basis, Z-cutoff
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.basis = 'backbone'
        self.cutoff = '1.0'
        self.zflag = '1'
        self.zcutoff = '-25.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_backbone_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_32(self):
        '''
        rna, one flexible region, all basis, no advanced options
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.reslow = '24'
        self.numcont = '7'
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0'
        self.numranges = '1'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        self.basis = 'all'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_all_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_all_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_33(self):
        '''
        rna, two flexible regions, heavy basis, no advanced options
        '''

        self.moltype = 'rna'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.numranges = '2'
        self.reslow = '24,47'
        self.numcont = '7,9' 
        self.lowres1 = '34'
        self.highres1 = '45'
        self.dtheta = '30.0,30.0'
        self.seed = '1,321'
        self.trials = '8'
        self.goback = '2'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_34(self):
        '''
        protein, five flexible regions, heavy basis, Z-cutoff
        '''

        self.zflag = '1'
        self.zcutoff = '-35.0'
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_5_heavy_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_5_heavy_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                          
    def test_35(self):
        '''
        protein, five flexible regions, heavy basis, all options except Z-cutoff
        '''

        self.lowrg = '44.0'
        self.highrg = '48.0'
        self.directedmc = '46.0'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'monomer_mc_constraints.txt')
        gui_mimic_monomer_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_5_heavy_nozcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_5_heavy_nozcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

          
                                                                                  


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

