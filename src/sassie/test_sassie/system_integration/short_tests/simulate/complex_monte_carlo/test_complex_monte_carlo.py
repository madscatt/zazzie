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

import sassie.sasmol.sasmol as sasmol
import sassie.simulate.monte_carlo.complex.gui_mimic_complex_monte_carlo as gui_mimic_complex_monte_carlo
#import gui_mimic_complex_monte_carlo as gui_mimic_complex_monte_carlo

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'simulate', 'complex_monte_carlo') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Complex_Monte_Carlo(MockerTestCase):

    '''
	NMR_DIHEDRAL is the module that contains the functions
	that are used to generate ensembles of structures by varying
	protein dihedral angles.  This particular version allows multiple
	flexible proteins in the presence of non-flexible proteins and
	nucleic acids.

	This module calls to C / Python extension modules to speed up
	calculations.


    INPUT:	
		PDB file with molecular structure data 
        number of trial attempts
        number of failed trial attempts before returning to previously-accepted structure
        temperature (K)
        total number of segments
        number of flexible segments
        overlap basis (heavy, backbone, all or user's choice of specific atom)

        For each flexible segment:
        molecule type (protein or rna)
        segment name
        number of flexible regions to vary
        residue range for each flexible region (pairs of hyphenated integers separated by commas; converted to original reslow and numcont format)
        maximum angle that each torsion can sample in a single move for each flexible region
        structure alignment range (pairs of hyphenated integers separated by commas; converted to original lowres1 and highres1 format)

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
        a.  one flexible segment
            i.  one flexible region
            ii. more than one flexible region
        b.  more than one flexible segment
            i.  one flexible region per segment
            ii. more than one flexible region per segment

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
    pdbfile:        string      name of input pdb file containing intial structure
    trials:         integer     number of Monte Carlo move attempts
    goback:         integer     number of failed Monte Carlo attempts before returning to previously accepted structure
    temp:           float       run temperature (K)
    nsegments:      integer     total number of segments   
    npsegments:     integer     number of segments containing flexible regions
    flpsegname:     string      names of segments with flexible regions (separated by commas if more than one)
    segbasis:       string      type of basis for overlap check ("all", "heavy", "backbone" or specific atom name, i.e., "CA")     
    seglow:         integer     low residue for (non-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
    seghigh:        integer     high residue for (no-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
    lowrg:          float       low Rg cutoff value if Advanced Input is chosen
    highrg:         float       high Rg cutoff value if Advanced Input is chosen
    zflag:          integer     enable zcutoff flag (0=no, 1=yes)
    zcutoff:        float       zcutoff value (discard structures with any z-axis coordinates less than this value)
    cflag:          integer     enable atomic constraint flag (0=no, 1=yes)
    confile:        string      name of file describing additional constraints to check before accepting a structure
    directedmc:     float       non-zero Rg value to guide Monte Carlo run; 0=no directed Monte Carlo (used if Advanced Input is chosen)

    psegvariables:              flexible segment variables
                    integer     number of flexible regions
                    float_array maximum angle that torsion can sample (in each flexible region)
                    int_array   low residue number for each flexible region
                    int_array   number of contiguous residues per flexible region (not enetered directly; parsed from entered residue range in GenApp)
                    string      molecule type ('protein' or 'rna')                                       


    Inputs not tested (options not currently implemented):

    psffilepath     string      path to psf file
    psffilename     string      psf file name
    parmfilepath    string      path to CHARMM parameter file
    parmfilename    string      name of CHARMM parameter file
    plotflag        integer     option to plot structure number vs Rg



    Test tree:

    project name
    input PDB file
    output DCD file

****************************************************
    non-flexible protein segment +
    single flexible segment with one flexible region
****************************************************
    protein            protein             protein             protein             rna             rna             rna                      
    heavy              backbone            all*                single atom         heavy           backbone       all*    

                                                            no advanced options 
                                                            low Rg cutoff
                                                            high Rg cutoff
                                                            directed MC
                                                            constraints
                                                            Z-cutoff
 
    *only run with no advanced options


**********************************************************
    non-flexible protein segment +
    single flexible segment with multiple flexible regions
**********************************************************
    protein            protein             protein             protein             rna             rna             rna                      
    heavy              backbone            all*                single atom         heavy           backbone       all*    

                                                            no advanced options 
                                                            low Rg cutoff
                                                            high Rg cutoff
                                                            directed MC
                                                            constraints
                                                            Z-cutoff
 
    *only run with no advanced options

*******************************************************
    multiple flexible segments with one flexible region
*******************************************************
 
    NOTE:  Since the basis atom options and advanced inputs were tested individually above, the following two options were chosen:

    protein:rna                     protein:rna                     
    6 flexible protein segments     6 flexible protein segments     
    non-flexible rna segment        flexible rna segment                           
    heavy                           backbone                           
    low Rg cutoff                   no advanced options                        
    high Rg cutoff
    directed MC

    '''

    module = 'complex_monte_carlo'

    def setUp(self):

       gui_mimic_complex_monte_carlo.test_variables(self, paths)


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
        ignoring self.dcdfile, which has a date stamp
        '''
        dirs_cmp = filecmp.dircmp(dir1, dir2,[self.dcdfile])
#        dirs_cmp = filecmp.dircmp(dir1, dir2)
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
        protein:protein, one flexible region, heavy atom basis, no advanced options
        '''

        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, heavy atom basis, low Rg cutoff
        '''

        self.lowrg = '42.0'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, heavy atom basis, high Rg cutoff
        '''

        self.highrg = '47.0'   
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, heavy atom basis, directed MC
        '''

        self.directedmc = '45.0' 
        gui_mimic_complex_monte_carlo.run_module(self)

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

        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'pai_vn_constraints.txt')
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, heavy atom basis, Z-cutoff
        '''

        self.segbasis = 'heavy'
        self.zflag = '1'
        self.zcutoff = '-18.0'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein one flexible region, backbone basis, no advanced options
        '''

        self.segbasis = 'backbone'   
        gui_mimic_complex_monte_carlo.run_module(self)

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

        self.segbasis = 'backbone' 
        self.lowrg = '42.0'    
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, backbone basis, high Rg cutoff
        '''

        self.segbasis = 'backbone'
        self.highrg = '47.0'    
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, backbone basis, directed MC
        '''

        self.segbasis = 'backbone'
        self.cutoff = '1.0' 
        self.directedmc = '45.0'    
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, backbone basis, constraints
        '''

        self.segbasis = 'backbone'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'pai_vn_constraints.txt')    
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, backbone basis, Z-cutoff
        '''

        self.segbasis = 'backbone'
        self.zflag = '1'
        self.zcutoff = '-19.0'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, all basis, no advanced options
        '''

        self.segbasis = 'all'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, atom basis, no advanced options
        '''

        self.segbasis = 'CA,CA'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, atom basis, low Rg cutoff
        '''

        self.segbasis = 'CA,CA'
        self.lowrg = '42.0'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, atom basis, high Rg cutoff
        '''

        self.segbasis = 'CA,CA'
        self.highrg = '47.0'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, atom basis, directed MC
        '''

        self.segbasis = 'CA,CA'
        self.directedmc = '45.0'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, atom basis, constraints
        '''

        self.segbasis = 'CA,CA'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'pai_vn_constraints.txt')         
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, one flexible region, atom basis, Z-cutoff
        '''

        self.segbasis = 'CA,CA'
        self.zflag = '1'
        self.zcutoff = '-19.0'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, heavy atom basis, no advanced options
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, heavy atom basis, low Rg cutoff
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.lowrg = '46.0'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, heavy atom basis, high Rg cutoff
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '50.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, heavy atom basis, directed MC
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.directedmc = '47.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, heavy atom basis, constraints
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'rna_protein_constraints.txt')
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, heavy atom basis, Z-cutoff
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.zflag = '1'
        self.zcutoff = '-60.0' 
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, backbone basis, no advanced options
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.segbasis = 'backbone'
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, backbone basis, low Rg cutoff
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.lowrg = '46.0'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.segbasis = 'backbone'
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
            module_data_path, self.runname, self.module, 'rna_heavy_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
        
    def test_28(self):
        '''
        protein:rna one flexible region, backbone basis, high Rg cutoff
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '50.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.segbasis = 'backbone'
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, backbone basis, directed MC
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.directedmc = '47.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.segbasis = 'backbone'
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, backbone basis, constraints
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'rna_protein_constraints.txt')
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.segbasis = 'backbone'
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, backbone basis, Z-cutoff
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.zflag = '1'
        self.zcutoff = '-60.0' 
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.segbasis = 'backbone'
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:rna one flexible region, all basis, no advanced options
        '''

        self.trials = '10'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.segbasis = 'all'
        self.flpsegname = 'RNA1'
        self.psegvariables = [['1', '30', '128', '2', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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
        protein:protein, multiple flexible regions, heavy atom basis, no advanced options
        '''

        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_34(self):
        '''
        protein:protein, multiple flexible regions, heavy atom basis, low Rg cutoff
        '''

        self.lowrg = '44.0'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]        
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_35(self):
        '''
        protein:protein, multiple flexible regions, heavy atom basis, high Rg cutoff
        '''

        self.highrg = '46.7'   
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
        
    def test_36(self):
        '''
        protein:protein, multiple flexible regions, heavy atom basis, directed MC
        '''

        self.directedmc = '45.0'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']] 
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_37(self):
        '''
        protein, multiple flexible regions, heavy atom basis, constraints
        '''

        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'pai_vn_constraints.txt')
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_38(self):
        '''
        protein:protein, multiple flexible regions, heavy atom basis, Z-cutoff
        '''

        self.zflag = '1'
        self.zcutoff = '-18.0'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_heavy_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_39(self):
        '''
        protein:protein multiple flexible regions, backbone basis, no advanced options
        '''

        self.segbasis = 'backbone'   
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_40(self):
        '''
        protein, multiple flexible regions, backbone basis, low Rg cutoff
        '''

        self.segbasis = 'backbone' 
        self.lowrg = '44.0'    
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_41(self):
        '''
        protein:protein, multiple flexible regions, backbone basis, high Rg cutoff
        '''

        self.segbasis = 'backbone'
        self.highrg = '46.7'    
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_42(self):
        '''
        protein:protein, multiple flexible regions, backbone basis, directed MC
        '''

        self.segbasis = 'backbone'
        self.cutoff = '1.0' 
        self.directedmc = '45.0'    
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_43(self):
        '''
        protein:protein, multiple flexible regions, backbone basis, constraints
        '''

        self.segbasis = 'backbone'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'pai_vn_constraints.txt')    
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_44(self):
        '''
        protein:protein, multiple flexible regions, backbone basis, Z-cutoff
        '''

        self.segbasis = 'backbone'
        self.zflag = '1'
        self.zcutoff = '-19.0'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_backbone_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_45(self):
        '''
        protein:protein, multiple flexible regions, atom basis, no advanced options
        '''

        self.segbasis = 'CA,CA'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_46(self):
        '''
        protein:protein, multiple flexible regions, atom basis, low Rg cutoff
        '''

        self.segbasis = 'CA,CA'
        self.lowrg = '44.0'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_47(self):
        '''
        protein:protein, multiple flexible regions, atom basis, high Rg cutoff
        '''

        self.segbasis = 'CA,CA'
        self.highrg = '46.7'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_48(self):
        '''
        protein:protein, multiple flexible regions, atom basis, directed MC
        '''

        self.segbasis = 'CA,CA'
        self.directedmc = '45.0'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
               
    def test_49(self):
        '''
        protein:protein, multiple flexible regions, atom basis, constraints
        '''

        self.segbasis = 'CA,CA'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'pai_vn_constraints.txt')         
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_50(self):
        '''
        protein:protein, multiple flexible regions, atom basis, Z-cutoff
        '''

        self.segbasis = 'CA,CA'
        self.zflag = '1'
        self.zcutoff = '-19.0'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_atom_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_51(self):
        '''
        protein:protein, multiple flexible regions, all basis, no advanced options
        '''

        self.segbasis = 'all'
        self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_all_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_2_all_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_52(self):
        '''
        protein:rna multiple flexible regions, heavy atom basis, no advanced options
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

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

    def test_53(self):
        '''
        protein:rna multiple flexible regions, heavy atom basis, low Rg cutoff
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.lowrg = '46.0'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
        
    def test_54(self):
        '''
        protein:rna multiple flexible regions, heavy atom basis, high Rg cutoff
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '56.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_55(self):
        '''
        protein:rna multiple flexible regions, heavy atom basis, directed MC
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.directedmc = '47.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                
    def test_56(self):
        '''
        protein:rna multiple flexible regions, heavy atom basis, constraints
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'rna_protein_constraints2.txt')
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                
    def test_57(self):
        '''
        protein:rna multiple flexible regions, heavy atom basis, Z-cutoff
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.zflag = '1'
        self.zcutoff = '-60.0' 
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_heavy_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_58(self):
        '''
        protein:rna multiple flexible regions, backbone basis, no advanced options
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.segbasis = 'backbone'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_59(self):
        '''
        protein:rna multiple flexible regions, backbone basis, low Rg cutoff
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.segbasis = 'backbone'
        self.seglow = '20'
        self.seghigh = '30'
        self.lowrg = '46.0'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_lowrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_lowrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
        
    def test_60(self):
        '''
        protein:rna multiple flexible regions, backbone basis, high Rg cutoff
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.segbasis = 'backbone'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '56.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_highrg', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_highrg')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_61(self):
        '''
        protein:rna multiple flexible regions, backbone basis, directed MC
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.segbasis = 'backbone'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.directedmc = '47.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                
    def test_62(self):
        '''
        protein:rna multiple flexible regions, backbone basis, constraints
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.segbasis = 'backbone'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.cflag = '1'
        self.confile = os.path.join(other_data_path,'rna_protein_constraints2.txt')
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_constraints', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_constraints')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                
    def test_63(self):
        '''
        protein:rna multiple flexible regions, backbone basis, Z-cutoff
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.segbasis = 'backbone'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.zflag = '1'
        self.zcutoff = '-60.0' 
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_zcutoff', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_backbone_zcutoff')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_64(self):
        '''
        protein:rna multiple flexible regions, all basis, no advanced options
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.segbasis = 'all'
        self.seglow = '20'
        self.seghigh = '30'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.flpsegname = 'RNA1'
        self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_all_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'rna_2_all_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_65(self):
        '''
        protein:rna 6 flexible protein regions, heavy basis, low Rg cutoff, high Rg cutoff, directed MC
        '''

        self.trials = '9'
        self.goback = '1'
        self.nsegments = '7'
        self.npsegments = '6'
        self.flpsegname = 'HFQ1,HFQ2,HFQ3,HFQ4,HFQ5,HFQ6'
        self.segbasis = 'heavy'
        self.seglow = '20,20,20,20,20,20'
        self.seghigh = '30,30,30,30,30,30'
        self.lowrg = '46.83'
        self.highrg = '46.89'
        self.directedmc = '46.85'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.psegvariables= [['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_6_1_heavy_lowrg_highrg_directedmc', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_6_1_heavy_lowrg_highrg_directedmc')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_66(self):
        '''
        protein:rna 6 flexible protein regions, 1 flexible rna region, backbone basis, no advanced options
        '''

        self.trials = '9'
        self.nsegments = '7'
        self.npsegments = '7'
        self.flpsegname = 'HFQ1,HFQ2,HFQ3,HFQ4,HFQ5,HFQ6,RNA1'
        self.segbasis = 'backbone'
        self.seglow = '20,20,20,20,20,20,20'
        self.seghigh = '30,30,30,30,30,30,30'
        self.highrg = '400.0'
        self.pdbfile = os.path.join(pdb_data_path,'rna_protein_complex.pdb')
        self.psegvariables= [['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '128', '1', 'rna']]
        self.seed = '1,321'
        gui_mimic_complex_monte_carlo.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.dcdfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'protein_6_1_rna_1_backbone_none', self.dcdfile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'protein_6_1_rna_1_backbone_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)
        if os.path.isfile('a.txt'):
            os.remove('a.txt')
            

if __name__=='__main__':
    main()

