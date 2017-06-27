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
import sassie.simulate.two_body_grid.gui_mimic_two_body_grid as gui_mimic_two_body_grid
#import gui_mimic_two_body_grid as gui_mimic_two_body_grid

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'simulate', 'two_body_grid') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Two_Body_Grid(MockerTestCase):

    '''
		TWO-BODY GRID is the program to move a molecule on a grid.
	
		Molecule 1 is the reference molecule.
		Molecule 2 is the molecule to be moved on the grid 
					The molecules must have "residue" fields


    INPUT:	
        run name
		PDB file for molecule 1 (reference structure)
        PDB file for molecule 2 (structure to be moved on the grid)
        output DCD file name
        accept supplied initial position for molecule 2 (0=no; use COM of molecule 2, 1=yes)
        initial position of molecule 2 (x,y,z)
        number of x,y,z moves (nx,ny,nz)
        dx,dy,dz step sizes (Angstroms)
        number of angular moves (nthetax,nthetay,nthetaz)
        dtheta (dthetax,dthetay,dthetaz) step sizes (degress)

        Advanced input:

        overlap basis atom (only CA is accepted)
        overlap cutoff distance
        low Rg cutoff
        high Rg cutoff
        enable zcutoff flag (0=no, 1=yes)
        zcutoff value (discard structures with any z-axis coordinates less than this value)
        enable atomic constraint flag (0=no, 1=yes)
        name of file describing additional constraints to check before accepting a structure

        Excluded residue input:

        number of segments in molecule 1 with residues excluded from overlap check        
        name of segments in molecule 1 with residues excluded from overlap check
        first amino acid residue per segment (molecule 1) 
        number of contiguous amino acid residues per segment (molecule 1)
        number of segments in molecule 2 with residues excluded from overlap check        
        name of segments in molecule 2 with residues excluded from overlap check
        first amino acid residue per segment (molecule 2) 
        number of contiguous amino acid residues per segment (molecule 2)

    OUTPUT:
        original PDB files for molecule 1 and molecule 2
        a DCD file with aligned coordinates
        a reference PDB file for the DCD file
        file containing x,y,z,thetax,thetay,thetaz coordinates of molecule 2

 
    Use cases:

    1.  initial position
        a.  use COM of molecule 2
        b.  use supplied initial position
    2.  advanced input options
        a.  low Rg cutoff
        b.  high Rg cutoff
        c.  Z-cutoff
        d.  contraints
    3.  exclude residues from overlap check
        a.  no residues excluded
        b.  exclude 


    Inputs tested:

    runname:        string      project name                          
    path:           string      input file path                 
    ofile  :        string      name of output dcd file containing accepted structures       
    pdbmol1:        string      name of pdb file containing molecule 1 reference structure
    pdbmol2:        string      name of pdb file containing molecule 2 structure to move
    accpos:         integer     accept supplied initial position for molecule 2 (0=no; use COM of molecule 2, 1=yes)
    position:       int_array   initial position of molecule 2
    trans:          int_arry    number of x,y,z, moves
    dtrans:         float_array dx,dy,dz step sizes (Angstroms)
    theta:          int_array   number of angular moves
    dtheta:         float_array dx,dy,dz step sizes (degrees)
    basis:          string      overlap basis atom (only CA is accepted)
    cutoff:         float       overlap cutoff distance 
    lowrg           float       low Rg cutoff value                 
    highrg          float       high Rg cutoff value                
    zflag           integer     enable zcutoff flag (0=no, 1=yes)
    zcutoff         float       zcutoff value (discard structures with any z-axis coordinates less than this value)
    cflag           integer     enable atomic constraint flag (0=no, 1=yes)
    confile         string      name of file describing additional constraints to check before accepting a structure
    nexsegments1    integer     number of segments in molecule 1 with residues excluded from overlap check
    nsegments1      string      name of segments in molecule 1 with residues excluded from overlap check
    reslow1         int_array   first amino acid residue per segment (molecule 1) 
    numcont1        int_array   number of contiguous amino acid residues per segment (molecule 1)
    nexsegments2    integer     number of segments in molecule 2 with residues excluded from overlap check
    nsegments2      string      name of segments in molecule 2 with residues excluded from overlap check
    reslow2         int_array   first amino acid residue per segment (molecule 2) 
    numcont2        int_array   number of contiguous amino acid residues per segment (molecule 2)


    Test tree:

    project name
    input PDB file
    output DCD file

*******************************************
    COM inital position for molecule 2
*******************************************

   no adv options    low Rg cutoff       high Rg cutoff     constraints         Z-cutoff            no adv options     no adv options     no adv options
   no excl residues  no excl residues    no excl residues   no excl residues    no excl residues    seg 1 excl res     seg 2 excl res     segs 1&2 excl res


 *******************************************
    supplied initial position for molecule 2
********************************************

   no adv options    low Rg cutoff       high Rg cutoff     constraints         Z-cutoff            no adv options     no adv options     no adv options
   no excl residues  no excl residues    no excl residues   no excl residues    no excl residues    seg 1 excl res     seg 2 excl res     segs 1&2 excl res

    '''

    module = 'two_body_grid'

    def setUp(self):

       gui_mimic_two_body_grid.test_variables(self, paths)


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
        ignoring files that have date stamps: self.ofile, *.sassie_json
        and *.sassie_log
        '''
        dirs_cmp = filecmp.dircmp(dir1, dir2,[self.ofile])    
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
        COM molecule 2 position, no advanced options, no excluded residues
        '''
  
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_2(self):
        '''
        COM molecule 2 position, low Rg cutoff, no excluded residues
        '''
  
        self.lowrg = '50.0'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_lowrg_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_lowrg_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_3(self):
        '''
        COM molecule 2 position, high Rg cutoff, no excluded residues
        '''
  
        self.highrg = '50.0'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_highrg_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_highrg_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_4(self):
        '''
        COM molecule 2 position, constraints, no excluded residues
        '''
  
        self.cfile = '1'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_constraints_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_constraints_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)


    def test_5(self):
        '''
        COM molecule 2 position, zcutoff, no excluded residues
        '''
  
        self.zflag = '1'
        self.zcutoff = '-90.0'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_zcutoff_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_zcutoff_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_6(self):
        '''
        COM molecule 2 position, no advanced options, seg 1 excluded residues
        '''
  
        self.nexsegments1 = '3'
        self.nsegments1 = 'INT1, INT2, LED1'
        self.reslow1 = '271, 271, 94'
        self.numcont1 = '18, 18, 31'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_seg1', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_seg1')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_7(self):
        '''
        COM molecule 2 position, no advanced options, seg 2 excluded residues
        '''
  
        self.nexsegments2 = '3'
        self.nsegments2 = 'INT3, INT4, LED2'
        self.reslow2 = '271, 271, 94'
        self.numcont2 = '18, 18, 31'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_seg2', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_seg2')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_8(self):
        '''
        COM molecule 2 position, no advanced options, seg1 & seg 2 excluded residues
        '''
  
        self.nexsegments1 = '3'
        self.nsegments1 = 'INT1, INT2, LED1'
        self.reslow1 = '271, 271, 94'
        self.numcont1 = '18, 18, 31'
        self.nexsegments2 = '3'
        self.nsegments2 = 'INT3, INT4, LED2'
        self.reslow2 = '271, 271, 94'
        self.numcont2 = '18, 18, 31'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_seg1seg2', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'com_none_seg1seg2')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)


    def test_9(self):
        '''
        supplied molecule 2 position, no advanced options, no excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_10(self):
        '''
        supplied molecule 2 position, lowrg, no excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        self.lowrg = '50.0'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_lowrg_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_lowrg_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_11(self):
        '''
        supplied molecule 2 position, highrg, no excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        self.highrg = '47.0'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_highrg_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_highrg_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
        
    def test_12(self):
        '''
        supplied molecule 2 position, constraints, no excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        self.cfile = '1'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_constraints_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_constraints_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_13(self):
        '''
        supplied molecule 2 position, zcutoff, no excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        self.zflag = '1'
        self.zcutoff = '-80.0'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_zcutoff_none', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_zcutoff_none')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_14(self):
        '''
        supplied molecule 2 position, no advanced options, seg 1 excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        self.nexsegments1 = '3'
        self.nsegments1 = 'INT1, INT2, LED1'
        self.reslow1 = '271, 271, 94'
        self.numcont1 = '18, 18, 31'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_seg1', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_seg1')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_15(self):
        '''
        supplied molecule 2 position, no advanced options, seg 2 excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        self.nexsegments2 = '3'
        self.nsegments2 = 'INT3, INT4, LED2'
        self.reslow2 = '271, 271, 94'
        self.numcont2 = '18, 18, 31'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_seg2', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_seg2')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_16(self):
        '''
        supplied molecule 2 position, no advanced options, seg1 & seg 2 excluded residues
        '''
  
        self.accpos='1'
        self.pos='-20, -20, -20'
        self.nexsegments1 = '3'
        self.nsegments1 = 'INT1, INT2, LED1'
        self.reslow1 = '271, 271, 94'
        self.numcont1 = '18, 18, 31'
        self.nexsegments2 = '3'
        self.nsegments2 = 'INT3, INT4, LED2'
        self.reslow2 = '271, 271, 94'
        self.numcont2 = '18, 18, 31'
        gui_mimic_two_body_grid.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        pdbfile = os.path.join(pdb_data_path,'iase1iase2.pdb')
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(pdbfile)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_seg1seg2', self.ofile)
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(pdbfile)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm correct files are in output directory '''

        outdirectory = os.path.join(self.runname, self.module)
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'supplied_none_seg1seg2')
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

