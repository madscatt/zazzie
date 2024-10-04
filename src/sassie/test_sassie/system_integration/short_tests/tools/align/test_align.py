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
import sassie.tools.gui_mimic_align as gui_mimic_align
#import gui_mimic_align as gui_mimic_align

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'tools', 'align') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Align(MockerTestCase):

    '''
    System integration test for align.py / sassie 1.0

    align.py overlaps molecules from a dcd/pdb file onto another molecule over a given basis. 

    Use cases:

    1.  Molecule 2 input file
        a.  input file is a PDB file
        b.  input file is a DCD file
    2.  Basis string (one for each molecule)
        a.  basis string is 'CA'
        b.  basis string is a different atom. There is currently no dictionary of acceptable basis atoms?
    3.  Extra basis statement ('None' or appropriate command.  Two examples of appropriate commands are given in the align docs.)
        a.  extra basis statement is 'None' 
        b.  extra basis statement is not 'None'
    4.  Z coordinate filter
        a.  zflag is 'True' (provide z cutoff value)
        b.  zflag is 'False' (z cutoff value is ignored)


    Inputs tested:

    runname:      string      project name                          
    path:         string      input/output filepath                 
    pdbmol1:      string      reference pdb (mol 1)                 : 2 reference pdb files (1 protein, 1 RNA)
    pdbmol2:      string      input pdb file (mol 2)                : 2 pdb files (1 protein, 1 RNA)
    infile:       string      input (pdb or dcd) filename (mol 2)   : 1 pdb file, 1 dcd file for each above 
    basis1:       string      basis for molecule 1                  : 2 basis atoms (CA for protein, P for RNA)
    basis2:       string      basis for molecule 2                  : 2 basis atoms (CA for protine, P for RNA)
    lowres1:      integer     low residue for overlap molecule 1    : 2 values (protein, RNA)
    highres1:     integer     high residue for overlap molecule 1   : 2 values (protein, RNA)
    lowres2:      integer     low residue for overlap molecule 2    : 2 values (protein, RNA)
    highres2:     integer     high residue for overlap molecule 2   : 2 values (protein, RNA)

    zflag:        boolean     flag for zcutoff value                : 2 values ('True' or 'False')
    zcutoff:      float       cutoff value for z coordinate         : any float

    Conditions:

    basis 1 must equal basis 2
    overlap regions for molecules 1 and 2 must have the same number of atoms
    zcutoff == True: frames containing ANY atoms with a z-value less than the cutoff are not written

    Inputs not tested (deprecated functionality):

    ebasis1:      string      extra basis statement molecule 1 
    ebasis2:      string      extra basis statement molecule 2 

    Test tree:


                                                *********************
                                                *   project name    *
                                                * input/output path *
                                                *********************
                                                *                   *
                                              *                       *
                                            *                           *
                            ********************                   *********************
                            *    mol 1: PDB    *                   *     mol 1: PDB    *
                            *    mol 2: PDB    *                   *     mol 2: DCD    *
                            ********************                   *********************
                                *            *                         *             *
                             *                  *                   *                   *
             **********************  *********************   *********************  *********************
             *       Protein      *  *       RNA         *   *       Protein     *  *       RNA         *
             *    basis 1: CA     *  *   basis 1: P      *   *    basis 1: CA    *  *   basis 1: P      *
             *    basis 2: CA     *  *   basis 2: P      *   *    basis 2: CA    *  *   basis 2: P      *
             *  lowres 1=lowres 2 *  * lowres 1=lowres2  *   * lowres 1=lowres 2 *  * lowres 1=lowres2  *
             * highres 1=highres 2*  *highres 1=highres 2*   *highres 1=highres 2*  *highres 1=highres 2*
             *                    *  *                   *   *  zcutoff tested   *  *                   *
             **********************  *********************   *********************  *********************



    '''

    module = 'align'

    def setUp(self):

       gui_mimic_align.test_variables(self, paths)


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

    def test_1(self):
        '''
        test dcd input and CA,CA basis
        '''

        gui_mimic_align.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbmol2)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'aligned_hiv1_gag_20_frames.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbmol2)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm output minmax file is correct '''
        outfile = os.path.join(self.runname, self.module,
                               self.ofile + '.minmax')
        correct_outfile = os.path.join(
            module_data_path, 'aligned_hiv1_gag_20_frames.dcd.minmax')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_2(self):
        '''
        test pdb input and CA,CA basis
        '''
        self.infile = os.path.join(pdb_data_path, 'hiv1_gag.pdb')
        self.ofile = 'aligned_hiv1_gag.pdb'

        gui_mimic_align.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'aligned_hiv1_gag.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm output minmax file is correct '''
        outfile = os.path.join(self.runname, self.module,
                               self.ofile + '.minmax')
        correct_outfile = os.path.join(
            module_data_path, 'aligned_hiv1_gag.pdb.minmax')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_3(self):
        '''
        test dcd input and P,P basis
        '''
        self.pdbmol1 = os.path.join(pdb_data_path, 'rna.pdb')
        self.pdbmol2 = os.path.join(pdb_data_path, 'rna.pdb')
        self.infile = os.path.join(dcd_data_path, 'rna_100_frames.dcd')
        self.ofile = 'aligned_rna_100_frames.dcd'
        self.basis1 = 'P'
        self.basis2 = 'P'
        self.lowres1 = '20'
        self.lowres2 = '20'
        self.highres1 = '60'
        self.highres2 = '60'

        gui_mimic_align.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbmol2)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'aligned_rna_100_frames.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbmol2)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm output minmax file is correct '''
        outfile = os.path.join(self.runname, self.module,
                               self.ofile + '.minmax')
        correct_outfile = os.path.join(
            module_data_path, 'aligned_rna_100_frames.dcd.minmax')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_4(self):
        '''
        test pdb input and P,P basis
        '''
        self.pdbmol1 = os.path.join(pdb_data_path, 'rna.pdb')
        self.pdbmol2 = os.path.join(pdb_data_path, 'rna.pdb')
        self.infile = os.path.join(pdb_data_path, 'rna.pdb')
        self.ofile = 'aligned_rna.pdb'
        self.basis1 = 'P'
        self.basis2 = 'P'
        self.lowres1 = '20'
        self.lowres2 = '20'
        self.highres1 = '60'
        self.highres2 = '60'

        gui_mimic_align.run_module(self)

        ''' confirm output pdb file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        correct_outfile = os.path.join(module_data_path, 'aligned_rna.pdb')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

        ''' confirm output minmax file is correct '''
        outfile = os.path.join(self.runname, self.module,
                               self.ofile + '.minmax')
        correct_outfile = os.path.join(
            module_data_path, 'aligned_rna.pdb.minmax')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)

    def test_5(self):
        '''
        test dcd input, CA,CA basis and zcutoff
        '''

        self.ofile = 'aligned_hiv1_gag_20_frames_zcutoff.dcd'
        self.zflag = 'True'
        self.zcutoff = '-66.0'
        
        gui_mimic_align.run_module(self)

        ''' confirm output dcd file is correct '''
        outfile = os.path.join(self.runname, self.module, self.ofile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbmol2)
        molecule.read_dcd(outfile)
        result_coor = molecule.coor()

        correct_outfile = os.path.join(
            module_data_path, 'aligned_hiv1_gag_20_frames_zcutoff.dcd')
        correct_molecule = sasmol.SasMol(0)
        correct_molecule.read_pdb(self.pdbmol2)
        correct_molecule.read_dcd(correct_outfile)
        correct_coor = correct_molecule.coor()

        self.assert_list_almost_equal(
            correct_coor, result_coor, self.precision)

        ''' confirm output minmax file is correct '''
        outfile = os.path.join(self.runname, self.module,
                               self.ofile + '.minmax')
        correct_outfile = os.path.join(
            module_data_path, 'aligned_hiv1_gag_20_frames_zcutoff.dcd.minmax')
        assert_equals(filecmp.cmp(outfile, correct_outfile), True)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

