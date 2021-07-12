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
import locale
import sassie.simulate.torsion_angle_md.gui_mimic_torsion_angle_md as gui_mimic_torsion_angle_md
#import gui_mimic_torsion_angle_md as gui_mimic_torsion_angle_md

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'torsion_angle_md') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_torsion_Angle_MD_Filter(MockerTestCase):

    '''
    System integration test for torsion_angle_md_filter.py / sassie 1.0

    Test to see whether torsion_angle_md_filter catches improper input.

    Inputs tested:

	runname:                        string      run name 
    infile:                         string      input pdb or dcd file name
    pdbfile:                        string      input (reference) pdb file name
    outfile:                        string      output dcd file name
    nsteps:                         integer     number of TAMD steps
    topfile:                        string      path and name of topology file
    parmfile:                       string      path and name of parameter file
    keepout:                        integer     keep output files (0==no, 1==yes)
    dcdfreq:                        integer     save individual dcd frequency
    charmexe                        string      path and name of charmm executable file
    temperature                     float       temperature
    rgforce                         float       rg force
    rgvalue                         float       rg value (Angstroms: ignored if rg force = 0.0)
    dna_segnames                    string      names of dsDNA segments
    number_of_flexible_segments     integer     number of flexible segments
    pretamd_min_steps               string      number of pre-TAMD minimization steps

    psegvariables:                              flexible segment variables
                                    string      flexible segment name for each flexible segment
                                    integer     number of flexible regions for each flexible segment
                                    int_array   low residue number for each flexible region
                                    int_array   number of contiguous residues per flexible regions
                                    string      molecule type ('protein', 'rna' or 'dna') for each flexible segment    

    Inputs not tested:

    path                            string      input path (is now part of input file name)
    poll_frequency                  float       time used in time.sleep command (not input by user)


    Use cases tested:

    1.  check if runname has incorrect character
    2.  check reference pdb file
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid   
            ii. PDB file isn't valid
    3.  check if input pdb or dcd file exists
        a)  file doesn't exist
        b)  file exists
            i.  file is not a valid pdb or dcd file
            ii. file is a valid pdb or dcd file
                A. file is not compatible with reference pdb file
                B. file is compatible with reference pdb file
    4.  check if nsteps is > 0
    5.  check if temperature is >= 0
    6.  check if rgvalue >= 0
    7.  check if rgforce >= 0
        a.  rgforce = 0
        b.  rgforce > 0
            i. rgvalue also must be > 0
    8.  check keepout is 0 or 1
    9.  check if pretamd_min_step >= 0
    10. check if number_flexible_segments >= 1    
    11. check dna_segnames
        a.  len(dna_segnames) > 0
            i. check if dna_segname is in PDB file
    12. check if topology file exists
    13. check if parameter file exists
    14. check if CHARMM executable file exists
    15. check if dcdfreq is <= nsteps
    16. check if dcdfreq is a divisor of nsteps
    17. check if pretamd_min_steps is an integer
    18. check flexible segment variables
        a. moltype is 'dna' when len(dna_segnames) > 0
        b. number of flexible segments >= 1 
        c. number of residue ranges = number of flexible regions
        d. moltype for each segment matches the PDB file
        e. number of ranges for each segment is >=1
        f. number of residue ranges matches the number of flexible regions
        g. segment name is in pdb file
        h. low residue can't include n-terminus
        i. low residue + number of contiguous residues can't include c-terminus
        j. pdb file contains the low and high flexible residues listed for each flexible region 
        k. residues in each flexible regions are in pdb file          
        l. residue ranges don't overlap

    '''

    def setUp(self):

        gui_mimic_torsion_angle_md.test_variables(self, paths)


    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_torsion_angle_md.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if pdbmol exists
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        pdbfile=os.path.join(self.path,self.pdbfile)
        expected_error = ['input pdb file, '+pdbfile+', does not exist']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if pdbmol is a valid pdb file
        '''
        self.pdbfile = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        pdbfile=os.path.join(self.path,self.pdbfile)
        expected_error = ['input pdb file, ' +
                          pdbfile + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test if infile exists
        '''
        self.infile = os.path.join(
            module_data_path, 'does_not_exist.dcd')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        infile=os.path.join(self.path,self.infile)
        expected_error = ['file : '+infile+' does not exist']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if infile is a valid pdb file
        '''
        self.infile = os.path.join(
            module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        infile=os.path.join(self.path,self.infile)
        expected_error = ['input trajectory file, '+infile+', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if infile (pdb) is compatible with reference pdb file
        '''
        self.infile = os.path.join(
            module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        infile=os.path.join(self.path,self.infile)
        expected_error = ['input pdb file '+self.pdbfile+' and pdb file '+self.infile+', are not compatible']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if infile (dcd) is compatible with reference pdb file
        '''
        self.infile = os.path.join(
            module_data_path, 'non_matching.dcd')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        infile=os.path.join(self.path,self.infile)
        expected_error = ['input pdb file '+self.pdbfile+' and dcd file '+self.infile+', are not compatible']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if nsteps > 0
        '''
        self.nsteps = '0'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['nsteps = 0?']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if temperature >= 0
        '''
        self.temperature = '-1'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['use a positive temperature, temperature = -1.0']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if rgvalue >= 0
        '''
        self.rgvalue = '-1'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['use a positive Rg value, rg value = -1.0']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if rgforce >= 0
        '''
        self.rgforce = '-1'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['use a positive Rg constraint force, rg force = -1.0']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if rgvalue > 0 when rgforce > 0
        '''
        self.rgforce = '5'
        self.rgvalue = '0'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Rg value must be > 0 when Rg constraint force is > 0, rg value = 0.0, rg force = 5.0']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if keepout = 0 or 1
        '''
        self.keepout = '2'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['option to keep output files needs to be 0 or 1 = 2']
        assert_equals(return_error, expected_error)
                        
    def test_13(self):
        '''
        test if pretamd_min_steps >= 0
        '''
        self.pretamd_min_steps = '-1'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['preliminary minimization steps needs to be >= 0 : -1']
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test if number_flexible_segments > 0
        '''
        self.number_flexible_segments = '0'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of flexible segments needs to be > 0 : 0']
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test if dna_segname is in pdb file
        '''
        self.dna_segnames = 'DNA1'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['dna segname : DNA1 is not in your PDB file']
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if 2nd dna_segname is in pdb file
        '''
        self.pdbfile = os.path.join(pdb_data_path,'c36_dsDNA60_min.pdb')
        self.infile = os.path.join(pdb_data_path,'c36_dsDNA60_min.pdb')
        self.dna_segnames = 'DNA1,DNA3'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['dna segname : DNA3 is not in your PDB file']
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test if topology file exists
        '''
        self.topfile = os.path.join(
            module_data_path, 'does_not_exist.inp')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['topology file does not exist : '+self.topfile]
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if parameter file exists
        '''
        self.parmfile = os.path.join(
            module_data_path, 'does_not_exist.inp')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['parameter file does not exist : '+self.parmfile]
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test if CHARMM executable file exists
        '''
        self.charmmexe = os.path.join(
            module_data_path, 'charmm.exe')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['CHARMM executable file does not exist : '+self.charmmexe]
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        test if dcdfreq <= number of steps
        '''
        self.dcdfreq = '11'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['dcdfreq must less than or equal to the number of TAMD steps: 11']
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test if dcdfreq is a divisor of the number of steps
        '''
        self.dcdfreq = '3'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['dcdfreq must be a divisor of number of TAMD steps: 3']
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test if pretamd_min_steps is an integer
        '''
        self.pretamd_min_steps = '-1.1'
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['preliminary energy minimization steps must be an integer']
        assert_equals(return_error, expected_error)

    def test_23(self):
        '''
        test if number of flexible regions >= 1
        '''
        self.all_snumranges=['0']
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['must have at least one flexible region defined for segment MA']
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test if number of residue ranges matches number of flexible regions
        '''
        self.all_snumranges=['2']
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of residue ranges does not match number of flexible regions']
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test if the molecule type matches the pdb file
        '''
        self.all_moltype=['rna']
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['molecule type rna provided for segment MA does not match the pdb file']
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test if the molecule type matches the pdb file
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.infile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.outfile = 'trunc2a.dcd'
        self.number_flexible_segments = '1'
        self.all_flexible_segnames=['TR2A']
        self.all_snumranges=['1']
        self.all_srlow=['24']       
        self.all_srnum=['6']
        self.all_moltype=['dna']
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            if self.all_moltype[i] == 'dna':
                self.dna_segnames += self.all_flexible_segnames[i] + ','
        if self.dna_segnames and self.dna_segnames[-1] ==',':
            self.dna_segnames = self.dna_segnames[:-1]
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        print 'dna_segnames in test: ', self.dna_segnames
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['molecule type dna provided for segment TR2A does not match the pdb file']
        assert_equals(return_error, expected_error)
                   
    def test_27(self):
        '''
        test if the molecule type matches the pdb file
        '''
        self.pdbfile = os.path.join(pdb_data_path,'ssDNA.pdb')
        self.infile = os.path.join(pdb_data_path,'ssDNA.pdb')
        self.outfile = 'ssDNA.dcd'
        self.number_flexible_segments = '1'
        self.all_flexible_segnames=['DNA1']
        self.all_snumranges=['1']
        self.all_srlow=['11']       
        self.all_srnum=['9']
        self.all_moltype=['rna']
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            if self.all_moltype[i] == 'dna':
                self.dna_segnames += self.all_flexible_segnames[i] + ','
        if self.dna_segnames and self.dna_segnames[-1] ==',':
            self.dna_segnames = self.dna_segnames[:-1]
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        print 'dna_segnames in test: ', self.dna_segnames
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['molecule type rna provided for segment DNA1 does not match the pdb file']
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test if segment name is in pdb file
        '''
        self.all_flexible_segnames=['GAG']
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['found no atoms using filter selection segname[i] == "GAG"', 'segment name not found in pdb file = GAG']
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test if low residue includes the n-terminus
        '''
        self.all_srlow=['1']
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['low residue can not include the n-terminus, reslow = 1']
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test if low residue plus number contiguous exceeds the number of amino acids-1
        '''
        self.all_srnum=['25']
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['your low residue plus number contiguous exceeds the number of amino acids-1 (138), reslow = 114 numcont = 25']
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test if pdb file contains the low residue
        '''
        self.infile = os.path.join(module_data_path,'missing_low_residue.pdb')
        self.pdbfile = os.path.join(module_data_path,'missing_low_residue.pdb')
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input pdb file, '+self.pdbfile+' does not have residue, 114 for segment name MA']
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test if pdb file contains the high residue
        '''
        self.infile = os.path.join(module_data_path,'missing_high_residue.pdb')
        self.pdbfile = os.path.join(module_data_path,'missing_high_residue.pdb')   
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input pdb file, '+self.pdbfile+' does not have residue, 134 for segment name MA']
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test if pdb file contains all residues in the flexible region
        '''
        self.infile = os.path.join(module_data_path,'missing_residue_in_range.pdb')
        self.pdbfile = os.path.join(module_data_path,'missing_residue_in_range.pdb')   
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input pdb file, '+self.pdbfile+' does not have residue, "130" for segment name MA, range = 1 : 139']
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test if residue ranges overlap
        '''
        self.all_snumranges=['2']
        self.all_srlow=['114,120']
        self.all_srnum=['7,14']
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        return_error = gui_mimic_torsion_angle_md.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['residue ranges overlap, low residue (i+1) = 120 high residue (i) = 121']
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
