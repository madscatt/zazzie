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
import sassie.simulate.monte_carlo.monomer.gui_mimic_monomer_monte_carlo as gui_mimic_monomer_monte_carlo
#import gui_mimic_monomer_monte_carlo as gui_mimic_monomer_monte_carlo

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'monomer_monte_carlo') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Monomer_Monte_Carlo_Filter(MockerTestCase):

    '''
    System integration test for generate_filter.py / sassie 1.0

    Test to see whether generate_filter catches improper input.

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


    Use cases tested:

    1.  check if runname has incorrect character
    2.  check input file path permissions 
        a.  no permission error
        b. permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed     
    3.  check pdbfile
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    4.  check if trials is > 0
    5.  check if temperature is >= 0
    6.  check if moltype is "protein" or "rna"
    7.  check if cutoff value is >= 0.001
    8.  check if zflag is 0 or 1
    9.  check if clflag is 0 or 1
    10.  check constraint file
        a.  file doesn't exist
        b.  file exists
    11. check constraint file parameters      
        a.  bad segment name in file
        b.  bad atom name in file
        c.  bad distance value in file
        d.  no distance value in file
        e.  COM or ATM type1 and type 2 in file
        f.  two type definintions in file
        g.  second resid1/resid2 value > first resid1/resid2 value
        h.  first resid1 value is in pdb file
        i.  second resid1 value is in pdb file
        j.  first resid2 value is in pdb file
        k.  second resid2 value is in pdb file
    12. check if pdb file contains the correct moltype
        a.  pdb file contains correct moltype
        b.  pdb file contains more than one moltype
        c.  pdb file contains the wrong moltype
    13. check if low Rg cutoff is higher than high Rg cutoff
    14. check if Rg cutoffs are > 0
        a.  low Rg cutoff is > 0
        b.  high Rg cutoff is > 0
    15. check if number of dtheta values matches the number of ranges
    16. check if number of low residue values matches the number of ranges
    17. check if the number of high residue values matches the number of ranges
    18. check if pdbfile has missing residue (numbers)
    19. check high and low residue values
        a.  pdbfile contains the low and high residues listed for each flexible region
        b.  number of contiguous residues is > 0
        c.  alignment and flexible regions don't overlap
        d.  low residue is lower than the n-terminal amino acid number  #NOT TESTED. This is redundant with lowres test 19a since aa number will not be in pdb file.
        e.  residue range doesn't exceed the number of amino acids      #NOT TESTED. This is redundant with test 19a since aa number will not be in pdb file.  
        f.  residue values increase from low to high
        g.  residue ranges don't overlap with alignment range 
    20. check alignment residue values
        a.  pdb file contains low and high residue amino acids
        b.  alignment range isn't too small (less than 3 points)
    21. check if directed Monte Carlo value is 0 or 1
    22. check overlap basis atoms
        a.  basis cannot be parsed                      #NOT TESTED  didn't find a string that can't be parsed that didn't already trigger another error
        b.  basis can be parsed
            i. check that atom name is in PDB file
            ii.check that atom has VDW paramters        #NOT TESTED  there are no atoms in vdw list that don't have vdw parameters
    '''

    def setUp(self):

        gui_mimic_monomer_monte_carlo.test_variables(self, paths)


    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_monomer_monte_carlo.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)


    def test_2(self):
        '''
        test if path exists
        '''
        self.path = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_monomer_monte_carlo.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path + '  [code = FalseFalseFalse]',
                          'path does not exist']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if directory has read permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_folder')
        ''' see if you can read the directory '''
#        print os.access('empty_folder', os.R_OK)

        self.path = os.path.join('./', 'empty_folder')
        return_error = gui_mimic_monomer_monte_carlo.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.path + '  [code = TrueFalseTrue]', 'read permission not allowed']
        assert_equals(return_error, expected_error)

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')        

    def test_4(self):
        '''
        test if directory has write permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder1')
        ''' see if you can write to the directory '''
        print os.access('empty_folder1', os.W_OK)
        ''' make the directory un-writeable'''
        os.system('chmod a-w empty_folder1')
        ''' see if you can write to the directory '''
        print os.access('empty_folder', os.W_OK)

        self.path = os.path.join('./', 'empty_folder1')
        return_error = gui_mimic_monomer_monte_carlo.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.path + '  [code = TrueTrueFalse]', 'write permission not allowed']
        assert_equals(return_error, expected_error)

        ''' make the directory writeable'''
        os.system('chmod a+w empty_folder1')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder1')        

    def test_5(self):
        '''
        test if pdbfile exists
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'does_not_exist!&@#X.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          self.pdbfile[3:] + ', does not exist']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if pdbfile is a valid pdb file
        '''
        self.pdbfile = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          self.pdbfile[3:] + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if trials is > 0
        '''
        self.trials = '0'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['trials = 0?']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if temperature >=0
        '''
        self.temp = '-1.0'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['use a positive temperature, temperature = -1.0']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if moltype is protein or rna
        '''
        self.moltype = 'rma'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['only protein and rna backbone dihedral move sets are defined, you entered : rma']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test for cutoff value >=0.001
        '''
        self.cutoff = '0.0005'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['use a larger cutoff value, cutoff = 0.0005']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if Z coordinate filter selection is 0 or 1
        '''
        self.zflag = '2'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR in Z coordinate filter selection: zflag == 0 for "no" and 1 for "yes", zflag = 2']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if atomic constraints selection is 0 or 1
        '''
        self.cflag = '2'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR in atomic constraints selection: cflag == 0 for "no" and 1 for "yes", cflag = 2']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if constraint file exists
        '''
        self.cflag = '1'
        self.confile = './does_not_exist.txt'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR in constraint filename selection: file : ./does_not_exist.txt does not exist']
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test for bad seg1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_seg1.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 segment GAG listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test for bad seg2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_seg2.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 segment GAG2 listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test for bad atom1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_atom1.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 atom name XA listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)             

    def test_17(self):
        '''
        test for bad atom2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_atom2.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 atom name ZA listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test for bad distance in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_distance.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 distance value is not appropriate: -100.0"]
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test for no distance in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'no_distance.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 no distance specified or error in line: COM"]
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        test for COM or ATM type1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_type1.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 TYPE1 is not valid (ATM OR COM): CON"]
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test for COM or ATM type2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_type2.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 TYPE2 is not valid (ATM OR COM): ATN"]
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test for two types COM and/or ATM in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'no_type2.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 Two type definitions are required for each constraint (ATM OR COM)"]
        assert_equals(return_error, expected_error)

    def test_23(self):
        '''
        test for second resid1 value equal or less than first
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_resid1.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid values in constraint file for constraint 1 are incorrect: second value is equal or less than first"]
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test for second resid2 value equal or less than first
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_resid2.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid values in constraint file for constraint 0 are incorrect: second value is equal or less than first"]
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test if first value in first resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid1_first.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 0 is not in segment GAG1"]
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test if second value in first resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid1_second.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 433 is not in segment GAG1"]
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test if first value in second resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid2_first.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 0 is not in segment GAG1"]
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test if second value in second resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid2_second.txt')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 435 is not in segment GAG1"]
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test for more than one moltype in pdb file
        '''
        self.pdbfile = os.path.join(module_data_path,'two_moltypes.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR: your PDB structure has more than one molecular type']
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test for wrong moltype in pdb file
        '''
        self.moltype='rna'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR: your PDB structure has been identified as protein but you entered rna']
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test for wrong moltype in pdb file
        '''
        self.moltype='protein'
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR: your PDB structure has been identified as rna but you entered protein']
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test for low Rg cutoff higher than high Rg cutoff
        '''
        self.lowrg = '401.0'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['low Rg cutoff is larger than high Rg cutoff, lowrg = 401.0 highrg = 400.0']
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test if low Rg cutoff > 0
        '''
        self.lowrg = '-1.0'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Rg cutoffs need to be >= zero, lowrg = -1.0 highrg = 400.0']
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test if high Rg cutoff > 0
        '''
        self.lowrg = '-5.0'
        self.highrg = '-1.0'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Rg cutoffs need to be >= zero, lowrg = -5.0 highrg = -1.0']
        assert_equals(return_error, expected_error)

    def test_35(self):
        '''
        test if number of dtheta values matches the number of ranges
        '''

        self.dtheta = '30.0'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of dtheta values does not match the number of ranges, dtheta = [30.0] numranges = 5']
        assert_equals(return_error, expected_error)

    def test_36(self):
        '''
        test if number of lowres values matches the number of ranges
        '''

        self.reslow = '123'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of low residue values does not match the number of ranges, lowres = [123] numranges = 5']
        assert_equals(return_error, expected_error)

    def test_37(self):
        '''
        test if number of highres values matches the number of ranges
        '''

        self.numcont = '22'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of contiguous residue values does not match the number of ranges, contiguous residues = [22] numranges = 5']
        assert_equals(return_error, expected_error)

    def test_38(self):
        '''
        test for missing resid in pdb file
        '''

        self.pdbfile = os.path.join(module_data_path,'missing_resid.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

#       This test doesn't catch missing residue #1 -- only a missing one after the beginning

        ''' check for file error '''
        expected_error = ['amino acid 4 is missing from pdbfile' + self.pdbfile]
        assert_equals(return_error, expected_error)

    def test_39(self):
        '''
        test for low residue number in pdb file
        '''

        self.pdbfile = os.path.join(module_data_path,'first_resid_2.pdb')
        self.reslow = '1,277,354,378,408'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input pdb file, ' + self.pdbfile + ' does not have low residue amino acid, "1" for segment number 0, range = 2 : 431']
        assert_equals(return_error, expected_error)

    def test_40(self):
        '''
        test for high residue number in pdb file
        '''

        self.pdbfile = os.path.join(module_data_path,'first_resid_2.pdb')
        self.reslow = '123,277,354,378,430'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input pdb file, ' + self.pdbfile + ' does not have residue amino acid, "435" for segment number 4, range = 2 : 431']
        assert_equals(return_error, expected_error)

    def test_41(self):
        '''
        test for number of contiguous residues >0
        '''

        self.numcont = '22,6,0,12,5' 
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of contiguous residues "0" should be greater than 0!']
        assert_equals(return_error, expected_error)

    def test_42(self):
        '''
        test if low residue values increase from low to high
        '''

        self.reslow = '123,277,379,378,408'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['low residue values must increase from low to high, reslow = [123, 277, 379, 378, 408]']
        assert_equals(return_error, expected_error)

    def test_43(self):
        '''
        test if residue ranges overlap with alignment range
        '''

        self.lowres1 = '284'
        self.highres1 = '380'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['alignment and flexible ranges should not overlap!']
        assert_equals(return_error, expected_error)

    def test_44(self):
        '''
        test if lowres1 is in pdb file
        '''

        self.lowres1 = '416'
        self.highres1 = '432'        
        self.pdbfile = os.path.join(module_data_path,'first_resid_2.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' + self.pdbfile + ' does not have high residue amino acid, 432, range = 2 : 431']
        assert_equals(return_error, expected_error)

    def test_45(self):
        '''
        test if highres1 is in pdb file
        '''

        self.lowres1 = '1'
        self.highres1 = '15'
        self.pdbfile = os.path.join(module_data_path,'first_resid_2.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' + self.pdbfile + ' does not have low residue amino acid, 1, range = 2 : 431']
        assert_equals(return_error, expected_error)

    def test_46(self):
        '''
        test if alignment basis is less than 3 points
        '''

        self.lowres1 = '416'
        self.highres1 = '418'        
        self.pdbfile = os.path.join(module_data_path,'first_resid_2.pdb')
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['alignment basis is too small (less than 3 points) or low residue > high residue']
        assert_equals(return_error, expected_error)

    def test_47(self):
        '''
        test if directedmc >=0
        '''

        self.directedmc = '-1.0'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['directed Monte Carlo needs to be 0 or a float > 0 (the "goal Rg") ... you entered: -1.0']
        assert_equals(return_error, expected_error)

    def test_45(self):
        '''
        test if basis atom is in PDB file 
        '''

        self.basis = 'CQ'
        return_error = gui_mimic_monomer_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['overlap basis atom CQ is not in your PDB file']
        assert_equals(return_error, expected_error)


               

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
