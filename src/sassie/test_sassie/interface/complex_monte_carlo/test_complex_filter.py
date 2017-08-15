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
import sassie.simulate.complex_monte_carlo.gui_mimic_complex_monte_carlo as gui_mimic_complex_monte_carlo
#import gui_mimic_complex_monte_carlo as gui_mimic_complex_monte_carlo

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'complex_monte_carlo') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Complex_Monte_Carlo_Filter(MockerTestCase):

    '''
    System integration test for complex_filter.py / sassie 1.0

    Test to see whether complex_filter catches improper input.

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

    psegvariables:  
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
    5.  check if goback is > 0
    6.  check if temperature is >= 0
    7.  check if zflag is 0 or 1        #NOTE:  zcutoff test is commented out in complex_filter.py
    8.  check if clflag is 0 or 1
    9.  check constraint file
        a.  file doesn't exist
        b.  file exists
    10. check constraint file parameters      
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
    11. check if directed Monte Carlo value is 0 or 1        
    12. check if low Rg cutoff is higher than high Rg cutoff
    13. check if Rg cutoffs are > 0
        a.  low Rg cutoff is > 0
        b.  high Rg cutoff is > 0
    14. check if number of segments is >= 1
    15. check if number of flexible segments is >= 1
    16. check if the number of flexible segments is <= the number of segments
    17. check that number of basis names is the same as the number of segments (for basis != 'all', 'backbone' or 'heavy')
    18. check if number of flexible segment names matches the number of flexible segments
    19. check if number of alignment low residues matches the number of flexible segments
    20. check if the number of alignment high residues matches the number of flexible segments
    21. check if each (flexible?) segment in the pdb file contains the correct moltype  #NOT TESTED  need to loop over segment names
    22. check overlap basis atoms
        a. check that atom name is in PDB file      #NOT TESTED  error handling is commented out in complex_filter.py
        b. check that atom has VDW paramters        #NOT TESTED  there are no atoms in vdw list that don't have vdw parameters
    23. check overlap in initial structure          #NOT TESTED  generates a warning only; program does not exit due to overlap in initial structure
    24. check if flexible segments are in PDB file
    25. check if total number of segments matches the number of segments in PDB file
    26. check if pdbfile has missing residue (numbers)  #NOT TESTED need to loop over segment names
    27. check flexible segment variables
        a. angle values are float types
        b. angle values are in the range 0.0 to 180.0
        c. number of ranges for each segment is an integer
        d. number of ranges for each segment is >=1
        e. low resid is an integer array
        f. number of contiguous residues is an integer array
        g. moltype for each segment matches the PDB file
        h. PDB file contains low and high alignment residues listed for each flexible region
        i. low alignment resid < high alignment resid               
        j. alignment range for each segment isn't too small (less than 3 points)   
        k. number of angle values matches the number of ranges
        l. number of low residue values matches the number of ranges
        m. number of contiguous residues matches the number of ranges
        n. PDB file contains the low and high flexible residues listed for each flexible region 
        o. number of contiguous residues is >= 0
        p. alignment and flexible regions don't overlap
        q. low residue can't include n-terminus (for numranges > 1)
        r. low residue values increase from low to high (for numranges > 1)
        s. residue ranges don't overlap (for numranges > 1)
        t. low residue + number of contiguous doesn't exceed number of amino acids-1 (for numranges > 1)
        u. low residue + number of contiguous doesn't exceed number of amino acids-1 (for numranges = 1)

    '''

    def setUp(self):

        gui_mimic_complex_monte_carlo.test_variables(self, paths)


    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_complex_monte_carlo.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)


    def test_2(self):
        '''
        test if path exists
        '''
        self.path = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_complex_monte_carlo.run_module(
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
        return_error = gui_mimic_complex_monte_carlo.run_module(
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
        self.path = os.path.join(module_data_path, 'no_write_permission')
        return_error = gui_mimic_complex_monte_carlo.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path +
                          '  [code = TrueTrueFalse]', 'write permission not allowed']
        assert_equals(return_error, expected_error)


    def test_5(self):
        '''
        test if pdbfile exists
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'does_not_exist!&@#X.pdb')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          self.pdbfile + ', does not exist']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if pdbfile is a valid pdb file
        '''
        self.pdbfile = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          self.pdbfile + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if trials is > 0
        '''
        self.trials = '0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['trials = 0?']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if goback is > 0
        '''
        self.goback = '0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['goback = 0?']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if temperature >=0
        '''
        self.temp = '-1.0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['use a positive temperature, temperature = -1.0']
        assert_equals(return_error, expected_error)               

    def test_10(self):
        '''
        test if Z coordinate filter selection is 0 or 1
        '''
        self.zflag = '2'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['zflag == 0 for "no" and 1 for "yes", zflag = 2']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if atomic constraints selection is 0 or 1
        '''
        self.cflag = '2'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['cflag == 0 for "no" and 1 for "yes", cflag = 2']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if constraint file exists
        '''
        self.cflag = '1'
        self.confile = './does_not_exist.txt'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['file : ./does_not_exist.txt does not exist']]
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test for bad seg1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_seg1.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 segment PAI listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test for bad seg2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_seg2.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 segment VN2 listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test for bad atom1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_atom1.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 atom name XA listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)             

    def test_16(self):
        '''
        test for bad atom2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_atom2.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 atom name ZA listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test for bad distance in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_distance.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 distance value is not appropriate: -100.0"]
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test for no distance in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'no_distance.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 no distance specified or error in line: COM"]
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test for COM or ATM type1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_type1.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 TYPE1 is not valid (ATM OR COM): CON"]
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        test for COM or ATM type2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_type2.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 TYPE2 is not valid (ATM OR COM): ATN"]
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test for two types COM and/or ATM in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'no_type2.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 Two type definitions are required for each constraint (ATM OR COM)"]
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test for second resid1 value equal or less than first
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_resid1.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid values in constraint file for constraint 1 are incorrect: second value is equal or less than first"]
        assert_equals(return_error, expected_error)

    def test_23(self):
        '''
        test for second resid2 value equal or less than first
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_resid2.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid values in constraint file for constraint 0 are incorrect: second value is equal or less than first"]
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test if first value in first resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid1_first.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 0 is not in segment PAI1"]
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test if second value in first resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid1_second.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 433 is not in segment PAI1"]
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test if first value in second resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid2_first.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 0 is not in segment VN1"]
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test if second value in second resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid2_second.txt')
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 135 is not in segment VN1"]
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test for low Rg cutoff higher than high Rg cutoff
        '''
        self.lowrg = '51.0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['low Rg cutoff is larger than high Rg cutoff, lowrg = 51.0 highrg = 50.0']
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test if low Rg cutoff > 0
        '''
        self.lowrg = '-1.0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Rg cutoffs need to be >= zero, lowrg = -1.0 highrg = 50.0']
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test if high Rg cutoff > 0
        '''
        self.lowrg = '-5.0'
        self.highrg = '-1.0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Rg cutoffs need to be >= zero, lowrg = -5.0 highrg = -1.0']
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test if number of segments >= 1
        '''
        self.nsegments = '0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of segments 0 should be equal/greater than 1!']
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test if number of flexible segments >= 1
        '''
        self.npsegments = '0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of flexible segments 0 should be equal/greater than 1!']
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test if the number of flexible segments <= the number of segments
        '''
        self.npsegments = '3'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of flexible segments 3 should be equal/less than the number of total segments 2!']
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test if the number of segment basis matches the number of segments (for basis != 'all', 'backbone' or 'heavy')
        '''
        self.segbasis = 'CA'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of segment basis does not match the number of segments: number of segbasis = 1 number of segments = 2',
   'segment overlap basis entries can be "heavy", "backbone", "all", or a comma delimited list of atom names ... one for each segment']
        assert_equals(return_error, expected_error)

    def test_35(self):
        '''
        test if the number of flexible segment names matches the number of flexible segments
        '''
        self.flpsegname = 'VN1,VN2'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of flexible segment names does not match the number of flexible segments: number of flexible segment names = 2 number of flexible segments = 1']
        assert_equals(return_error, expected_error)

    def test_36(self):
        '''
        test if the number of alignment low residues matches the number of flexible segments
        '''
        self.seglow = '1,1'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of alignment low residues does not match the number of flexible segments: number of alignment low residues = 2 number of  flexible segments = 1']
        assert_equals(return_error, expected_error)

    def test_37(self):
        '''
        test if the number of alignment high residues matches the number of flexible segments
        '''
        self.seghigh = '30,30'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of alignment high residues does not match the number of flexible segments: number of alignment high residues = 2 number of  flexible segments = 1']
        assert_equals(return_error, expected_error)

    def test_38(self):
        '''
        test if single flexible segment is in PDB file 
        '''
        self.flpsegname = 'VN2'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The flexible segment name "VN2" is not found in the pdb file!']
        assert_equals(return_error, expected_error)

    def test_39(self):
        '''
        test if flexible segments are in PDB file (first segment name is in file; second segment name isn't in file)
        '''
        self.npsegments = '2'
        self.seglow = '1,1'
        self.seghigh = '30,30'
        self.flpsegname = 'VN1,VN2'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The flexible segment name "VN2" is not found in the pdb file!']
        assert_equals(return_error, expected_error)

    def test_40(self):
        '''
        test if total number of segments is equal to number of segments in PDB file
        '''
        self.nsegments = '3'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the total number of segments entered does NOT match the number of segments in the pdb file']
        assert_equals(return_error, expected_error)

    def test_41(self):
        '''
        test if directedmc >=0
        '''
        self.directedmc = '-1.0'
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['directed Monte Carlo needs to be 0 or a float > 0 (the "goal Rg") ... you entered: -1.0']
        assert_equals(return_error, expected_error)

    def test_42(self):
        '''
        test if number of ranges is an integer type
        '''
        self.psegvariables= [['1.0', '30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The number of ranges "1.0" for flexible segment number 1 in the flexible segment input fields should be an integer type!']]
        assert_equals(return_error, expected_error)

    def test_43(self):
        '''
        test if number of ranges is >= 1
        '''
        self.psegvariables= [['0', '30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The number of ranges "0" for flexible segment number 1 in the flexible segment input fields should be equal/greater than 1!']]
        assert_equals(return_error, expected_error)

    def test_44(self):
        '''
        test if the angle value is a float type
        '''
        self.psegvariables= [['1', '3o', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The angle value "3o" should be a float type!']]
        assert_equals(return_error, expected_error)

    def test_45(self):
        '''
        test if the angle value is between 0 and 180
        '''
        self.psegvariables= [['1', '190', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The angle value "190" should be in the range of (0.0,180.0)!']]
        assert_equals(return_error, expected_error)

    def test_46(self):
        '''
        test if the low resid is an integer array
        '''
        self.psegvariables= [['1', '30', '40.0', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The low resid "40.0" for flexible segment number 1 in the flexible segment input fields should be an integer array!']]
        assert_equals(return_error, expected_error)

    def test_47(self):
        '''
        test if the number of contiguous residues is an integer array
        '''
        self.psegvariables= [['1', '30', '40', '89.0', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The number of contiguous residues "89.0" for flexible segment number 1 in the flexible segment input fields should be an integer array!']]
        assert_equals(return_error, expected_error)               

    def test_48(self):
        '''
        test if the molecule type provided for the flexible segment matches that in the PDB file
        '''
        self.psegvariables= [['1', '30', '40', '89', 'rna']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The molecule type "rna" provided for flexible segment number 1 in the flexible segment input fields does not match the pdb file!']]
        assert_equals(return_error, expected_error)

    def test_49(self):
        '''
        test if the flexible residue in the input PDB file has low alignment residue
        '''
        self.pdbfile = os.path.join(module_data_path,'missing_resid.pdb')
        self.psegvariables= [['1', '30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['input pdb file, ' + str(self.pdbfile) + ' does not have low alignment amino acid residue, 1, range = 2 : 130']]
        assert_equals(return_error, expected_error)

    def test_50(self):
        '''
        test if the flexible residue in the input PDB file has high alignment residue
        '''
        self.pdbfile = os.path.join(module_data_path,'missing_resid1.pdb')
        self.psegvariables= [['1', '30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['input pdb file, ' + str(self.pdbfile) + ' does not have high alignment amino acid residue, 39, range = 1 : 130']]
        assert_equals(return_error, expected_error)

    def test_51(self):
        '''
        test if low alignment residue < high alignment residue
        '''
        self.seglow = '20'
        self.seghigh = '1'
        self.psegvariables= [['1', '30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['alignment basis is too small (less than 3 points) or low residue > high residue']]
        assert_equals(return_error, expected_error)

    def test_52(self):
        '''
        test if alignment range > 3
        '''
        self.seglow = '1'
        self.seghigh = '3'
        self.psegvariables= [['1', '30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['alignment basis is too small (less than 3 points) or low residue > high residue']]
        assert_equals(return_error, expected_error)

    def test_53(self):
        '''
        test number of angle values matches the number of ranges
        '''
        self.psegvariables= [['1', '30,30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['the number of dtheta values does not match the number of ranges, dtheta = [30.0, 30.0] numranges = 1']]
        assert_equals(return_error, expected_error)

    def test_54(self):
        '''
        test number of low residue values matches the number of ranges
        '''
        self.psegvariables= [['1', '30', '40,130', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['the number of low residue values does not match the number of ranges, lowres = [40, 130] numranges = 1']]
        assert_equals(return_error, expected_error)

    def test_55(self):
        '''
        test number of contiguous residues matches the number of ranges
        '''
        self.psegvariables= [['1', '30', '40', '89,2', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['the number of contiguous residues does not match the number of ranges, contiguous = [89, 2] numranges = 1']]
        assert_equals(return_error, expected_error)

    def test_56(self):
        '''
        test if low flexible residue is in PDB file
        '''
        self.psegvariables= [['1', '30', '131', '2', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['Input pdb file, ' + str(self.pdbfile) + ' does not have low residue amino acid, "131" for segment number 1, range = 1 : 130']]
        assert_equals(return_error, expected_error)

    def test_57(self):
        '''
        test if low+contiguous flexible residue is in PDB file
        '''
        self.pdbfile = os.path.join(module_data_path,'missing_resid2.pdb')
        self.psegvariables= [['1', '30', '40', '89', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['Input pdb file, ' + str(self.pdbfile) + ' does not have low+contiguous residue amino acid, "129" for segment number 1, range = 1 : 130']]
        assert_equals(return_error, expected_error)

    def test_58(self):
        '''
        test if number of contiguous residues is >=0
        '''
        self.psegvariables= [['1', '30', '40', '-1', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The number of contiguous residues "-1" should be greater than 0!']]
        assert_equals(return_error, expected_error)

    def test_59(self):
        '''
        test if alignment and flexible ranges overlap
        '''
        self.psegvariables= [['1', '30', '39', '90', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['alignment and flexible ranges should not overlap!']]
        assert_equals(return_error, expected_error)

    def test_60(self):
        '''
        test if low residue includes the n-terminus
        '''
        self.seglow = '80'
        self.seghigh = '90'
        self.psegvariables= [['2', '30,30', '1,10', '10,10', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['low residue can not include the n-terminus, reslow = [1, 10]']]
        assert_equals(return_error, expected_error)

    def test_61(self):
        '''
        test if low residue values increase from low to high
        '''
        self.seglow = '80'
        self.seghigh = '90'
        self.psegvariables= [['2', '30,30', '15,10', '10,10', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['low residue values must increase from low to high, reslow = 15']]
        assert_equals(return_error, expected_error)

    def test_62(self):
        '''
        test if residue ranges overlap
        '''
        self.seglow = '80'
        self.seghigh = '90'
        self.psegvariables= [['2', '30,30', '2,10', '10,10', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['low residue values plus number contiguous overlap, reslow = 2 numcont = 10']]
        assert_equals(return_error, expected_error)

    def test_63(self):
        '''
        test if low residue plus number contiguous exceeds the number of amino acids-1 (numranges > 1)
        '''

        self.psegvariables= [['2', '30,30', '40,60', '18,70', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['your low residue plus number contiguous exceeds the number of amino acids-1 (129), reslow = 60 numcont = 70']]
        assert_equals(return_error, expected_error)

    def test_64(self):
        '''
        test if low residue plus number contiguous exceeds the number of amino acids-1 (numranges = 1)
        '''

        self.psegvariables= [['1', '30', '40', '90', 'protein']]
        return_error = gui_mimic_complex_monte_carlo.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [['your low residue plus number contiguous exceeds the number of amino acids-1 (129), reslow = 40 numcont = 90']]
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
