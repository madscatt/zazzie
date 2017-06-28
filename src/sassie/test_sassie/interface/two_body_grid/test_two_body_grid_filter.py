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
import sassie.simulate.two_body_grid.gui_mimic_two_body_grid as gui_mimic_two_body_grid
#import gui_mimic_two_body_grid as gui_mimic_two_body_grid

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'two_body_grid') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Two_Body_Grid_Filter(MockerTestCase):

    '''
    System integration test for two_body_grid_filter.py / sassie 1.0

    Test to see whether two_body_grid_filter catches improper input.

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


    Use cases tested:

    1.  check if runname has incorrect character
    2.  check input file path permissions 
        a.  no permission error
        b. permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed     
    3.  check pdbfile1
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    4.  check pdbfile2
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid            
    5.  check if accpos is 0 or 1
    6.  check if pos consists of three values
    7.  check if trans consists of three values
    8.  check if trans consists of three values > 0
    9.  check if dtrans consists of three values
    10. check if theta consists of three values
    11. check if dtheta consists on three values
    12. check if dtheta consists of three values > 0
    13. check for 'CA' overlap basis
    14. check if cutoff value is >= 1.0
    15. check if zflag is 0 or 1
    16. check if clflag is 0 or 1
    17. check constraint file
        a.  file doesn't exist
        b.  file exists
    18. check constraint file parameters      
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
    19. check if low Rg cutoff is higher than high Rg cutoff
    20. check if Rg cutoffs are >= 0
        a.  low Rg cutoff is >= 0
        b.  high Rg cutoff is >= 0
    21. check if nexsegments1 and nexsegments2 are >= 0
        a.  nexsegments1 >= 0
        b.  nexsegments2 >= 0
    22. check if number nsegment1 values matches nexsegments1
    23. check if number nsegment2 values matches nexsegments2
    24. check if the number of reslow1 and numcont1 values matches nexsegments1
        a.  len(reslow1) = nexsegments1
        b.  len(numcont1) = nexsegments1
    25. check if the number of reslow2 and numcont2 values matches nexsegments2
        a.  len(reslow2) = nexsegments2
        b.  len(numcont2) = nexsegments2

    Use cases not tested:

    1.  check if pdbfile1 and pdbfile2 have missing residue (numbers)
    2.  check high and low residue values for pdbfile1 and pdbfile2
        a.  pdbfile contains the low and high residues listed for each segment
        b.  number of contiguous residues is > 0

    '''

    def setUp(self):

        gui_mimic_two_body_grid.test_variables(self, paths)


    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_two_body_grid.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if path exists
        '''
        self.path = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_two_body_grid.run_module(
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
        return_error = gui_mimic_two_body_grid.run_module(
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
        return_error = gui_mimic_two_body_grid.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path +
                          '  [code = TrueTrueFalse]', 'write permission not allowed']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if pdbmol1 exists
        '''
        self.pdbmol1 = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        pdbfile1=self.path+'/'+self.pdbmol1
        expected_error = ['file : '+pdbfile1+' does not exist']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if pdbmol1 is a valid pdb file
        '''
        self.pdbmol1 = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        pdbfile1=self.path+'/'+self.pdbmol1
        expected_error = ['input pdb file, ' +
                          pdbfile1[3:] + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if pdbmol2 exists
        '''
        self.pdbmol2 = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        pdbfile2=self.path+'/'+self.pdbmol2
        expected_error = ['file : '+pdbfile2+' does not exist']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if pdbmol2 is a valid pdb file
        '''
        self.pdbmol2 = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        pdbfile2=self.path+'/'+self.pdbmol2
        expected_error = ['input pdb file, ' +
                          pdbfile2[3:] + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if accpos is 0 or 1
        '''
        self.accpos = '2'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['accept supplied position needs to be (0==no or 1==yes) : 2']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if pos has a length of 3
        '''
        self.pos = '-20,-20'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['three float values are required for initial position (x,y,z) : [-20.0, -20.0]']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if trans has a length of 3
        '''
        self.trans = '2,2'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['three int values are required for number of x,y,z moves : [2, 2]']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if trans consists of 3 values > 0
        '''
        self.trans = '0,1,1'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['you must specifiy at least ONE translational "move" for each axis : [0, 1, 1]']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if trans consists of 3 values > 0
        '''
        self.trans = '1,0,1'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['you must specifiy at least ONE translational "move" for each axis : [1, 0, 1]']
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test if trans consists of 3 values > 0
        '''
        self.trans = '1,1,0'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['you must specifiy at least ONE translational "move" for each axis : [1, 1, 0]']
        assert_equals(return_error, expected_error)
              
    def test_15(self):
        '''
        test if dtrans has a length of 3
        '''
        self.dtrans = '2,2'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['three float values are required for dx,dy,dz step sizes : [2.0, 2.0]']
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if theta has a length of 3
        '''
        self.theta = '2,2'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['three int values are required for theta angular moves : [2, 2]']
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test if theta consists of 3 values > 0
        '''
        self.theta = '0,1,1'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['you must specifiy at least ONE angular "move" for each axis : [0, 1, 1]']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if theta consists of 3 values > 0
        '''
        self.theta = '1,0,1'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['you must specifiy at least ONE angular "move" for each axis : [1, 0, 1]']
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test if theta consists of 3 values > 0
        '''
        self.theta = '1,1,0'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['you must specifiy at least ONE angular "move" for each axis : [1, 1, 0]']
        assert_equals(return_error, expected_error)


    def test_20(self):
        '''
        test if dtheta has a length of 3
        '''
        self.dtheta = '2'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['three float values are required for dtheta (x,y,z) step sizes : [2.0]']
        assert_equals(return_error, expected_error)
        
    def test_21(self):
        '''
        test if basis is 'CA'
        '''
        self.basis = 'heavy'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['only "CA" is accepted as a basis']
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test if cutoff >= 1
        '''
        self.cutoff = '0.8'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['use a larger cutoff value, cutoff = 0.8']
        assert_equals(return_error, expected_error)
                               
    def test_23(self):
        '''
        test if Z coordinate filter selection is 0 or 1
        '''
        self.zflag = '2'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR in Z coordinate filter selection: zflag == 0 for "no" and 1 for "yes", zflag = 2']
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test if atomic constraints selection is 0 or 1
        '''
        self.cflag = '2'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR in atomic constraints selection: cflag == 0 for "no" and 1 for "yes", cflag = 2']
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test if constraint file exists
        '''
        self.cflag = '1'
        self.confile = './does_not_exist.txt'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['ERROR in constraint filename selection: file : ./does_not_exist.txt does not exist']
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test for bad seg1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_seg1.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 segment INT0 listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test for bad seg2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_seg2.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 segment INT5 listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test for bad atom1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_atom1.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 atom name XA listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)             

    def test_29(self):
        '''
        test for bad atom2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_atom2.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 atom name ZA listed in constraint file is not in your PDB file"]
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test for bad distance in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_distance.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 distance value is not appropriate: -100.0"]
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test for no distance in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'no_distance.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 no distance specified or error in line: COM"]
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test for COM or ATM type1 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_type1.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 TYPE1 is not valid (ATM OR COM): CON"]
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test for COM or ATM type2 in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_type2.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 2 TYPE2 is not valid (ATM OR COM): ATN"]
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test for two types COM and/or ATM in constraint file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'no_type2.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : LINE 1 Two type definitions are required for each constraint (ATM OR COM)"]
        assert_equals(return_error, expected_error)


    def test_35(self):
        '''
        test for second resid1 value equal or less than first
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_resid1.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid values in constraint file for constraint 1 are incorrect: second value is equal or less than first"]
        assert_equals(return_error, expected_error)

    def test_36(self):
        '''
        test for second resid2 value equal or less than first
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'bad_resid2.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid values in constraint file for constraint 0 are incorrect: second value is equal or less than first"]
        assert_equals(return_error, expected_error)

    def test_37(self):
        '''
        test if first value in first resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid1_first.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 0 is not in segment INT2"]
        assert_equals(return_error, expected_error)

    def test_38(self):
        '''
        test if second value in first resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid1_second.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 290 is not in segment INT1"]
        assert_equals(return_error, expected_error)

    def test_39(self):
        '''
        test if first value in second resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid2_first.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 0 is not in segment INT4"]
        assert_equals(return_error, expected_error)

    def test_40(self):
        '''
        test if second value in second resid range is in pdb file
        '''
        self.cflag = '1'
        self.confile = os.path.join(module_data_path,'missing_resid2_second.txt')
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = [" : resid 290 is not in segment INT3"]
        assert_equals(return_error, expected_error)

    def test_41(self):
        '''
        test for low Rg cutoff higher than high Rg cutoff
        '''
        self.lowrg = '401.0'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['low Rg cutoff is larger than high Rg cutoff, lowrg = 401.0 highrg = 400.0']
        assert_equals(return_error, expected_error)

    def test_42(self):
        '''
        test if low Rg cutoff > 0
        '''
        self.lowrg = '-1.0'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Rg cutoffs need to be >= zero, lowrg = -1.0 highrg = 400.0']
        assert_equals(return_error, expected_error)

    def test_43(self):
        '''
        test if high Rg cutoff > 0
        '''
        self.lowrg = '-5.0'
        self.highrg = '-1.0'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Rg cutoffs need to be >= zero, lowrg = -5.0 highrg = -1.0']
        assert_equals(return_error, expected_error)

    def test_44(self):
        '''
        test if number of nexsegments1 >= 0
        '''

        self.nexsegments1 = '-1'
        self.nsegments1 = 'INT1,INT2,LED1'
        self.reslow1 = '271, 271, 94'
        self.numcont1 = '18, 18, 31'        
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of excluded segments needs to be >= 0 (mol1) : -1']
        assert_equals(return_error, expected_error)

    def test_45(self):
        '''
        test if number of nexsegments2 >= 0
        '''

        self.nexsegments2 = '-1'
        self.nsegments2 = 'INT3,INT4,LED2'
        self.reslow2 = '271, 271, 94'
        self.numcont2 = '18, 18, 31'        
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of excluded segments needs to be >= 0 (mol2) : -1']
        assert_equals(return_error, expected_error)

    def test_46(self):
        '''
        test if number of nsegments1 
        '''

        self.nexsegments1 = '3'
        self.nsegments1 = 'INT1,INT2,LED1,LED2'
        self.reslow1 = '271, 271, 94'
        self.numcont1 = '18, 18, 31'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of segment names does not match number of excluded segments (mol1) : INT1,INT2,LED1,LED2']
        assert_equals(return_error, expected_error)

    def test_47(self):
        '''
        test if number of nsegments2 
        '''

        self.nexsegments2 = '3'
        self.nsegments2 = 'INT3,INT4,LED1,LED2'
        self.reslow2 = '271, 271, 94'
        self.numcont2 = '18, 18, 31'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of segment names does not match number of excluded segments (mol2) : INT3,INT4,LED1,LED2']
        assert_equals(return_error, expected_error)

    def test_48(self):
        '''
        test if number of lowres values matches the number of nexsegments1
        '''

        self.nexsegments1 = '3'
        self.nsegments1 = 'INT1,INT2,LED1'
        self.reslow1 = '2,5,5,2'
        self.numcont1 = '18, 18, 31'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of low residue values does not match the number of excluded segments (mol1), lowres1 = [2, 5, 5, 2] nexsegments1 = 3']
        assert_equals(return_error, expected_error)

    def test_49(self):
        '''
        test if number of lowres values matches the number of nexsegments2
        '''

        self.nexsegments2 = '3'
        self.nsegments2 = 'INT3,INT4,LED2'
        self.reslow2 = '2,5,5,2'
        self.numcont2 = '18, 18, 31'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of low residue values does not match the number of excluded segments (mol2), lowres2 = [2, 5, 5, 2] nexsegments2 = 3']
        assert_equals(return_error, expected_error) 

    def test_50(self):
        '''
        test if number of contiguous residue values matches the number of nexsegments1
        '''

        self.nexsegments1 = '3'
        self.nsegments1 = 'INT1,INT2,LED1'
        self.reslow1 = '271,271,94'
        self.numcont1 = '18,18,31,10'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of contiguous residues does not match the number of excluded segments (mol1), numcont1 = [18, 18, 31, 10] nexsegments1 = 3']
        assert_equals(return_error, expected_error)

    def test_51(self):
        '''
        test if number of contiguous residue values matches the number of nexsegments2
        '''


        self.nexsegments2 = '3'
        self.nsegments2 = 'INT3,INT4,LED2'
        self.reslow2 = '271,271,94'
        self.numcont2 = '18,18,31,10'
        return_error = gui_mimic_two_body_grid.run_module(self, test_filter=True)

        ''' check for file error '''
        expected_error = ['the number of contiguous residues does not match the number of excluded segments (mol2), numcont2 = [18, 18, 31, 10] nexsegments2 = 3']
        assert_equals(return_error, expected_error)       

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
