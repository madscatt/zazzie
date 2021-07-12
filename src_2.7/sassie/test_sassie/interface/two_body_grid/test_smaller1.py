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
    3.  check pdbmol1
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    4.  check pdbmol2
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

    1.  check if pdbmol1 and pdbmol2 have missing residue (numbers)
    2.  check high and low residue values for pdbmol1 and pdbmol2
        a.  pdbfile contains the low and high residues listed for each segment
        b.  number of contiguous residues is > 0

    '''

    def setUp(self):

        gui_mimic_two_body_grid.test_variables(self, paths)

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
