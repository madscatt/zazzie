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
import sassie.simulate.monomer_monte_carlo.gui_mimic_monomer_monte_carlo as gui_mimic_monomer_monte_carlo
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
        d.  low residue is lower than the n-terminal amino acid number  #This is redundant with lowres test 19a since aa number will not be in pdb file.  Not tested.
        e.  residue range doesn't exceed the number of amino acids      #This is redundant with test 19a since aa number will not be in pdb file.  Not tested.
        f.  residue values increase from low to high
        g.  residue ranges don't overlap with alignment range 
    20. check alignment residue values
        a.  pdb file contains low and high residue amino acids
        b.  alignment range isn't too small (less than 3 points)
    21. check if directed Monte Carlo value is >=0
    22. check overlap basis atoms
        a.  basis cannot be parsed                      #NOT YET TESTED  didn't find a string that can't be parsed that didn't already trigger another error
        b.  basis can be parsed
            i. check that atom name is in PDB file
            ii.check that atom has VDW paramters        #NOT TESTED  there are no atoms in vdw list that don't have vdw parameters
   
   
    '''

    def setUp(self):

        gui_mimic_monomer_monte_carlo.test_variables(self, paths)


    def test_49(self):
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
