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
    4.  check if trajectory file has nonzero number of frames
    5.  check if nsteps is > 0
    6.  check if temperature is >= 0
    7.  check if rgvalue >= 0
    8.  check if rgforce >= 0
        a.  rgforce = 0
        b.  rgforce > 0
            i. rgvalue also must be > 0
    9.  check keepout is 0 or 1
    10. check if pretamd_min_step >= 0
    11. check if number_flexible_segments >= 1    
    12  check dna_segnames
        a.  len(dna_segnames) > 0
            i. check if dna_segname is in PDB file
    13. check if topology file exists
    14. check if parameter file exists
    15. check if CHARMM executable file exists
    16. check if dcdfreq is <= nsteps
    17. check if dcdfreq is a divisor of nsteps
    18. check if pretamd_min_steps is an integer
    19. check flexible segment variables
        a. moltype is 'dna' when len(dna_segnames) > 0
        b. number of flexible segments >= 1 
        c. number of residue ranges = number of flexible regions
        d. number of ranges for each segment is >=1
        e. segment name is in pdb file
        f. low residue can't include n-terminus
        g. low residue + number of contiguous residues can't include c-terminus
        h. pdb file contains the low and high flexible residues listed for each flexible region 
        i. residues in each flexible regions are in pdb file          
        j. residue ranges don't overlap



    '''

    def setUp(self):

        gui_mimic_torsion_angle_md.test_variables(self, paths)


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


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
