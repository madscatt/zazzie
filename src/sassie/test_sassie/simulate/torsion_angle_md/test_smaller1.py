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
import locale
import shutil
import numpy
import multiprocessing

import sasmol.sasmol as sasmol
import sassie.simulate.torsion_angle_md.gui_mimic_torsion_angle_md as gui_mimic_torsion_angle_md
#import gui_mimic_torsion_angle_md as gui_mimic_torsion_angle_md

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'simulate', 'torsion_angle_md') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_torsion_angle_md(MockerTestCase):

    '''
	TAMD is the module that contains the functions
    	that are used to run a series of tamd dynamics calculations
	on a set of structures in a supplied pdb/dcd file.

        REFERENCES:

        J. Chen et al.
        Journal of Computational Chemistry  26  1565-1578  (2005)

        W. Zhang et al.
        Journal of Molecular Graphics and Modeling  73  179-190  (2017)        

        INPUT:  variable descriptions:
              
        reference PDB file
        input PDB or DCD file
        number of preliminary minimization steps
        number of TAMD steps
        poll frequency
        path and name of topology file
        path and name of parameter file
        path and name of charmm executable file
        output (DCD) file
        keep output files (0==no, 1==yes)
        temperature (K)
        Rg force
        Rg value (Angstroms)
        number of flexible segments
        segvariables:   flexible segment names, number of ranges, low residues, contiguous residues, molecule type
        dsDNA segment names (if any)
        
        Advanced input:

        frequency to save individual DCD files 

        
        OUTPUT:
                
        files stored in ~/run_name/torsion_angle_md directory:
                
        original PDB file with molecular structure data
        tamd_output.pdb (PDB file created by TAMD)
        tamd_output.psf (PSF file created by TAMD)
        TAMD tree file
        TAMD restart file
        DCD file containing TAMD trajectories (size depends on DCD write frequency)
        tamd_dyn_00001.dcd file (trajectory files -- if keep output files is chosen)
        min_00001.out (output files) 
        input file for TAMD run (temp.inp)
        temp_0.pdb (temporary PDB files)
        


     Use cases:

    1.  input
        a.  one frame input file
            i.  input file is a PDB file
            ii.  input file is a DCD file
        b.  multiple frames in input file
            i.  input file is a PDB file
            ii.  input file is a DCD file   

    2.  flexible segments
        a.  one flexible segment
            i.  one flexible region
            ii. more than one flexible region
        b.  more than one flexible segment
            i.  one flexible region per segment
            ii. more than one flexible region per segment

    3.  molecule type
        a.  protein
        b.  rna
        c.  dna

    4.  options:
        a.  Rg force and Rg value !=0
        b.  keep output files

        

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

    path                            string      input path (not in variables)
    poll_frequency                  float       time used in time.sleep command (not input by user)

    Test tree:

    project name

***************************
    one flexible segment
***************************

    reference PDB     reference PDB      reference PDB      reference PDB      reference PDB    reference PDB    *reference PDB file
    input PDB         input PDB          input PDB          input DCD          input DCD        input DCD        *trajectory file
    output DCD        output DCD         output DCD         output DCD         output DCD       output DCD       *output file
    protein           rna                ssdna              protein            rna              ssdna

                                                   no options                                                         
                                                   keep output files
                                                   Rg force and Rg value != 0
                                                   multiple flexible regions

    
********************************************************************************
    multiple flexible segments (protein and dsdna only -- these tests take time)
********************************************************************************

    reference PDB        reference PDB        reference PDB        reference PDB     *reference PDB file
    input PDB            input PDB            input DCD            input DCD         *trajectory file
    output DCD           output DCD           output DCD           output DCD        *output file
    protein              dsdna                protein              dsdna

                                                   no options                                                         
                                                   keep output files
                                                   Rg force and Rg value != 0
                                                   multiple flexible regions per segment
                                                                                                                                                            

****************************************************************************************
    multiple frames in input file (two cases chosen -- these tests take several minutes) 
****************************************************************************************

    reference PDB                   reference PDB
    input DCD                       input PDB
    output DCD                      output DCD
    protein/dsdna                   protein/protein
    keep output files               multiple flexible regions per segment
    Rg force and Rg value != 0      

output dcd and pdb files will not be the same; testing for completion and correct number of frames in output DCD file only
    '''


    module = 'torsion_angle_md'

    def setUp(self):

       gui_mimic_torsion_angle_md.test_variables(self, paths)


    def check_for_completion(self):
        '''
        tests for completion of run without an error
        looks for junk.out file in current directory, which won't exist unless there was an error
        '''

        test = os.path.isfile('junk.out')
        if(test):
            return False
        return True

    def check_charmm_input_file(self, testfile, txt):
        '''
        tests for "set rgforce" in temp.inp file
        used when Rg force is > 0
        '''

        def check(testfile, txt):
            with open(testfile) as dataf:
                return any(txt in line for line in dataf)

        if check(testfile, txt):
            return True
        else:
            return False


    def test_41(self):
        '''
        multiple frame DCD input, multiple flexible segments, protein/dsdna, multiple flexible regions per segment
        '''

        self.pdbfile = os.path.join(pdb_data_path,'c36_w601_ncp_min.pdb')
        self.infile = os.path.join(dcd_data_path,'c36_w601_ncp_2_frames.dcd')
        self.outfile = 'c36_w601_ncp_2_frames.dcd'
        self.number_flexible_segments = '3'
        self.all_flexible_segnames=['1H3','DNA1','DNA2']
        self.all_snumranges=['2','2','2']
        self.all_srlow=['2,10','16,151','194,329']       
        self.all_srnum=['3,30','9,9','9,9']
        self.all_moltype=['protein','dna','dna']
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            if self.all_moltype[i] == 'dna':
                self.dna_segnames += self.all_flexible_segnames[i] + ','
        if self.dna_segnames and self.dna_segnames[-1] ==',':
            self.dna_segnames = self.dna_segnames[:-1]
        self.psegvariables = []
        for i in xrange(locale.atoi(self.number_flexible_segments)):
            self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
        print 'dna_segnames in test: ', self.dna_segnames
        gui_mimic_torsion_angle_md.run_module(self)

        '''check for completion'''
        assert_equals(self.check_for_completion(),True)

        '''check for proper number of frames in output DCD file'''
        outfile = os.path.join(self.runname, self.module, self.outfile)
        molecule = sasmol.SasMol(0)
        molecule.read_pdb(self.pdbfile)
        molecule.read_dcd(outfile)
        dcdfile = molecule.open_dcd_read(outfile)
        nf = dcdfile[2]
        print 'number of frames in output DCD file: ', nf
        assert_equals(nf,2)               


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__=='__main__':
    main()

