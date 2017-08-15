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
import string
import sassie.tools.contrast_calculator.gui_mimic_contrast_calculator as gui_mimic_contrast_calculator
#import gui_mimic_contrast_calculator as gui_mimic_contrast_calculator

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'interface', 'contrast_calculator') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Contrast_Calculator_Filter(MockerTestCase):

    '''
    System integration test for contrast_calculator_filter.py / sassie 2.0

    Test to see whether contrast_calculator_filter catches improper input.

    Inputs tested:

            runname:                            project name
            inpath:                             path name for input files (pdb or sequence)
            outfile:                            output filename
            numfiles:                           number of input files (protein, rna or dna)
            solute_conc:                        concentration of solute
            d2ostep:                            step in fraction D2O used when calculating SLD, contrast and I(0)
            fexchp:                             fraction of exchangeable hydrogen atoms that exchange for protein components
            fexchn:                             fraction lf exchangeable hydrogen atoms that exchange for nucleic acid (rna, dna) components
            seqfiles:                           names of sequence or pdb files
            numunits:                           number of times the sequence from seqfile is used in the protein, rna and/or dna complex
            fracdeut:                           fraction of deuteration for the subunit represented in the seqfile
            moltype:                            type of molecule represented in the seqfile (protein, rna, or dna)
            isFasta:                            indicates whether the file is a FASTA file (0=PDB, 1=FASTA)
            plotflag:                           flag to indicate whether to plot data (0=no plot, 1=plot)
            numsolv:                            number of non-water components in the solvent
            solv_comp:                          chemical formula representing the non-water solvent component
            solv_conc:                          molar concentration of the non-water solvent component
            number_of_chemicals:                number of non-protein (rna or dna) components in the solute
            formula_array:                      chemical formulas for the solute components
            number_exchangeable_hydrogens:      the number of exchangeable hydrogen atoms in the solute component
            fraction_exchangeable_hydrogens:    fraction of exchangeable hydrogens that exchange in the solute compoent
            mass_density:                       mass density of solute component   


    Use cases tested:

    1.  check numfiles
        a.  numfiles >= 0
        b.  numfiles < 0
    2.  check numsolv
        a.  numsolv >= 0
        b.  numsolv < 0
    3.  check number_of_chemicals
        a.  number_of_chemicals >=0
        b.  number_of_chemicals < 0    
    4.  check contrast variables
        a.  numfiles + numformlas >= 1
        b.  runname is not blank
        c.  outfile name is not blank
        e.  solute concentration is > 0
        f.  D2O step is between 1 and 100
        g.  D2O step can divide 100% evenly
        h.  fraction of exchangeable protein hydrogens is between 0.0 and 1.0
        i.  fraction of exchangeable nucleic acid hydrogens is between 0.0 and 1.0         
    5.  check input file path permissions        
        a.  no permission error
        b.  permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed
    6.  check input file                            
        a.  input file doesn't exist
            i.  input file name is blank
            ii. input file name is not blank
        b.  input file exists
            i.  input file is not a valid PDB file
            ii. input file is a valid PDB file
                a.  PDB file doesn't contain non-protein, rna or dna residues
                b.  PDB file doesn't contain residues that don't match the molecule type
            iii.input file is a valid FASTA file
                a.  FASTA file is not empty
                b.  FASTA file doesn't contain residues that don't match the molecule type
    7.  check molecule type
        a.  molecule type is protein, rna or dna
        b.  molecule type is not protein, rna or dna
    8.  check number of units
        a.  number of units is an integer
        b.  number of units is not an integer
        c.  number of units > 1
        d.  number of units < 1
    9.  check deuteration fraction
        a.  deuteration fraction is a number
        b.  deuteration fraction is not a number
        c.  deuteration fraction is between 0.0 and 1.0
        d.  deuteration fraction is not between 0.0 and 1.0
    10. check number of exchangeable hydrogens
        a.  number of exchangeable hydrogens is in integer
        b.  number of exchangeable hydrogens is not an integer
        c.  number of exchangeable hydrogens >= 0
        d.  number of exchangeable hydrogens < 0
    11. check fraction of exchangeable hydrogens  
        a. fraction of exchangeable hydrogens is a number
        b. fraction of exchangeable hydrogens is not a number
        c. fraction of exchangeable hydrogens is between 0.0 and 1.0
        d. fraction of exchangeable hydrogens is not between 0.0 and 1.0
    12. check mass density
        a.  mass density is a number
        b.  mass density is not a number
        c.  mass density > 0
        d.  mass density <=0
    13. check solvent component concentration
        a.  solvent component concentration is a number
        b.  solvent component concentration is not a number
        c.  solvent component concentration > 0
        d.  solvent component concentration <= 0  

             
    Not tested:

    numfiles, numsolv and number_of_chemicals are integers   (caught in input filter)
    solute_conc is a number (caugt in input filter)
    D2O step is an integer (caught in input filter)
    fraction of exchangeable protein and nucleic acid hydrogens are numbers (caught in input filter) 
    ivariables exist (for sassie_gui?)

    '''

    def setUp(self):

        gui_mimic_contrast_calculator.test_variables(self, paths)


        
    def test_19(self):
        '''
        test if input file exists (second input file of two doesn't exist)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1']        
        self.inpath = other_data_path
        self.seqfiles = ['protein_sequence.txt', 'non_existent.txt']
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input pdb file "' +
                          os.path.join(self.inpath,self.seqfiles[1]) + '" does not exist.']
        assert_equals(return_error, expected_error)



    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)
            
if __name__ == '__main__':
    main()
