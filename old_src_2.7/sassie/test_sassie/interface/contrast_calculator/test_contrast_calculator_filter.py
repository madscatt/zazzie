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


    def test_1(self):
        '''
        test if numfiles is >= 0
        '''
        self.numfiles = '-1'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['The number of input files must be greater than or equal to zero.']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if numsolv is >= 0
        '''
        self.numsolv = '-1'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['The number of solvent components must be greater than or equal to zero.']
        assert_equals(return_error, expected_error)

    def test_3(self):
        '''
        test if number_of_chemicals is >= 0
        '''
        self.number_of_chemicals = '-1'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['The number of additional components must be greater than or equal to zero.']
        assert_equals(return_error, expected_error)

    def test_4(self):
        '''
        test if numfiles, numsolv and number_of_chemicals is >= 0
        '''
        self.number_of_chemicals = '-1'
        self.numsolv = '-1'
        self.numfiles = '-1'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['The number of input files must be greater than or equal to zero.',
        'The number of solvent components must be greater than or equal to zero.',
        'The number of additional components must be greater than or equal to zero.']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if numfiles + number_of_chemicals >= 1
        '''
        self.numfiles = '0'
        self.number_of_chemicals = '0'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ["You must enter either one file (PDB or Fasta) or one chemical formula"]
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if runname is blank
        '''
        self.runname = ''
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['Please enter a project name']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if outfile is blank
        '''
        self.outfile = ''
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['Please enter an output filename']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if solute concentration > 0
        '''
        self.solute_conc = '-1'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['Solute concentration must be greater than 0.']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if d2ostep is between 0 and 100
        '''
        self.d2ostep = '110'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['D2O step must be an integer between 1 and 100']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if d2ostep divides 100 evenly
        '''
        self.d2ostep = '3'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['D2O step size must divide 100% evenly (1%, 2%, 5%, etc.)']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if fexchp is between 0.0 and 1.0
        '''
        self.fexchp = '1.1'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['The fraction of exchangeable protein hydrogens must be between 0.0 and 1.0.']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if fexch is between 0.0 and 1.0
        '''
        self.fexchn = '1.1'
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['The fraction of exchangeable nucleic acid hydrogens must be between 0.0 and 1.0.']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if input file path exists
        '''
        self.inpath = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_contrast_calculator.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['Permission error in input file path ' + self.inpath + ':  [code = FalseFalseFalse]',
                          'Path does not exist.']
        assert_equals(return_error, expected_error)

    def test_14(self):
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

        self.inpath = os.path.join('./', 'empty_folder')
        return_error = gui_mimic_contrast_calculator.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['Permission error in input file path ' +
                          self.inpath + ':  [code = TrueFalseTrue]', 'Read permission not allowed.']
        assert_equals(return_error, expected_error)

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')        

    def test_15(self):
        '''
        test if directory has write permission
        '''

        ''' make a directory '''
        os.system('mkdir empty_folder1')
        ''' see if you can write to the directory '''
#        print os.access('empty_folder1', os.W_OK)
        ''' make the directory un-writeable'''
        os.system('chmod a-w empty_folder1')
        ''' see if you can write to the directory '''
#        print os.access('empty_folder', os.W_OK)

        self.inpath = os.path.join('./', 'empty_folder1')
        return_error = gui_mimic_contrast_calculator.run_module(
            self, file_check=True)
#        print 'return_error: ', return_error

        ''' check for path error '''
        expected_error = ['Permission error in input file path ' +
                          self.inpath + ':  [code = TrueTrueFalse]', 'Write permission not allowed.']
#        print 'expected_error: ', expected_error
        assert_equals(return_error, expected_error)

        ''' make the directory writeable'''
        os.system('chmod a+w empty_folder1')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder1')        

        
    def test_16(self):
        '''
        test if directory has write permission
        '''
        self.inpath = os.path.join(module_data_path, 'no_write_permission')
        return_error = gui_mimic_contrast_calculator.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['Permission error in input file path ' + self.inpath +
                          ':  [code = TrueTrueFalse]', 'Write permission not allowed.']
        assert_equals(return_error, expected_error)

    def test_17(self):
        '''
        test if input file is not blank
        '''
        self.inpath = module_data_path
        self.seqfiles = ['']
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['No input file has been specified on line 1.']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if input file exists (first input file of two doesn't exist)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1']        
        self.inpath = other_data_path
        self.seqfiles = ['non_existent.txt', 'dna_sequence.txt']
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input pdb file "' +
                          os.path.join(self.inpath,self.seqfiles[0]) + '" does not exist.']
        assert_equals(return_error, expected_error)
        
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

    def test_20(self):
        '''
        test if input file is a valid PDB file (first of two is not valid)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'protein']
        self.isFasta = ['0', '0']        
        self.inpath = module_data_path
        self.seqfiles = ['not_valid.pdb','skp_trimer.pdb']
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file "'+os.path.join(self.inpath,self.seqfiles[0])+'" is not in valid PDB format. FASTA?']
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test if input file is a valid PDB file (second of two is not valid)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'protein']
        self.isFasta = ['0', '0']        
        self.inpath = module_data_path
        self.seqfiles = ['skp_trimer.pdb','not_valid.pdb']
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file "'+os.path.join(self.inpath,self.seqfiles[1])+'" is not in valid PDB format. FASTA?']
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test if input PDB file has valid molecule type (first of two has no molecule type)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['', 'protein']
        self.isFasta = ['0', '0'] 
        self.seqfiles = ['skp_trimer.pdb', 'ompA']               
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 1: Select a molecule type']
        assert_equals(return_error, expected_error)

    def test_23(self):
        '''
        test if input PDB file has valid molecule type (second of two has no molecule type)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', '']
        self.isFasta = ['0', '0'] 
        self.seqfiles = ['skp_trimer.pdb', 'ompA.pdb']               
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 2: Select a molecule type']
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test for empty FASTA file (first of two is empty)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['empty_fasta.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The FASTA file "'+os.path.join(self.inpath,self.seqfiles[0])+'" appears to be empty.']
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test for empty FASTA file (second of two is empty)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'empty_fasta.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The FASTA file "'+os.path.join(self.inpath,self.seqfiles[1])+'" appears to be empty.']
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test for FASTA file with invalid molecule type (first of two has no molecule type)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The molecule type for input file 1 must be DNA, RNA or protein.']
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test for FASTA file with invalid molecule type (second of two has no molecule type)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', '']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The molecule type for input file 2 must be DNA, RNA or protein.']
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test FASTA file for non-rna residue 
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'rna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'wrong_residue.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 2 contains the non-rna residue T.']
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test FASTA file for non-protein residue 
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'protein']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'wrong_residue.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 2 contains the non-protein residue U.']
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test FASTA file for non-dna residue 
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'wrong_residue.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 2 contains the non-dna residue U.']
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test for non-integer number of units (first of two is not an integer)
        '''
        self.numfiles = '2'
        self.numunits = ['1.5', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of units of each molecule for input file 1 must be an integer.']
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test for non-integer number of units (second of two is not an integer)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1.5']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of units of each molecule for input file 2 must be an integer.']
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test if number of units > 1 (first of two is < 1)
        '''
        self.numfiles = '2'
        self.numunits = ['0', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of units of each molecule for input file 1 must be at least 1.']
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test if number of units > 1 (second of two is < 1)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '0']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of units of each molecule for input file 2 must be at least 1.']
        assert_equals(return_error, expected_error)

    def test_35(self):
        '''
        test if deuteration fraction is a number (first of two is not a number)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['f', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The deuteration fraction of input file 1 must be a number.']
        assert_equals(return_error, expected_error)

    def test_36(self):
        '''
        test if deuteration fraction is a number (second of two is not a number)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['1', 'f']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The deuteration fraction of input file 2 must be a number.']
        assert_equals(return_error, expected_error)

    def test_37(self):
        '''
        test if deuteration fraction is between 0.0 and 1.0 (first of two is not between 0.0 and 1.0)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['2', '1']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The deuteration fraction of input file 1 must be between 0.0 and 1.0']
        assert_equals(return_error, expected_error)

    def test_38(self):
        '''
        test if deuteration fraction is between 0.0 and 1.0 (second of two is not between 0.0 and 1.0)
        '''
        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['1', '2']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['1', '1'] 
        self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']              
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The deuteration fraction of input file 2 must be between 0.0 and 1.0']
        assert_equals(return_error, expected_error)

    def test_39(self):
        '''
        test if the number of exchangeable hydrogens is an integer (first of two is not an integer)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12.1', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of exchangeable hydrogens for additional component 1 must be an integer.']
        assert_equals(return_error, expected_error)

    def test_40(self):
        '''
        test if the number of exchangeable hydrogens is an integer (second of two is not an integer)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5.1']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of exchangeable hydrogens for additional component 2 must be an integer.']
        assert_equals(return_error, expected_error)

    def test_41(self):
        '''
        test if the number of exchangeable hydrogens is >= 0 (first of two is less than 0)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['-1', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of exchangeable hydrogens for additional component 1 must be greater than or equal to zero.']
        assert_equals(return_error, expected_error)

    def test_42(self):
        '''
        test if the number of exchangeable hydrogens is >= 0 (second of two is less than 0)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '-1']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The number of exchangeable hydrogens for additional component 2 must be greater than or equal to zero.']
        assert_equals(return_error, expected_error)

    def test_43(self):
        '''
        test if the fraction of exchangeable hydrogens is a number (first of two is not a number)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['0.(5', '0.45']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The fraction of exchangeable hydrogens for additional component 1 must be a number.']
        assert_equals(return_error, expected_error)

    def test_44(self):
        '''
        test if the fraction of exchangeable hydrogens is a number (second of two is not a number)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.4%']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The fraction of exchangeable hydrogens for additional component 2 must be a number.']
        assert_equals(return_error, expected_error)

    def test_45(self):
        '''
        test if the fraction of exchangeable hydrogens is a between 0.0 and 1.0 (first of two is not between 0.0 and 1.0)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['1.95', '0.45']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The fraction of exchangeable hydrogens for additional component 1 must be between 0.0 and 1.0.']
        assert_equals(return_error, expected_error)

    def test_46(self):
        '''
        test if the fraction of exchangeable hydrogens is a between 0.0 and 1.0 (second of two is not between 0.0 and 1.0)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '-0.45']
        self.mass_density = ['1.1', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The fraction of exchangeable hydrogens for additional component 2 must be between 0.0 and 1.0.']
        assert_equals(return_error, expected_error)

    def test_47(self):
        '''
        test if the mass density is a number (first of two is not a number)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The mass density of additional component 1 must be a number.']
        assert_equals(return_error, expected_error)

    def test_48(self):
        '''
        test if the mass density (second of two is not a number)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['1.1', '1,3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The mass density of additional component 2 must be a number.']
        assert_equals(return_error, expected_error)

    def test_49(self):
        '''
        test if the mass density is > 0 (first of two is not > 0)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['0', '1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The mass density of additional component 1 must be greater than 0.']
        assert_equals(return_error, expected_error)

    def test_50(self):
        '''
        test if mass density > 0 (second of two is not > 0)
        '''

        self.number_of_chemicals = '2'
        self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']   
        self.number_exchangeable_hydrogens = ['12', '5']
        self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
        self.mass_density = ['1.1', '-1.3']              
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The mass density of additional component 2 must be greater than 0.']
        assert_equals(return_error, expected_error)

    def test_51(self):
        '''
        test if the solvent component concentration is a number (first of two is not a number)
        '''

        self.numsolv = '2'
        self.solv_comp = ['NaCl','KCl']
        self.solv_conc = ['0,15','0.05'] 
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The concentration of solvent component 1 must be a number.']
        assert_equals(return_error, expected_error)

    def test_52(self):
        '''
        test if the solvent component concentration is a number (second of two is not a number)
        '''

        self.numsolv = '2'
        self.solv_comp = ['NaCl','KCl']
        self.solv_conc = ['0.15',''] 
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The concentration of solvent component 2 must be a number.']
        assert_equals(return_error, expected_error)

    def test_53(self):
        '''
        test if the solvent component concentration is > 0 (first of two is not >0)
        '''

        self.numsolv = '2'
        self.solv_comp = ['NaCl','KCl']
        self.solv_conc = ['-0.15','0.05'] 
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The concentration of solvent component 1 must be greater than 0.']
        assert_equals(return_error, expected_error)

    def test_54(self):
        '''
        test if the solvent component concentration is > 0 (second of two is not >0)
        '''

        self.numsolv = '2'
        self.solv_comp = ['NaCl','KCl']
        self.solv_conc = ['0.15','0.0'] 
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['The concentration of solvent component 2 must be greater than 0.']
        assert_equals(return_error, expected_error)

    def test_55(self):
        '''
        test PDB file for non-protein residue
        '''

        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'protein']
        self.isFasta = ['0', '0'] 
        self.seqfiles = ['skp_trimer.pdb', 'c36_dsDNA60_min.pdb']               
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 2 contains residues that do not match the moltype "protein"']
        assert_equals(return_error, expected_error)

    def test_56(self):
        '''
        test PDB file for non-DNA residue
        '''

        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'dna']
        self.isFasta = ['0', '0'] 
        self.seqfiles = ['skp_trimer.pdb', 'trunc2a_min.pdb']               
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 2 contains residues that do not match the moltype "dna"']
        assert_equals(return_error, expected_error)

    def test_57(self):
        '''
        test PDB file for non-RNA residue
        '''

        self.numfiles = '2'
        self.numunits = ['1', '1']
        self.fracdeut = ['0', '0']
        self.moltype = ['protein', 'rna']
        self.isFasta = ['0', '0'] 
        self.seqfiles = ['skp_trimer.pdb', 'ompA.pdb']               
        self.inpath = module_data_path
        return_error = gui_mimic_contrast_calculator.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['Input file 2 contains residues that do not match the moltype "rna"']
        assert_equals(return_error, expected_error)


    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)
            
if __name__ == '__main__':
    main()
