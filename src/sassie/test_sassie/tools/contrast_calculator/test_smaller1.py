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
import shutil
import numpy
import multiprocessing

import sasmol.sasmol as sasmol
import sassie.tools.contrast_calculator.gui_mimic_contrast_calculator as gui_mimic_contrast_calculator
#import gui_mimic_contrast_calculator as gui_mimic_contrast_calculator

import filecmp
from itertools import ifilter, izip
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'tools', 'contrast_calculator') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Contrast_Calculator(MockerTestCase):

    '''
    System integration test for contrast_calculator.py / sassie 2.0

    CONTRAST_CALCULATOR is the module that calculates SLD, I(0) and contrast v. %D2O 
	for a given molecule or complex based only on sequence or chemical formula.  It will
	read the sequence from a FASTA or PDB file.


    INPUT:

            runname:                            project name
            inpath:                             path name for input files (pdb or sequence)
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


    OUTPUT:
            outfile:              	output filename prefix

            outfile_contrast.txt
            outfile_sld.txt
            outfile_izero.txt

    Use cases:

    1.  Input file
        a.  input file is a PDB file
            i.  protein
            ii. DNA
            iii.RNA
        b.  input file is a FASTA file
            i.  protein
            ii. DNA
            iii.RNA
        c.  no input file (no protein, DNA or RNA components)
    2.  Additional components
        a.  there are no additional components
        b.  there are additional components defined by chemical formula, number of exchangeable hydrogens
                fraction of exchangeable hydrogens, mass density
    3.  Non-water solvent components
        a.  there are no non-water solvent components
        b.  there are non-water solvent components defined by chemical formula, concentration

    Options:
        One component
        More than one component

    Inputs tested:

            runname:                            string      project name
            inpath:                             string      path name for input files (pdb or sequence)
            outfile:                            string      output filename
            numfiles:                           integer     number of input files (protein, rna or dna)
            solute_conc:                        float       concentration of solute
            d2ostep:                            integer     step in fraction D2O (non-zero, 1 to 100, must divide 100 evenly)
            fexchp:                             float       fraction of exchangeable hydrogen atoms that exchange for protein components (0.0 to 1.0)
            fexchn:                             float       fraction lf exchangeable hydrogen atoms that exchange for nucleic acid components (0.0 to 1.0)
            seqfiles:                           string      names of sequence or pdb files
            numunits:                           integer     number of times the sequence from seqfile is used in the protein, rna and/or dna complex
            fracdeut:                           float       fraction of deuteration for the subunit represented in the seqfile (0.0 to 1.0)
            moltype:                            string      type of molecule represented in the seqfile (protein, rna, or dna)
            isFasta:                            integer     indicates whether the file is a FASTA file (0=PDB, 1=FASTA)
            numsolv:                            integer     number of non-water components in the solvent
            solv_comp:                          string      chemical formulas representing the non-water solvent components
            solv_conc:                          float       molar concentration of the non-water solvent component
            number_of_chemicals:                integer     number of non-protein (rna or dna) components in the solute
            formula_array:                      string      chemical formulas for the solute components
            number_exchangeable_hydrogens:      integer     the number of exchangeable hydrogen atoms in the solute component
            fraction_exchangeable_hydrogens:    float       fraction of exchangeable hydrogens that exchange in the solute compoent (0.0 to 1.0)
            mass_density:                       float       mass density of solute component (non-zero)  


    Test tree:

    project name
    input path (if applicable)

**************************
    one input file
**************************
    input PDB       input FASTA         input PDB       input FASTA     input PDB       input FASTA          
    protein         protein             DNA             DNA             RNA             RNA                          
                                no additinal components
                                no non-water solvent components           
    
    input PDB       input FASTA         input PDB       input FASTA     input PDB       input FASTA          
    protein         protein             DNA             DNA             RNA             RNA                          
                                no additinal components
                                two non-water solvent components                   

    input PDB       input FASTA         input PDB       input FASTA     input PDB       input FASTA          
    protein         protein             DNA             DNA             RNA             RNA                          
                                two additinal components
                                no non-water solvent components   

    input PDB       input FASTA         input PDB       input FASTA     input PDB       input FASTA          
    protein         protein             DNA             DNA             RNA             RNA                          
                                two additinal components
                                two non-water solvent components  

**************************
    no input files
**************************

    chemical formula                    chemical formula                    two chemical formulas               two chemical formulas
    no non-water solvent components     two non-water solvent components    no non-water solvent components     two non-water solvent components


**********************************
    multiple input files
**********************************

    input PDB       input FASTA     input FASTA
    protein         protein         protein
    input PDB       input FASTA     input PDB       
    protein         DNA             DNA     
            no additinal components
            no non-water solvents

    input PDB       input FASTA     input FASTA
    protein         protein         protein
    input PDB       input FASTA     input PDB       
    protein         DNA             DNA     
            two additinal components
            two non-water solvents
           
    

    '''

    module = 'contrast_calculator'

    def setUp(self):

       gui_mimic_contrast_calculator.test_variables(self, paths)



#   methods needed to check the output files ignoring the first line that contains the date
    def ignore_date(self,line):
        if line.startswith('#Date'):
            return False # ignore it
        return True

    def check_files(self,file1,file2):
        with open(file1) as f1, open(file2) as f2:
            f1 = ifilter(self.ignore_date, f1)
            f2 = ifilter(self.ignore_date, f2)
            value = all(x == y for x, y in izip(f1, f2))
        return value

    def test_1(self):
        '''
        test PDB input, protein, no additional components, no non-water solvent components
        '''

        self.inpath = pdb_data_path
        self.seqfiles = ['hiv1_gag.pdb']
        self.isFasta = ['0']  

        gui_mimic_contrast_calculator.run_module(self)

        ''' confirm output output files are correct '''
        testfile = os.path.join(self.runname, self.module,self.outfile+'_contrast.txt')
        correct_testfile = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_noad_nonw',self.outfile+'_contrast.txt')
        assert_equals(self.check_files(testfile,correct_testfile), True)                       
        testfile1 = os.path.join(self.runname, self.module,self.outfile+'_sld.txt')
        correct_testfile1 = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_noad_nonw',self.outfile+'_sld.txt')
        assert_equals(self.check_files(testfile1,correct_testfile1), True)
        testfile2 = os.path.join(self.runname, self.module,self.outfile+'_izero.txt')
        correct_testfile2 = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_noad_nonw',self.outfile+'_izero.txt')
        assert_equals(self.check_files(testfile2,correct_testfile2), True)        
                       

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

