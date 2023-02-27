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
import string
import shutil
import numpy
import multiprocessing

import sasmol.sasmol as sasmol
import sassie.calculate.sld_mol.gui_mimic_sld_mol as gui_mimic_sld_mol
#import gui_mimic_sld_mol as gui_mimic_sld_mol

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'calculate', 'sld_mol') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Sld_Mol(MockerTestCase):

    '''
    System integration test for sld_mol.py / sassie 1.0

    SLD_MOL is the module that calculates the scattering length density profile
    from a dcd/pdb file

    This method compares an experimentally derived SLD profile with heavy atom
    distribution from a pdb or dcd file containing protein structure(s).
    It performs a fit allowing a normalization factor, z-pos & constant shift.
    The SLD profile is convolved with a Gaussian of user defined width to mimic instrument
    resolution and roughness and outputs a text file with frame number and fit_error.
    
        INPUT:  variable descriptions:

            runname:            project name                                          
            pdbfile:            reference PDB file
            dcdfile:            input filename (DCD or PDB)
            expdatafile:        experimental SLD data file name
            outputfile:         output file name 
            runtype:            0 = average SLD over all structures, 1 = best fit SLD for each individual structure
            bulk_sld:           SLD for bulk solvent
            xon:                scattering type: 0 = neutron, 1 = x-ray
            num_deut_regions:   number of fully deuterated regions in molecule
            deut_low_res:       low residue number(s) for deuterated region(s)
            deut_high_res:      high residue number(s) for deuterated region(s)
            sldfit:             0 = no optimization of z0 and A0, 1 = optimize z0 and A0
            sldoffset:          offset correction to experimental SLD
            dbin:               bin width for z values
            width:              SLD width for Gaussian smoothing
            zfit0:              z0 value (anchoring position) for calculated SLDs initial estimate
            zfitmin:            minimum z value used during optimization
            zfitmax:            maximum z value used during optimization
            zevalmin:           minimum z to evaluate experimental SLD during error calculation
            zevalmax:           maximum z to evaluate experimental SLD during error calculation
            A0:                 fraction surface coverage for calculated SLDs initial estimate
            Amin:               minimum A0 value used during optimization
            Amax:               maximum A0 value used during optimization


        OUTPUT:

            files stored in ./"runname"/sld_mol/ directory:

            outputfile:         output file containing z0, A0 and fit-error for each calculated SLD
            average_sld.txt:    file containing average SLD
            sldfile*.txt:       files containing individual SLDs for each frame (runtype = 1 only)
            bestworstfile.txt:  file containing the filenames and error values for the best- and worst-fitting structures


    Use cases:

    1.  runtype = 0; sldfit = 0
        a. input file is DCD
        b. input file is PDB

    2.  runtype = 0; sldfit = 1
        a. input file is DCD
        b. input file is PDB

    3.  runtype = 1; sldfit = 0
        a. input file is DCD
        b. input file is PDB

    4.  runtype = 1; sldfit = 1     ##not recommended but allowed by software so test for at least one case
        a. input file is DCD
        b. input file is PDB

    Selection options (apply to all cases above):
        a.  sldoffset
        b.  xon
        c.  num_deut_regions
            i.  single deuterated region
            ii. multiple deuterated regions


    Inputs tested:

    runname:            string      project name                                          
    path:               string      input file path
    pdbfile:            string      input PDB file
    dcdfile:            string      input DCD or PDB file
    expdatafile:        string      experimental SLD file
    outputfile:         string      output file containing z0, A0 and fit-error for each SLD
    runtype:            integer     0 = average sld over all structures, 1 = best fit sld for each individual structure
    bulk_sld:           float       SLD for bulk solvent
    xon:                integer     0 = neutron, 1 = x-ray
    num_deut_regions:   integer     number of deuterated regions in molecule
    deut_low_res:       int_array   low residue number(s) for deuterated region(s)
    deut_high_res:      int_array   high residue number(s) for deuterated region(s)
    sldfit:             integer     0 = no optimization of z0 and A0, 1 = optimize z0 and A0
    sldoffset:          float       offset to experimental SLD
    dbin:               float       bin size for z values
    width:              float       bin width for Gaussian smoothing
    zfit0:              float       z0 value for calculated SLDs
    zfitmin:            float       minimum z value used during optimization
    zfitmax:            float       maximum z value used during optimization
    zevalmin:           float       minimum z to evaluate experimental SLD
    zevalmax:           float       maximum z to evaluate experimental SLD
    A0:                 float       fraction surface coverage for calculated SLDs
    Amin:               float       minimum A0 value used during optimization
    Amax:               float       maximum A0 value used during optimization
    plotflag:           integer     flag for plotting data (0: no plot; 1: matplotlib; 2: gnuplot)



    Test tree:

    project name
    input/output path

*****************************
    runtype = 0; sldfit = 0
*****************************
    reference PDB           reference PDB       reference PDB       reference PDB       reference PDB       reference PDB
    input DCD               input PDB           input DCD           input PDB           input DCD           input PDB
    no deut regions         no deut regions     single deut region  single deut region  mult deut regions   multiple deut regions 

                                                    sldoffset
                                                    xon*        *for no deut regions only

  
 
    '''

    module = 'sld_mol'

    def setUp(self):

       gui_mimic_sld_mol.test_variables(self, paths)


    def assert_list_almost_equal(self, a, b, places=5):
        if (len(a) != len(b)):
            raise TypeError
        else:
            for i in range(len(a)):
                if isinstance(a[i], (int, float, numpy.generic)):
                    if (numpy.isnan(a[i]) and numpy.isnan(b[i])):
                        continue
                    self.assertAlmostEqual(a[i], b[i], places)
                else:
                    self.assert_list_almost_equal(a[i], b[i], places)


    def check_lines(self, fname, expected_lines):
        '''
        compares number of expected lines with number of actual lines in a file
        '''
        lines = 0
        with open(fname) as f:
            for line in f:
                lines = lines + 1
#        print 'lines: ', lines
        if (lines == expected_lines):
            return True
        else:
            return False        


    def test_43(self):
        '''
        test runtype=1, sldfit = 1, input DCD, no deut regions
        '''

        self.sldfit = '1'
        self.runtype = '1'
        gui_mimic_sld_mol.run_module(self)

        ''' check for completion by confirming output files have the correct number of lines'''
        avgfile = os.path.join(self.runname, self.module, 'average_sld.txt')
        avgfile_lines = 328
        assert_equals(self.check_lines(avgfile, avgfile_lines), True)
        sld1file = os.path.join(self.runname, self.module, 'sldfile_00001.txt')
        sld1file_lines = 328
        assert_equals(self.check_lines(sld1file, sld1file_lines), True)        
        sld2file = os.path.join(self.runname, self.module, 'sldfile_00002.txt')
        sld2file_lines = 328
        assert_equals(self.check_lines(sld2file, sld2file_lines), True)     
        bwfile = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        bwfile_lines = 4
        assert_equals(self.check_lines(bwfile, bwfile_lines), True)
        resultsfile = os.path.join(self.runname, self.module, 'results_f12.dat')
        resultsfile_lines = 3
        assert_equals(self.check_lines(resultsfile, resultsfile_lines), True)

    def test_44(self):
        '''
        test runtype=1, sldfit = 1, input PDB, no deut regions, sldoffset
        '''

        self.sldfit = '1'
        self.runtype = '1'
        self.sldoffset = '10.0'
        self.zevalmin = '17.5'     
        self.zevalmax = '189.0'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')        
        gui_mimic_sld_mol.run_module(self)

        ''' check for completion by confirming output files have the correct number of lines'''
        avgfile = os.path.join(self.runname, self.module, 'average_sld.txt')
        avgfile_lines = 328
        assert_equals(self.check_lines(avgfile, avgfile_lines), True)
        sld1file = os.path.join(self.runname, self.module, 'sldfile_00001.txt')
        sld1file_lines = 328
        assert_equals(self.check_lines(sld1file, sld1file_lines), True)        
        sld2file = os.path.join(self.runname, self.module, 'sldfile_00002.txt')
        sld2file_lines = 328
        assert_equals(self.check_lines(sld2file, sld2file_lines), True)     
        bwfile = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        bwfile_lines = 4
        assert_equals(self.check_lines(bwfile, bwfile_lines), True)
        resultsfile = os.path.join(self.runname, self.module, 'results_f12.dat')
        resultsfile_lines = 3
        assert_equals(self.check_lines(resultsfile, resultsfile_lines), True)                
        
    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

