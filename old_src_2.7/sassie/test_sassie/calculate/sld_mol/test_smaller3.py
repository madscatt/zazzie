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

    def test_29(self):
        '''
        test runtype=1, sldfit = 0, input DCD, no deut regions
        '''

        self.precision = 5
        self.runtype = '1'
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_30(self):
        '''
        test runtype=1, sldfit = 0, input DCD, one deut region
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '1'
        self.z0 = '-1.1'
        self.A0 = '0.12'
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_31(self):
        '''
        test runtype=1, sldfit = 0, input DCD, mult deut regions
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '2'
        self.lowres = '1,150'
        self.highres = '145,200'
        self.z0 = '-1.1'
        self.A0 = '0.12'
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_32(self):
        '''
        test runtype=1, sldfit = 0, input DCD, no deut regions, sldoffset
        '''

        self.precision = 5
        self.runtype = '1'
        self.z0 = '1.7'
        self.zevalmin = '17.5'     
        self.zevalmax = '189.0'    
        self.sldoffset = '10.0'          
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_sldoffset', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_sldoffset', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_sldoffset', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_sldoffset', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_sldoffset', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_33(self):
        '''
        test runtype=1, sldfit = 0, input DCD, one deut region, sldoffset
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '1'
        self.z0 = '3.8'
        self.A0 = '0.12'
        self.zevalmin = '17.5'     
        self.zevalmax = '189.0'    
        self.sldoffset = '10.0'                 
        gui_mimic_sld_mol.run_module(self)

        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one_sldoffset', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one_sldoffset', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one_sldoffset', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one_sldoffset', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_one_sldoffset', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)
        
    def test_34(self):
        '''
        test runtype=1, sldfit = 0, input DCD, mult deut regions, sldoffset
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '2'
        self.lowres = '1,150'
        self.highres = '145,200'
        self.z0 = '3.8'
        self.A0 = '0.12'
        self.zevalmin = '17.5'     
        self.zevalmax = '189.0'    
        self.sldoffset = '10.0'         
        gui_mimic_sld_mol.run_module(self)


        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult_sldoffset', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult_sldoffset', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult_sldoffset', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult_sldoffset', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_mult_sldoffset', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_35(self):
        '''
        test runtype=1, sldfit = 0, input PDB, no deut regions
        '''

        self.precision = 5
        self.runtype = '1'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')
        gui_mimic_sld_mol.run_module(self)


        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_36(self):
        '''
        test runtype=1, sldfit = 0, input PDB, one deut region
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '1'
        self.z0 = '-1.1'
        self.A0 = '0.12'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')
        gui_mimic_sld_mol.run_module(self)

        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)
        
    def test_37(self):
        '''
        test runtype=1, sldfit = 0, input PDB, mult deut regions
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '2'
        self.lowres = '1,150'
        self.highres = '145,200'
        self.z0 = '-1.1'
        self.A0 = '0.12'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)
        
    def test_38(self):
        '''
        test runtype=1, sldfit = 0, input PDB, no deut regions, sldoffset
        '''

        self.precision = 5
        self.runtype = '1'
        self.z0 = '1.7'
        self.zevalmin = '17.5'     
        self.zevalmax = '189.0'    
        self.sldoffset = '10.0'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')          
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_sldoffset', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_sldoffset', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_sldoffset', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_sldoffset', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_sldoffset', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_39(self):
        '''
        test runtype=1, sldfit = 0, input PDB, one deut region, sldoffset
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '1'
        self.z0 = '3.8'
        self.A0 = '0.12'
        self.zevalmin = '17.5'     
        self.zevalmax = '189.0'    
        self.sldoffset = '10.0'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')                         
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one_sldoffset', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one_sldoffset', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one_sldoffset', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one_sldoffset', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_one_sldoffset', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

        
    def test_40(self):
        '''
        test runtype=1, sldfit = 0, input PDB, mult deut regions, sldoffset
        '''

        self.precision = 5
        self.runtype = '1'
        self.numdregions = '2'
        self.lowres = '1,150'
        self.highres = '145,200'
        self.z0 = '3.8'
        self.A0 = '0.12'
        self.zevalmin = '17.5'     
        self.zevalmax = '189.0'    
        self.sldoffset = '10.0'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')                 
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult_sldoffset', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult_sldoffset', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult_sldoffset', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult_sldoffset', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)
        ''' compare bestworstfile.txt '''
        bestworst_file = os.path.join(self.runname, self.module, 'bestworstfile.txt')
        corr_bestworst_file = os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_mult_sldoffset', 'bestworstfile.txt')
        assert_equals(filecmp.cmp(bestworst_file, corr_bestworst_file), True)

    def test_41(self):
        '''
        test runtype=1, sldfit = 0, input DCD, no deut regions, xon
        '''

        self.precision = 5
        self.runtype = '1'
        self.xon = '1'
        gui_mimic_sld_mol.run_module(self)
        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_xon', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_xon', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_xon', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_dcd_none_xon', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)

    def test_42(self):
        '''
        test runtype=1, sldfit = 0, input PDB, no deut regions, xon
        '''

        self.precision = 5
        self.runtype = '1'
        self.xon = '1'
        self.dcdfile = os.path.join(pdb_data_path, 'f12.pdb')                         
        gui_mimic_sld_mol.run_module(self)

        
        ''' confirm values in output files are correct to within 5 decimal places '''
        outfile = open(os.path.join(self.runname, self.module, self.outputfile), 'r').readlines()
        outz = []
        outa = []
        outerr = []
        for i in range(len(outfile)):
            lin = string.split(outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                outz.append(float(lin[1]))
                outa.append(float(lin[2]))
                outerr.append(float(lin[3]))
        correct_outfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_xon', self.outputfile), 'r').readlines()
        corroutz = []
        corrouta = []
        corrouterr = []
        for i in range(len(correct_outfile)):
            lin = string.split(correct_outfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
                corroutz.append(float(lin[1]))
                corrouta.append(float(lin[2]))
                corrouterr.append(float(lin[3]))            
        self.assert_list_almost_equal(corroutz, outz, self.precision)
        self.assert_list_almost_equal(corrouta, outa, self.precision)
        self.assert_list_almost_equal(corrouterr, outerr, 3)        
        avgfile = open(os.path.join(self.runname, self.module, 'average_sld.txt'), 'r').readlines()
        avgz = []
        avgsld = []
        for i in xrange(len(avgfile)):
            lin = string.split(avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                avgz.append(float(lin[0]))
                avgsld.append(float(lin[1]))        
        corr_avgfile = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_xon', 'average_sld.txt'), 'r').readlines()
        corravgz = []
        corravgsld = []
        for i in xrange(len(corr_avgfile)):
            lin = string.split(corr_avgfile[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corravgz.append(float(lin[0]))
                corravgsld.append(float(lin[1]))        
        self.assert_list_almost_equal(corravgz, avgz, self.precision)
        self.assert_list_almost_equal(corravgsld, avgsld, self.precision)
        sldfile1 = open(os.path.join(self.runname, self.module, 'sldfile_00001.txt'), 'r').readlines()
        sld1z = []
        sld1 = []
        for i in xrange(len(sldfile1)):
            lin = string.split(sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld1z.append(float(lin[0]))
                sld1.append(float(lin[1]))        
        corr_sldfile1 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_xon', 'sldfile_00001.txt'), 'r').readlines()
        corrsld1z = []
        corrsld1 = []
        for i in xrange(len(corr_sldfile1)):
            lin = string.split(corr_sldfile1[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld1z.append(float(lin[0]))
                corrsld1.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld1z, sld1z, self.precision)
        self.assert_list_almost_equal(corrsld1, sld1, self.precision)
        sldfile2 = open(os.path.join(self.runname, self.module, 'sldfile_00002.txt'), 'r').readlines()
        sld2z = []
        sld2 = []
        for i in xrange(len(sldfile2)):
            lin = string.split(sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                sld2z.append(float(lin[0]))
                sld2.append(float(lin[1]))        
        corr_sldfile2 = open(os.path.join(
            module_data_path, self.runname, self.module, '10_pdb_none_xon', 'sldfile_00002.txt'), 'r').readlines()
        corrsld2z = []
        corrsld2 = []
        for i in xrange(len(corr_sldfile2)):
            lin = string.split(corr_sldfile2[i])    
            if(lin[0][0] != "#" and len(lin) >= 2):
                corrsld2z.append(float(lin[0]))
                corrsld2.append(float(lin[1]))        
        self.assert_list_almost_equal(corrsld2z, sld2z, self.precision)
        self.assert_list_almost_equal(corrsld2, sld2, self.precision)                
        
    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

