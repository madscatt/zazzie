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
import sassie.analyze.density_plot.gui_mimic_density_plot as gui_mimic_density_plot
#import gui_mimic_density_plot as gui_mimic_density_plot

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'analyze', 'density_plot') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}

class Test_Density_Plot(MockerTestCase):

    '''
    System integration test for density_plot.py / sassie 1.0

    DENSITY_PLOT is the function to read in variables from GUI input and compare
    generate three-dimensional volumetric data files using the GAUSSIAN file
    format. 

    INPUT:  variable descriptions:

            runname:        project name
            path:           path name for input files
            dcdfile:        input trajectory filename (pdb or dcd)
            pdbfile:        reference pdb name
            xlength:        x boxlength
            ylength:        y boxlength
            zlength:        z boxlength
            gridsp:         grid spacing (angstroms)            
            nsegments:	    number of segments
            segvariables:   number of flexible regions, high and low residues, basis string, segment name
            save_occupancy: save the unweighted raw cube data ('Y' or 'N')
            equalweights:   use equalweights (1=yes) or weights from file (0=no)
            weightsfile:    filename containing weights per structure


    OUTPUT:
            ofile:              	output filename prefix

            files stored in ~/runname/filter directory:

            *_complete.cube:	Gaussian volumetric cube file of all basis atoms
            *_region_X.cube:	Gaussian volumetric cube file of basis atoms in region X


    The output density will be normalized as follows against the maximum density value from the composite map including all atoms:
    rho(norm)[i][j][k] = rho[i][j][k]*100.0/max(rho)
    where i=1,2,...,Nsegments, j=1,2,...,Nregions, and k=1,2,...,Ngridpoints

    Use cases:

    1.  Input trajectory file
        a.  input file is a PDB file
        b.  input file is a DCD file
    2.  Single Segment
        a.  single flexible region
        b.  multiple flexible regions
    3.  Multiple Segments
        a.  single flexible region
        b.  multiple flexible regions

    Selection option (apply to all cases above):
        a.  weight file
        b.  save occupancy   

    Inputs tested:

    runname:                string      project name
    path:                   string      path name for input files    
    pdbfile:                string      input reference pdb file
    dcdfile:                string      input trajectory file (pdb or dcd)
    ofile:                  string      output file name prefix                         
    xlength:                float       x boxlength
    ylength                 float       y boxlength
    zlength:                float       z boxlength             
    gridsp:                 float       grid spacing (angstroms)             
    nsegments:              integer     number of segments   
    save_occupancy          string      save the unweighted raw cube data (only 'Y' or 'N' allowed, case-sensitive)
    equalweights            integer     use equalweights (1=yes) or weights from file (0=no)
    weightsfile:            string      filename containing weights per structure

    segvariables:           variables for each segment (listed below)  
    number of flexible regions:         string              number of flexible regions per segment
    high and low residies:              integer array       high and low residue numbers for each flexible region
    basis string:                       string              basis atom for each segment
    segment name:                       string              segment name for each segment



    Test tree:

    project name
    input/output path

**************************
    single segment 
**************************
    reference PDB           reference PDB           reference PDB           reference PDB       *reference PDB file
    input PDB               input DCD               input PDB               input DCD           *trajectory file
    single region           single region           multiple regions        multiple regions    *number of flexible regions
                                           save occupancy
                                           no weight file
    
    reference PDB           reference PDB           reference PDB           reference PDB       *reference PDB file
    input PDB               input DCD               input PDB               input DCD           *trajectory file
    single region           single region           multiple regions        multiple regions    *number of flexible regions
                                           no save occupancy
                                           no weight file                                     

    reference PDB           reference PDB           reference PDB           reference PDB       *reference PDB file
    input PDB               input DCD               input PDB               input DCD           *trajectory file
    single region           single region           multiple regions        multiple regions    *number of flexible regions
                                           save occupancy
                                           weight file
    
    reference PDB           reference PDB           reference PDB           reference PDB       *reference PDB file
    input PDB               input DCD               input PDB               input DCD           *trajectory file
    single region           single region           multiple regions        multiple regions    *number of flexible regions
                                           no save occupancy
                                           weight file                                     


**********************************
    multiple segments
**********************************
Since single region, multiple regions, save occupancy and weight file options were tested above,
the following tests were chosen for multiple segments:

    reference PDB           reference PDB       *reference PDB file
    input PDB               input DCD           *trajectory file
    multiple regions        multiple regions    *number of flexible regions
                       no save occupancy
                       no weight file


    reference PDB           reference PDB       *reference PDB file
    input PDB               input DCD           *trajectory file
    multiple regions        multiple regions    *number of flexible regions
                       save occupancy
                       weight file                                     
 
    '''

    module = 'density_plot'

    def setUp(self):

       gui_mimic_density_plot.test_variables(self, paths)


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

    def check_dir_trees_equal(self, dir1, dir2):
        '''
        compares directories recursively as well as files within them,
        ignoring *.sassie_json and *.sassie_log files
        these files must be ignored since they have date stamps in their names
        '''

        dirs_cmp = filecmp.dircmp(dir1, dir2)
#        if len(dirs_cmp.left_only)>0 or len(dirs_cmp.right_only)>0 or \
#            len(dirs_cmp.funny_files)>0:
        if len(dirs_cmp.right_only)>0 or len(dirs_cmp.funny_files)>0:
            return False
        (_, mismatch, errors) =  filecmp.cmpfiles(
            dir1, dir2, dirs_cmp.common_files)
        if len(mismatch)>0 or len(errors)>0:
            return False
        for common_dir in dirs_cmp.common_dirs:
            new_dir1 = os.path.join(dir1, common_dir)
            new_dir2 = os.path.join(dir2, common_dir)
            if not self.check_dir_trees_equal1(new_dir1, new_dir2):
                return False
        return True


    def test_1(self):
        '''
        test PDB input, single region, save occupancy, no weight file
        '''

        self.save_occupancy = 'Y'
        self.dcdfile = os.path.join(pdb_data_path,'hiv1_gag_20_frames.pdb')   

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_occ_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                       
    def test_2(self):
        '''
        test DCD input, single region, save occupancy, no weight file
        '''

        self.save_occupancy = 'Y'  

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_occ_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_3(self):
        '''
        test PDB input, multiple regions, save occupancy, no weight file
        '''

        self.save_occupancy = 'Y'  
        self.dcdfile = os.path.join(pdb_data_path,'trunc2a_20_frames.pdb')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_mult_occ_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                       
    def test_4(self):
        '''
        test dcd input, multiple regions, save occupancy, no weight file
        '''

        self.save_occupancy = 'Y'  
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_mult_occ_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_5(self):
        '''
        test PDB input, single region, no save occupancy, no weight file
        '''

        self.dcdfile = os.path.join(pdb_data_path,'hiv1_gag_20_frames.pdb')   

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                       
    def test_6(self):
        '''
        test DCD input, single region, no save occupancy, no weight file
        '''
 
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_7(self):
        '''
        test PDB input, multiple regions, no save occupancy, no weight file
        '''

        self.dcdfile = os.path.join(pdb_data_path,'trunc2a_20_frames.pdb')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_mult_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_8(self):
        '''
        test dcd input, multiple regions, no save occupancy, no weight file
        '''
 
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_mult_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_9(self):
        '''
        test PDB input, single region, save occupancy, weight file
        '''

        self.save_occupancy = 'Y'
        self.equalweights = '0'
        self.dcdfile = os.path.join(pdb_data_path,'hiv1_gag_20_frames.pdb')
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt')   

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_occ_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                       
    def test_10(self):
        '''
        test DCD input, single region, save occupancy, weight file
        '''

        self.save_occupancy = 'Y' 
        self.equalweights = '0'         
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt') 
        
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_occ_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_11(self):
        '''
        test PDB input, multiple regions, save occupancy, weight file
        '''

        self.save_occupancy = 'Y'
        self.equalweights = '0'
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt')
        self.dcdfile = os.path.join(pdb_data_path,'trunc2a_20_frames.pdb')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_mult_occ_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_12(self):
        '''
        test DCD input, multiple regions, save occupancy, weight file
        '''

        self.save_occupancy = 'Y' 
        self.equalweights = '0'         
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt') 
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 
                
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_mult_occ_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_13(self):
        '''
        test PDB input, single region, no save occupancy, weight file
        '''

        self.equalweights = '0'
        self.dcdfile = os.path.join(pdb_data_path,'hiv1_gag_20_frames.pdb')
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt')   

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)
                       
    def test_14(self):
        '''
        test DCD input, single region, no save occupancy, weight file
        '''

        self.equalweights = '0'         
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt') 
        
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_15(self):
        '''
        test PDB input, multiple regions, no save occupancy, weight file
        '''

        self.equalweights = '0'
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt')
        self.dcdfile = os.path.join(pdb_data_path,'trunc2a_20_frames.pdb')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_mult_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_16(self):
        '''
        test DCD input, multiple regions, no save occupancy, weight file
        '''

        self.equalweights = '0'         
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt') 
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')        
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']] 
                
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_mult_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_17(self):
        '''
        test PDB input, multiple segments, multiple regions, no save occupancy, no weight file
        '''

        self.nsegments = '2'
        self.dcdfile = os.path.join(pdb_data_path,'pai_vn_20_frames.pdb')
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')        
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CA', u'VN1']]

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_mseg_mult_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_18(self):
        '''
        test DCD input, multiple segments, multiple regions, no save occupancy, no weight file
        '''

        self.nsegments = '2'
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')        
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CA', u'VN1']]
                
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_mseg_mult_nowt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

                                                             
    def test_19(self):
        '''
        test PDB input, multiple segments, multiple regions, save occupancy, weight file
        '''

        self.equalweights = '0'
        self.save_occupancy = 'Y'
        self.nsegments = '2'
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt')
        self.dcdfile = os.path.join(pdb_data_path,'pai_vn_20_frames.pdb')
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')        
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CA', u'VN1']]

        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'pdb_mseg_mult_occ_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)

    def test_20(self):
        '''
        test DCD input, multiple segments, multiple regions, save occupancy, weight file
        '''

        self.equalweights = '0'         
        self.save_occupancy = 'Y'
        self.nsegments = '2'
        self.weightsfile = os.path.join(other_data_path,'weights_file_density_plot.txt')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')        
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CA', u'VN1']]
                
        gui_mimic_density_plot.run_module(self)

        ''' confirm output dcd file is correct '''
        outdirectory = os.path.join(self.runname, self.module)
        print 'outdirectory: ', outdirectory
        correct_outdirectory = os.path.join(
            module_data_path, self.runname, self.module, 'dcd_mseg_mult_occ_wt')
        print 'correct outdirectory: ', correct_outdirectory
        assert_equals(self.check_dir_trees_equal(outdirectory,correct_outdirectory), True)


                       
    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)


if __name__=='__main__':
    main()

