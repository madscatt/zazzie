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
import sassie.analyze.density_plot.gui_mimic_density_plot as gui_mimic_density_plot
#import gui_mimic_density_plot as gui_mimic_density_plot

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
    __file__)), '..', '..', 'data', 'interface', 'density_plot') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Density_Plot_Filter(MockerTestCase):

    '''
    System integration test for density_plot_filter.py / sassie 1.0

    Test to see whether density_plot_filter catches improper input.

    Inputs tested:

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


    Use cases tested:

    1.  check runname
        a.  check for invalid characters in runname    
    2.  check input file path permissions 
        a.  no permission error
        b.  permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed
    3.  check input PDB file
        a.  PDB file doesn't exist
        b.  PDB file exists
            i.  PDB file is valid
            ii. PDB file isn't valid
    4.  check input trajectory file
        a.  file doesn't exist
        b.  file exists
            i.  file name doesn't start with "." and must end in ".pdb" or ".dcd"
            ii. file is .pdb (value = 0)
                1.  PDB file is valid
                2.  PDB file isn't valid
            iii.file is .dcd (value = 1)
                1.  DCD file is valid
                2.  DCD file isn't valid
    5.  check if input PDB and trajectory files are compatible
        a.  trajectory file is a DCD file
            i.  input PDB and DCD are compatible
            ii. input PDB and DCD are not compatible
        b.  trajectory file is a PDB file
            i.  input PDB and PDB are compatible
            ii. input PDB and PDB are not compatible
    6.  check box size (xlength, ylength, zlength)
        a.  size is > 0
        b.  size is <= 0
    7.  check grid spacing
        a.  grid spacing is > 0
        b.  grid spacing is <=0
    8.  check equalweights
        a.  equalweights = 0 or 1
        b.  equalweights is not = 0 or 1
    9.  check save_occupancy
        a.  save_occupancy is "Y" or "N"
        b.  save_occupancy is not "Y" or "N"
    10. check weight file
        a.  weight file doesn't exist
        b.  weight file column 3 can only have 0.0 or 1.0
        c.  weights can't all be 0
        d.  weight file column 1 can only have positive integers
        e.  weight file must have structure numbers less than or equal to the number of frames in the trajectory file
        f.  weight file length must be the same as the number of frames in the trajectory file
        g.  weight file can't contain duplicate structure numbers
        h.  unknown error encountered                
    11. check segment variables
        a.  number of segments in segvariables is equal to the value of nsegments
        b.  number of segments in segvariables matches the number of segments in input PDB file
        c.  number of ranges is an integer
        d.  number of ranges is >=1
        e.  low resid is an integer array
        f.  high resid is an integer array
        g.  segment name(s) are found in PDB file
        h.  segment basis is found in PDB file
        i.  number of low/high residues matches the number of ranges
        j.  resid statistics can be obtained for segment
        k.  low resid is present in segment in PDB file
        l.  high resid is present in segment in PDB file
        m.  low resid < high resid               

    '''

    def setUp(self):

        gui_mimic_density_plot.test_variables(self, paths)


    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)

    def test_2(self):
        '''
        test if path exists
        '''
        self.path = os.path.join(module_data_path, 'non_existent_path')
        return_error = gui_mimic_density_plot.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path + '  [code = FalseFalseFalse]',
                          'path does not exist']
        assert_equals(return_error, expected_error)

    def test_3(self):
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

        self.path = os.path.join('./', 'empty_folder')
        return_error = gui_mimic_density_plot.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.path + '  [code = TrueFalseTrue]', 'read permission not allowed']
        assert_equals(return_error, expected_error)

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')        

    def test_4(self):
        '''
        test if directory has write permission
        '''
        self.path = os.path.join(module_data_path, 'no_write_permission')
        return_error = gui_mimic_density_plot.run_module(
            self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path +
                          '  [code = TrueTrueFalse]', 'write permission not allowed']
        assert_equals(return_error, expected_error)

    def test_5(self):
        '''
        test if input PDB file exists
        '''
        self.pdbfile = os.path.join(
            module_data_path, 'does_not_exist.pdb')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          os.path.basename(self.pdbfile) + ', does not exist']                 
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if input PDB file is a valid file
        '''
        self.pdbfile = os.path.join(module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb file, ' +
                          os.path.basename(self.pdbfile) + ', is not a valid pdb file']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if input trajectory file exists
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'does_not_exist.dcd')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' +
                          self.dcdfile + ' does not exist']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if input trajectory file names doesn't start with "." and ends in ".pdb" or ".dcd"
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'test_dcd.fil')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input trajectory filename cannot start with "." and must end with .pdb or .dcd : test_dcd.fil']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if input trajectory file is a valid dcd file
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'not_valid.dcd')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file, ' +
                          os.path.basename(self.dcdfile) + ', is not a valid dcd or pdb file']
        assert_equals(return_error, expected_error)


    def test_10(self):
        '''
        test if input trajectory file is a valid pdb file
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'not_valid.pdb')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input file, ' +
                          os.path.basename(self.dcdfile) + ', is not a valid dcd or pdb file']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if "dcd" trajectory file has the same number of atoms as input PDB
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'non_matching.dcd')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb ' + os.path.basename(self.pdbfile) + ' and dcd file ' +
                          os.path.basename(self.dcdfile) + ' are not compatible']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if "pdb" trajectory file has the same number of atoms as input PDB
        '''
        self.dcdfile = os.path.join(
            module_data_path, 'non_matching.pdb')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['input pdb ' + os.path.basename(self.pdbfile) + ' and dcd file ' +
                          os.path.basename(self.dcdfile) + ' are not compatible']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if xlength <= 0
        '''
        self.xlength = '-1.0'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for SAS type error '''
        expected_error = ['boxlengths need to be > 0, lengths = ' +
                     str(self.xlength) + ',' + str(self.ylength) + ',' + str(self.zlength)]
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test if ylength <= 0
        '''
        self.ylength = '-1.0'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for SAS type error '''
        expected_error = ['boxlengths need to be > 0, lengths = ' +
                     str(self.xlength) + ',' + str(self.ylength) + ',' + str(self.zlength)]
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test if zlength <= 0
        '''
        self.zlength = '0.0'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for SAS type error '''
        expected_error = ['boxlengths need to be > 0, lengths = ' +
                     str(self.xlength) + ',' + str(self.ylength) + ',' + str(self.zlength)]
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if grid spacing <= 0
        '''
        self.gridsp = '0.0'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for SAS type error '''
        expected_error = ['grid spacing needs to be > 0, grid spacing = ' + str(self.gridsp)]
        assert_equals(return_error, expected_error)
        
    def test_17(self):
        '''
        test if equalweights = 0 or 1
        '''
        self.equalweights = '2'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for SAS type error '''
        expected_error = ['equalweights == 0 for "no" and 1 for "yes", equalweights = ' + str(self.equalweights)]
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if save occupancy = "Y" or "N"  (case sensitive)
        '''
        self.save_occupancy = 'y'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for SAS type error '''
        expected_error = ['save occupancy data needs to be either "Y" or "N": you entered: ' + self.save_occupancy]
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test if weight file exists
        '''

        self.equalweights = "0"
        self.weightsfile = 'non_existent_file.txt'
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['file : ' + self.weightsfile + ' does not exist', 'weight file not readable or does not exist']
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        test if length of weight file is equal to the number of frames in the trajectory file 
        '''

        self.equalweights = "0"
        self.weightsfile = os.path.join(module_data_path,'wrong_number_of_weight_values.txt')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['number of lines 19 in weightsfile, ' + os.path.basename(self.weightsfile) + ', does not match the number of frames in input file, ' + os.path.basename(self.dcdfile) + '; number of frames = 20']
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test for unknown error encountered when reading weight file
        '''

        self.equalweights = "0"
        self.weightsfile = os.path.join(
            module_data_path, 'weird_weight_file.txt')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ["unable to open and read your weight file : " + os.path.basename(self.weightsfile)]
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test for all zeros in column 3 of weight file
        '''

        self.equalweights = "0"
        self.weightsfile = os.path.join(
            module_data_path, 'all_zeros_weight_file.txt')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ["all weights in weights file are zero : " + os.path.basename(self.weightsfile)]
        assert_equals(return_error, expected_error)

    def test_23(self):
        '''
        test for numbers other than 0 or 1 in column 3 of weight file
        '''

        self.equalweights = "0"
        self.weightsfile = os.path.join(
            module_data_path, 'bad_column_three_weight_file.txt')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['weight file column 3 can only have 0 or 1 : 5.000000 was found']
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test for non-positive integer in weight file
        '''

        self.equalweights = "0"
        self.weightsfile = os.path.join(
            module_data_path, 'negative_number_weight_file.txt')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['weight file column 1 can only have positive integers : "-7" was found in the weight file, '
             + os.path.basename(self.weightsfile)]
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test if there are structure numbers greater than number of frames in trajectory file in weight file
        '''

        self.equalweights = "0"
        self.weightsfile = os.path.join(
            module_data_path, 'number_too_large_weight_file.txt')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = ['there are 20 frames in input file ' + os.path.basename(self.dcdfile) +
                ': "25" was found in weight file, ' + os.path.basename(self.weightsfile)]
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test if the number of segments in segvariables is equal to the value of nsegments
        '''

        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CA', u'VN1']] 
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)
#        print return_error

        ''' check for file error '''
        expected_error = [['The number of segments in segvariables is not equal to the value of nsegments!']]
#        print expected_error
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test if the number of segments in segvariables matches the number of segments in the input PDB file
        '''

        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The number of segments in segvariables is not equal to that from the pdb file!']]
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test if the number of ranges is an integer
        '''

        self.segvariables = [[u'1.5', '6', '123', u'CA', u'GAG']]
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'The number of ranges "1.5" for segment number 0 in the segment input fields should be an integer type!']]
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test if the number of ranges is >=1 (one segment)
        '''

        self.segvariables = [[u'0', '6', '123', u'CA', u'GAG']]
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'The number of ranges "0" for segment number 0 in the segment input fields should be equal/greater than 1!']]
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test if number of ranges >= (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'0', '1,40', '39,130', u'CA', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'The number of ranges "0" for segment number 1 in the segment input fields should be equal/greater than 1!']]
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test if low resid is an integer array (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'3', '1,31.5,47', '30,46,80', u'P', u'TR2A']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The low resid "1,31.5,47" for segment number 0 in the segment input fields should be an integer array!']]
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test if low resid is an integer array (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1.5,40', '39,130', u'CA', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The low resid "1.5,40" for segment number 1 in the segment input fields should be an integer array!']]
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test if high resid is an integer array (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'3', '1,31,47', '30,46,80.0', u'P', u'TR2A']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The high resid "30,46,80.0" for segment number 0 in the segment input fields should be an integer array!']]
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test if high resid is an integer array (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130.2', u'CA', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['The high resid "39,130.2" for segment number 1 in the segment input fields should be an integer array!']]
        assert_equals(return_error, expected_error)

    def test_35(self):
        '''
        test if segment name is in PDB file (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'NONE']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'The segment name "NONE" is not found in the pdb file!']]
        assert_equals(return_error, expected_error)

    def test_36(self):
        '''
        test if segment name is found in PDB file (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CA', u'VM1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'The segment name "VM1" is not found in the pdb file!']]
        assert_equals(return_error, expected_error)

    def test_37(self):
        '''
        test if basis is in PDB file (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'3', '1,31,47', '30,46,80', u'B', u'TR2A']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'The segment basis "B" is not found in the pdb file!']]
        assert_equals(return_error, expected_error)

    def test_38(self):
        '''
        test if basis is found in PDB file (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CC', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'The segment basis "CC" is not found in the pdb file!']]
        assert_equals(return_error, expected_error)

    def test_39(self):
        '''
        test if number of low/high residues matches the number of ranges (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'2', '1,31,47', '30,46,80', u'P', u'TR2A']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['the number of low/high residue input does not match the number of ranges!',
            'arlow = [[1, 31, 47]] : arhigh = [[30, 46, 80]]','len(arlow) = 1 : len(arhigh) = 1',
            'numranges = 2']]
        assert_equals(return_error, expected_error)

    def test_40(self):
        '''
        test if number of low/high residues matches the number of ranges (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'3', '1,40', '39,130', u'CA', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['the number of low/high residue input does not match the number of ranges!', 
             'arlow = [[1], [1, 40]] : arhigh = [[379], [39, 130]]', 'len(arlow) = 2 : len(arhigh) = 2', 
             'numranges = 3']]
        assert_equals(return_error, expected_error)

    def test_41(self):
        '''
        test if number of low resid is in PDB file (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'3', '1,31,81', '30,46,85', u'P', u'TR2A']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['input pdb file, trunc2a_min.pdb does not have low residue amino acid, 81, range = 1 : 80']]
        assert_equals(return_error, expected_error)

    def test_42(self):
        '''
        test if low resid is in PDB file (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,131', '39,132', u'CA', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['input pdb file, pai_vn_start.pdb does not have low residue amino acid, 131, range = 1 : 130']]
        assert_equals(return_error, expected_error)

    def test_43(self):
        '''
        test if number of high resid is in PDB file (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'3', '1,31,47', '30,46,85', u'P', u'TR2A']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['input pdb file, trunc2a_min.pdb does not have high residue amino acid, 85, range = 1 : 80']]
        assert_equals(return_error, expected_error)

    def test_44(self):
        '''
        test if high resid is in PDB file (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,132', u'CA', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [['input pdb file, pai_vn_start.pdb does not have high residue amino acid, 132, range = 1 : 130']]
        assert_equals(return_error, expected_error)

    def test_45(self):
        '''
        test if low resid > high resid (one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'trunc2a_min.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'trunc2a_20_frames.dcd')
        self.segvariables = [[u'3', '1,46,47', '30,31,80', u'P', u'TR2A']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'the low-resid: 46, is greater than/equal to the high-resid: 31 in segname: TR2A']]
        assert_equals(return_error, expected_error)

    def test_46(self):
        '''
        test if low resid > high resid (more than one segment)
        '''
        self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
        self.dcdfile = os.path.join(dcd_data_path,'pai_vn_20_frames.dcd')
        self.nsegments = '2'
        self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '39,40', '1,130', u'CA', u'VN1']]   
        return_error = gui_mimic_density_plot.run_module(
            self, test_filter=True)

        ''' check for file error '''
        expected_error = [[u'the low-resid: 39, is greater than/equal to the high-resid: 1 in segname: VN1']]
        assert_equals(return_error, expected_error)
 

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)
            
if __name__ == '__main__':
    main()
