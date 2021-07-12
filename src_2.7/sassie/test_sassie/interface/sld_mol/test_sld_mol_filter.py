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
import sassie.calculate.sld_mol.gui_mimic_sld_mol as gui_mimic_sld_mol
#import gui_mimic_sld_mol as gui_mimic_sld_mol

import filecmp
from unittest import main
from nose.tools import assert_equals
from mocker import Mocker, MockerTestCase

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'interface', 'sld_mol') + os.path.sep

paths = {'pdb_data_path' : pdb_data_path, 'dcd_data_path' : dcd_data_path, 'other_data_path' : other_data_path, 'module_data_path' : module_data_path}


class Test_Sld_Mol_Filter(MockerTestCase):

    '''
    System integration test for sld_mol.py / sassie 1.0

    SLD_MOL is the module that calculates the scattering length density profile
    from a dcd/pdb file

    This method compares an experimentally derived SLD profile with heavy atom
    distribution from a pdb or dcd file containing protein structure(s).
    It performs a fit allowing a normalization factor, z-pos & constant shift.
    The SLD profile is convolved with a Gaussian of user defined width to mimic instrument
    resolution and roughness and outputs a text file with frame number and fit_error.

    Inputs tested:

    runname:            string      project name                                          
    path:               string      input file path
    pdbfile:            string      input PDB file
    dcdfile:            string      input DCD file
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

    Use cases tested:

    1.  check if runname has incorrect character
    2.  check input file path permissions 
        a.  no permission error
        b. permission error
            i.   path doesn't not exist
            ii.  read permission not allowed
            iii. write permission not allowed 
    3.  check if runtype is 0 or 1         
    4.  check if xon is 0 or 1
    5.  check if plotflag is 0, 1 or 2
    6.  check if number of deuterated regions is >= 0
    7.  check if number of low deuterated values matches the number of deuterated regions
    8.  check if number of high deuterated values matches the number of deuterated regions
    9.  check if experimental SLD file exists
    10. check for repeated z-values in experimental SLD
    11. check if zevalmin is less than experimental SLD + sldoffset
    12. check if zevalmax is greater than experimental SLD + sldoffset
    13. check if input dcd file exists
    14. check if reference pdb file exists
    15. check if pdb and dcd files are compatible
    16. check if dcd trajectory file is a valid dcd file
    17. check if pdb trajectory file is compatible with reference pdb file
    18. check if dbin is between 0 and 1
    19. check if width is between 0 and 5
    20. check if sldfit is between 0 and 1
    21. check if A0, Amin, Amax > 0 and A0, Amin, Amax < 1
    22. check if Amin >= Amax
    23. check if Amin > A0 or Amax < A0
    24. check if zevalmin < zfitmin
    25. check if zfitmin > zfitmax
    26. check if zfitmin > zfit0 or zfitmax < zfit0
    27. check if zevalmin >= zevalmax
        
    '''

    def setUp(self):

        gui_mimic_sld_mol.test_variables(self, paths)

    def test_1(self):
        '''
        test if runname has incorrect character
        '''
        self.runname = 'run_&'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file or path : run_& has incorrect character : &']
        assert_equals(return_error, expected_error)


    def test_2(self):
        '''
        test if path exists
        '''
        self.path = os.path.join(module_data_path,'non_existent_path')
        return_error = gui_mimic_sld_mol.run_module(self, file_check=True)

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
        print os.access('empty_folder', os.R_OK)
        ''' make the directory un-readable'''
        os.system('chmod a-r empty_folder')
        ''' see if you can read the directory '''
        print os.access('empty_folder', os.R_OK)

        self.path= os.path.join('./','empty_folder')
        return_error = gui_mimic_sld_mol.run_module(self, file_check=True)

        ''' check for path error '''
        expected_error = ['permission error in input file path ' + self.path + '  [code = TrueFalseTrue]', 'read permission not allowed']
        assert_equals(return_error, expected_error)  

        ''' make the directory readable'''
        os.system('chmod a+r empty_folder')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder')

    def test_4(self):
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

        self.path = os.path.join('./', 'empty_folder1')
        return_error = gui_mimic_sld_mol.run_module(
            self, file_check=True)
#        print 'return_error: ', return_error

        ''' check for path error '''
        expected_error = ['permission error in input file path ' +
                          self.path + '  [code = TrueTrueFalse]', 'write permission not allowed']
#        print 'expected_error: ', expected_error
        assert_equals(return_error, expected_error)

        ''' make the directory writeable'''
        os.system('chmod a+w empty_folder1')
        ''' remove the directory '''
        os.system('rm -Rf empty_folder1')        

    def test_5(self):
        '''
        test if runtype is 0 or 1
        '''
        self.runtype = '2'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['run type entered needs to be either 0 or 1 : 2']
        assert_equals(return_error, expected_error)

    def test_6(self):
        '''
        test if xon is 0 or 1
        '''
        self.xon = '2'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['scattering type (xon) entered needs to be either 0 or 1 : 2']
        assert_equals(return_error, expected_error)

    def test_7(self):
        '''
        test if plotflag is 0, 1 or 2
        '''
        self.plotflag = '3'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['plotflag needs to be 0, 1, or 2 : 3']
        assert_equals(return_error, expected_error)

    def test_8(self):
        '''
        test if number of deuterated regions is >= 0
        '''
        self.numdregions = '-1'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['number of deuterated regions needs to be >= 0 : -1']
        assert_equals(return_error, expected_error)

    def test_9(self):
        '''
        test if number of number of low deuterated values matches the number of deuterated regions
        '''
        self.numdregions = '1'
        self.lowres = '1,10'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['number of low deuterated values does not match the number of regions: len(deut_low_res) = 2 num_deut_regions = 1']
        assert_equals(return_error, expected_error)

    def test_10(self):
        '''
        test if number of high deuterated values matches the number of deuterated regions
        '''
        self.numdregions = '1'
        self.highres = '20,50'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['number of high deuterated values does not match the number of regions: len(deut_high_res) = 2 num_deut_regions = 1']
        assert_equals(return_error, expected_error)

    def test_11(self):
        '''
        test if expdatafile exists
        '''

        self.expdatafile = 'non_existent.dat'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file : non_existent.dat does not exist']
        assert_equals(return_error, expected_error)

    def test_12(self):
        '''
        test if there are repeated z-values in experimental SLD
        '''

        self.expdatafile = os.path.join(module_data_path,'repeated_z_value.dat')
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['repeated z-value in experimental data: 11.871\n']
        assert_equals(return_error, expected_error)

    def test_13(self):
        '''
        test if zevalmin < min (z_values) + sldoffset
        '''

        self.zevalmin = '2.5'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['minimum evaluation value is less than experimental data + sldoffset: minz = 6.431 : zevalmin = 2.5']
        assert_equals(return_error, expected_error)

    def test_14(self):
        '''
        test if zevalmax > max (z_values) + sldoffset
        '''

        self.zevalmax = '190.0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['maximum evaluation value is greater than experimental data + sldoffset: maxz = 179.881 : zevalmax = 190.0']
        assert_equals(return_error, expected_error)

    def test_15(self):
        '''
        test if dcdfile exists
        '''

        self.dcdfile = 'non_existent.dcd'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file : non_existent.dcd does not exist']
        assert_equals(return_error, expected_error)

    def test_16(self):
        '''
        test if pdbfile exists
        '''

        self.pdbfile = 'non_existent.dcd'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['file : non_existent.dcd does not exist']
        assert_equals(return_error, expected_error)
       
    def test_17(self):
        '''
        test if pdb and dcd files are compatible
        '''

        self.dcdfile = os.path.join(module_data_path,'non_matching.dcd')
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['pdbfile '+self.pdbfile+' and dcdfile '+self.dcdfile+' are not compatible (different number of atoms)']
        assert_equals(return_error, expected_error)

    def test_18(self):
        '''
        test if dcd file is a valid dcd file
        '''

        self.dcdfile = os.path.join(module_data_path,'not_valid.dcd')
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['input file '+self.dcdfile+' is not a valid pdb or dcd file or it does not exist']
        assert_equals(return_error, expected_error)

    def test_19(self):
        '''
        test if pdb trajectory file and reference pdb files are compatible
        '''

        self.dcdfile = os.path.join(module_data_path,'non_matching.pdb')
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['reference pdb file '+self.pdbfile+' and input pdb file '+self.dcdfile+' are not compatible']
        assert_equals(return_error, expected_error)

    def test_20(self):
        '''
        test if dbin is between 0 and 1
        '''

        self.dbin = '1.1'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['bin width needs to be greater than zero and less than 1.0 angstroms: 1.1']
        assert_equals(return_error, expected_error)

    def test_21(self):
        '''
        test if width is between 0 and 5
        '''

        self.width = '-0.5'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['smoothing width needs to be greater than zero and less than 5.0 angstroms: -0.5']
        assert_equals(return_error, expected_error)

    def test_22(self):
        '''
        test if sldfit is between 0 and 1
        '''

        self.sldfit = '2'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['fit SLD needs to be either 0 or 1 :2']
        assert_equals(return_error, expected_error)

    def test_23(self):
        '''
        test if A0 > 0
        '''

        self.A0 = '0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage has to be a positive value less than or equal to 1.0: A0: 0.0 Amin :0.1 Amax: 0.5']
        assert_equals(return_error, expected_error)

    def test_24(self):
        '''
        test if Amin > 0
        '''

        self.Amin = '0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage has to be a positive value less than or equal to 1.0: A0: 0.27 Amin :0.0 Amax: 0.5']
        assert_equals(return_error, expected_error)

    def test_25(self):
        '''
        test if Amax > 0
        '''

        self.Amax = '0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage has to be a positive value less than or equal to 1.0: A0: 0.27 Amin :0.1 Amax: 0.0']
        assert_equals(return_error, expected_error)

    def test_26(self):
        '''
        test if A0 > 1
        '''

        self.A0 = '1.1'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage has to be a positive value less than or equal to 1.0: A0: 1.1 Amin :0.1 Amax: 0.5']
        assert_equals(return_error, expected_error)

    def test_27(self):
        '''
        test if Amin > 1
        '''

        self.Amin = '2.0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage has to be a positive value less than or equal to 1.0: A0: 0.27 Amin :2.0 Amax: 0.5']
        assert_equals(return_error, expected_error)

    def test_28(self):
        '''
        test if Amax > 1
        '''

        self.Amax = '1.5'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage has to be a positive value less than or equal to 1.0: A0: 0.27 Amin :0.1 Amax: 1.5']
        assert_equals(return_error, expected_error)

    def test_29(self):
        '''
        test if Amin >= Amax
        '''

        self.Amin = '0.5'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage maximum has to be greater than surface coverage minumum: Amin :0.5 Amax: 0.5']
        assert_equals(return_error, expected_error)

    def test_30(self):
        '''
        test if Amin > A0
        '''

        self.Amin = '0.3'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage must be between Amin and Amax: A0 :0.27 Amin: 0.3 Amax: 0.5']
        assert_equals(return_error, expected_error)

    def test_31(self):
        '''
        test if Amax < A0
        '''

        self.Amax = '0.25'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['surface coverage must be between Amin and Amax: A0 :0.27 Amin: 0.1 Amax: 0.25']
        assert_equals(return_error, expected_error)

    def test_32(self):
        '''
        test if zevalmin < zfitmin
        '''

        self.zmin = '20.0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['error evaluation minimum has to be greater than zfit minimum: zevalmin = 7.5 zevalmax = 179.0']
        assert_equals(return_error, expected_error)

    def test_33(self):
        '''
        test if zfitmin > zfitmax
        '''

        self.zmax = '-20.0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['zfit maxiumum has to be greater than zfit minimum: zfitmin = -10.0 zfitmax = -20.0']
        assert_equals(return_error, expected_error)

    def test_34(self):
        '''
        test if zfitmin > zfit)
        '''

        self.zmin = '-2.0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['zfit value must be between zfit maxiumum and zfit minimum: zfit0 :-3.5 zfitmin: -2.0 zfitmax: 10.0']
        assert_equals(return_error, expected_error)

    def test_35(self):
        '''
        test if zfitmax < zfit)
        '''

        self.zmax = '-4.0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['zfit value must be between zfit maxiumum and zfit minimum: zfit0 :-3.5 zfitmin: -10.0 zfitmax: -4.0']
        assert_equals(return_error, expected_error)
        
    def test_36(self):
        '''
        test if zevalmin >= zevalmax
        '''

        self.zevalmax = '6.0'
        return_error = gui_mimic_sld_mol.run_module(
            self, test_filter=True)

        ''' check for value error '''
        expected_error = ['error evaluation maximum has to be greater than error evaluation minumum: zevalmin = 7.5 zevalmax = 6.0']
        assert_equals(return_error, expected_error)
            

    def tearDown(self):
        if os.path.exists(self.runname):
            shutil.rmtree(self.runname)

if __name__ == '__main__':
    main()
