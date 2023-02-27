'''
Driver method to run the sld_mol module
'''

import sys
import string
import os
import shutil
import time

import sassie.calculate.sld_mol.sld_mol as sld_mol
#import sld_mol as sld_mol
import sassie.interface.input_filter as input_filter
import sassie.interface.sld_mol.sld_mol_filter as sld_mol_filter
#import sld_mol_filter as sld_mol_filter
import multiprocessing


def user_variables(self, **kwargs):

    #### BEGIN USER EDIT
    #### BEGIN USER EDIT
    #### BEGIN USER EDIT

    self.runname = 'run_0'
    self.path = './'
    self.pdbfile = 'hiv1_gag.pdb'
    self.dcdfile = 'f12.dcd'
#    self.dcdfile = 'f12.pdb' 
    self.expdatafile = 'sld_gag.dat'
    self.outputfile = 'results_f12.dat'    
    self.runtype = '0'          #0 = average sld over all structures, 1 = best fit sld for each individual structure
    self.bulksld = '-0.567e-06'
    self.xon = '0'
#    self.numdregions = '2'
#    self.lowres = '1,150'
#    self.highres = '145,200'
    self.numdregions = '1'
    self.lowres = '1'
    self.highres = '200'
    self.dbin = '0.5'           #heavy atom profile
    self.width = '2.5'          #Gaussian smoothing
#    self.z0 = '3.8'
#    self.z0 = '1.7'
    self.z0 = '-3.6'
#    self.z0 = '-1.1'
#    self.z0 = '0.0'
    self.zmin = '-10.0'
    self.zmax = '10.0'
    self.zevalmin = '7.5'     #z-range to evaluate the error function
    self.zevalmax = '179.0'
#    self.zevalmin = '17.5'     #z-range to evaluate the error function
#    self.zevalmax = '189.0'    
    self.sldoffset = '0.0'
#    self.sldoffset = '10.0'    #offset to experimental SLD; zevalmin and zevalmax must be adjusted to accommodate offset    
#    self.A0 = '0.29'
#    self.A0 = '0.27'
#    self.A0 = '0.12'
    self.A0 = '0.23'
    self.Amin = '0.1'
    self.Amax = '0.5'
    self.plotflag = '2'         #0 == NO PLOT, 1 == matplotlib, 2 == Gnuplot
    self.sldfit = '0'

    self.testflag = False

    #### END USER EDIT
    #### END USER EDIT
    #### END USER EDIT


def test_variables(self, paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the gui_mimic_align class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    other_data_path = paths['other_data_path']
    module_data_path = paths['module_data_path']

    self.runname = 'run_0'
    self.path = ''
    self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag.pdb')
    self.dcdfile = os.path.join(dcd_data_path,'f12.dcd')
    self.expdatafile = os.path.join(other_data_path,'sld_gag.dat')
    self.outputfile = 'results_f12.dat'
    self.runtype = '0'
    self.bulksld = '-0.567e-06'
    self.xon = '0'
    self.numdregions = '0'
    self.lowres = '1'
    self.highres = '200'
    self.dbin = '0.5'
    self.width = '2.5'
    self.z0 = '-3.5'
    self.zmin = '-10.0'
    self.zmax = '10.0'
    self.zevalmin = '7.5'
    self.zevalmax = '179.0'
    self.sldoffset = '0.0'    
    self.A0 = '0.27'
    self.Amin = '0.1'
    self.Amax = '0.5'
    self.plotflag = '0'
    self.sldfit = '0'
    
    self.precision = 3
    self.testflag = True

def run_module(self, **kwargs):

    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables={}

    svariables['runname']		= (self.runname,'string')
    svariables['path']		= (self.path,'string')
    svariables['pdbfile']		= (self.pdbfile,'string')
    svariables['dcdfile']		= (self.dcdfile,'string')
    svariables['expdatafile']	= (self.expdatafile,'string')
    svariables['outputfile']	= (self.outputfile,'string')
    svariables['runtype']		= (self.runtype,'int')
    svariables['bulk_sld']		= (self.bulksld,'float')
    svariables['xon']		= (self.xon,'int')
    svariables['num_deut_regions']	= (self.numdregions,'int')
    svariables['deut_low_res']	= (self.lowres,'int_array')
    svariables['deut_high_res']	= (self.highres,'int_array')
    svariables['sldfit']		= (self.sldfit,'int')
    svariables['sldoffset']		= (self.sldoffset,'float')
    svariables['dbin']		= (self.dbin,'float')
    svariables['width']		= (self.width,'float')
    svariables['zfit0']		= (self.z0,'float')
    svariables['zfitmin']		= (self.zmin,'float')
    svariables['zfitmax']		= (self.zmax,'float')
    svariables['zevalmin']		= (self.zevalmin,'float')
    svariables['zevalmax']		= (self.zevalmax,'float')
    svariables['A0']		= (self.A0,'float')
    svariables['Amin']		= (self.Amin,'float')
    svariables['Amax']		= (self.Amax,'float')
    svariables['plotflag']		= (self.plotflag,'int') 


    error,self.variables=input_filter.type_check_and_convert(svariables)

    if(len(error)>0):
        print 'error = ',error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = sld_mol_filter.check_sld_mol(self.variables)
    except:
            error = sld_mol_filter.check_sld_mol(self.variables,no_file_check="true")
    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass

    runname=self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_sld_mol = sld_mol.sld_mol()
    this_sld_mol.main(self.variables, txtQueue)

class gui_mimic_sld_mol():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'sld_mol'

    def __init__(self, test, paths):

        if not test:
            user_variables(self)
        else:
            test_variables(self, paths)

        run_module(self)


if __name__ == '__main__':

    test = False  # option to run with test variables not implemented in 1.0.
    paths = None


# We are thinking of defining the install path so the gui mimic can be run from anywhere as long as it is called from that particular python
# That way, the test files will always be available to the user.
    if test:
        pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
        dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
        other_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'other_common') + os.path.sep
        module_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'interface', 'chi_square_filter') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_sld_mol(test, paths)
    print "time used: ", time.time() - start



