'''
Driver method to run the em_to_sas module
'''

import sys
import string
import os
import shutil
import time

import sassie.calculate.em_to_sas.em_to_sas as em_to_sas
#import em_to_sas as em_to_sas
import sassie.interface.input_filter as input_filter
import sassie.interface.em_to_sas.em_to_sas_filter as em_to_sas_filter
#import em_to_sas_filter as em_to_sas_filter
import multiprocessing

def user_variables(self, **kwargs):

    #### BEGIN USER EDIT
    #### BEGIN USER EDIT
    #### BEGIN USER EDIT

    
    self.runname = 'run_0'
    self.emfiletype = '0'
    self.inputpath = './'
    self.emdensityfile = 'gag_complete.cube'
    self.threshold = '0.01'
    self.pdbfile = 'test.pdb'
    self.sasfile = 'test.iq'
    self.npoints = '51'
    self.qmax = '0.3'
    self.plotflag = '1'

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
    self.emfiletype = '1'
    self.inputpath = ''
    self.emdensityfile = os.path.join(other_data_path,'emd_1180.mrc')
    self.threshold = '1.4'
    self.pdbfile = 'test.pdb'
    self.sasfile = 'test.iq'
    self.npoints = '51'
    self.qmax = '0.3'
    self.plotflag = '0'

    self.testflag = True

def run_module(self, **kwargs):

    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables={}

    svariables['runname'] = (self.runname,'string')
    svariables['emfiletype'] = (self.emfiletype,'int')
    svariables['inputpath'] = (self.inputpath,'string')
    svariables['emdensityfile'] = (self.emdensityfile,'string')
    svariables['threshold'] = (self.threshold,'float')
    svariables['pdbfile'] = (self.pdbfile,'string')
    svariables['sasfile'] = (self.sasfile,'string')
    svariables['npoints'] = (self.npoints,'int')
    svariables['qmax'] = (self.qmax,'float')
    svariables['plotflag'] = (self.plotflag,'int')


    error, self.variables = input_filter.type_check_and_convert(svariables)

    if(len(error)>0):
        print 'error = ',error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = em_to_sas_filter.check_em_to_sas(self.variables)
    except:
            error = em_to_sas_filter.check_em_to_sas(self.variables,no_file_check="true")
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
    this_em_to_sas = em_to_sas.em_to_sas
    this_em_to_sas.main(self.variables, txtQueue)

class gui_mimic_em_to_sas():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'em_to_sas'

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
    run_gui = gui_mimic_em_to_sas(test, paths)
    print "time used: ", time.time() - start


