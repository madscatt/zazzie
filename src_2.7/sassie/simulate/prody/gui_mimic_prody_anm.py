import sys
import string
import os
import shutil
import time

import sassie.simulate.prody.prody_anm as prody_anm
#import prody_anm as prody_anm
import sassie.interface.input_filter as input_filter
import sassie.interface.prody.prody_filter as prody_filter
#import prody_filter as prody_filter
import multiprocessing

def user_variables(self, **kwargs):

    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT

    self.runname = 'run_0'
    self.pdbfile = 'hivr.pdb'
    self.number_modes = '5'  # number of normal modes to compute
    self.number_conformations_samp = '50'  # number of conformations to generate by random sampling of modes
    self.number_steps_traverse = '10'  # number of steps to tranverse each mode in both diretcions
    self.rmsd_conformations_samp = '1.0'  # average RMSD of randomly sampled conformations with respect to initial conformation
    self.rmsd_traverse = '1.5'  # maximum RMSD of conformations for trajectory from traversed mode with respect to initial conformation
    self.advanced_usage = '0'  # advanced usage option: 0=no; 1=yes
    self.advanced_usage_cmd = ' '  # user supplied ProDy command if advanced_usage=1

    self.testflag = False

    # END USER EDIT
    # END USER EDIT
    # END USER EDIT

def test_variables(self, paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    other_data_path = paths['other_data_path']
    module_data_path = paths['module_data_path']

    self.runname = 'run_0'
    self.pdbfile = os.path.join(pdb_data_path,'hivr.pdb')
    self.number_modes = '5'  
    self.number_conformations_samp = '50'  
    self.number_steps_traverse = '10'  
    self.rmsd_conformations_samp = '1.0'  
    self.rmsd_traverse = '1.5'  
    self.advanced_usage = '0'  
    self.advanced_usage_cmd = ' '  

    self.precision = 3
    self.testflag = True


def run_module(self, **kwargs):

    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['runname'] = (self.runname, 'string')
    svariables['pdbfile'] = (self.pdbfile, 'string')
    svariables['number_modes'] = (self.number_modes, 'int')
    svariables['number_conformations_samp'] = (self.number_conformations_samp, 'int')
    svariables['number_steps_traverse'] = (self.number_steps_traverse, 'int')
    svariables['rmsd_conformations_samp'] = (self.rmsd_conformations_samp, 'float')
    svariables['rmsd_traverse'] = (self.rmsd_traverse, 'float')
    svariables['advanced_usage'] = (self.advanced_usage, 'int')
    svariables['advanced_usage_cmd'] = (self.advanced_usage_cmd, 'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)

#    print 'variables: ', self.variables
    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = prody_filter.check_prody(self.variables)
    except:
            error = prody_filter.check_prody(self.variables, no_file_check="true")

#    print 'error: ', error
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

    runname = self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_prody = prody_anm.prody_anm()
    this_prody.main(self.variables, txtQueue)

class gui_mimic_prody():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'prody'

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
            __file__)), '..', '..', 'data', 'simulate', 'energy_minimization') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_prody(test, paths)
    print "time used: ", time.time() - start



    
