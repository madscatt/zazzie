'''
Driver method to run the energy minimization module
'''

import sys
import string
import os
import shutil
import time

import sassie.util.sasconfig as sasconfig
#import sasconfig as sasconfig
import sassie.interface.input_filter as input_filter
import sassie.simulate.energy_minimization.energy_minimization as energy_minimization
#import energy_minimization as energy_minimization
import sassie.interface.energy_minimization.minimize_filter as minimize_filter
#import minimize_filter as minimize_filter
import multiprocessing


def user_variables(self, **kwargs):

    ### BEGIN USER EDIT ###
    ### BEGIN USER EDIT ###
    ### BEGIN USER EDIT ###

    self.runname = 'run_0'
    self.infile = 'ten_mer.pdb'
    self.pdbfile = 'ten_mer2.pdb'
    self.outfile = 'min_ten_mer.dcd'
    self.nsteps = '100'
    self.resparmfile = os.path.join(
        sasconfig.__bin_path__, 'toppar', 'par_all27_prot_na.inp')
    self.psffile = 'ten_mer.psf'
    self.ncpu = '2'
    self.keepout = '1'
    self.dcdfreq = '20'
    self.infiletype = 'pdb'  # Not used; determined by file name suffix
    self.md = '1'
    self.mdsteps = '20'
    self.dielect = '80'
    self.temperature = '300.0'

    self.use_external_input_file = False
    self.external_input_file = 'external_input.inp'

    self.velocity_restart_file = "False"
    self.extended_system_restart_file = "False"

    self.testflag = False

    ### END USER EDIT ###
    ### END USER EDIT ###
    ### END USER EDIT ###


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
    self.infile = os.path.join(pdb_data_path, 'ten_mer.pdb')
    self.pdbfile = os.path.join(pdb_data_path, 'ten_mer.pdb')
    self.outfile = 'min_ten_mer.dcd'
    self.nsteps = '100'
    self.resparmfile = os.path.join(
        sasconfig.__bin_path__, 'toppar', 'par_all27_prot_na.inp')
    self.psffile = os.path.join(other_data_path, 'ten_mer.psf')
    self.ncpu = '2'
    self.keepout = '0'
    self.dcdfreq = '20'
    self.infiletype = 'pdb'  # not used; determined by file name suffix
    self.md = '0'
    self.mdsteps = '20'
    self.dielect = '80'
    self.temperature = '300.0'

    self.use_external_input_file = False
    self.external_input_file = 'external_input_2.inp'

    self.velocity_restart_file = "False"
    self.extended_system_restart_file = "False"

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
    svariables['infile'] = (self.infile, 'string')
    svariables['pdbfile'] = (self.pdbfile, 'string')
    svariables['outfile'] = (self.outfile, 'string')
    svariables['nsteps'] = (self.nsteps, 'int')
    svariables['resparmfile'] = (self.resparmfile, 'string')
    svariables['psffile'] = (self.psffile, 'string')
    svariables['ncpu'] = (self.ncpu, 'int')
    svariables['keepout'] = (self.keepout, 'int')
    svariables['dcdfreq'] = (self.dcdfreq, 'int')
    svariables['infiletype'] = (self.infiletype, 'string')
    svariables['md'] = (self.md, 'int')
    svariables['mdsteps'] = (self.mdsteps, 'int')
    svariables['dielect'] = (self.dielect, 'float')
    svariables['temperature'] = (self.temperature, 'float')

    svariables['use_external_input_file'] = (
        self.use_external_input_file, 'boolean')
    svariables['external_input_file'] = (self.external_input_file, 'string')

    svariables['velocity_restart_file'] = (
        self.velocity_restart_file, 'string')
    svariables['extended_system_restart_file'] = (
        self.extended_system_restart_file, 'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)

#    print 'variables: ', self.variables
    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = minimize_filter.check_minimize(self.variables)
    except:
        error = minimize_filter.check_minimize(
            self.variables, no_file_check="true")

    print 'error: ', error
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
    this_energy_minimization = energy_minimization.energy_minimization()
    this_energy_minimization.main(self.variables, txtQueue)

#   NOTE:  getting broken pipe error (not always, but most of the time, 2.0 only, after epilogue finishes, "time used" also appears before error msg):
#   File "/opt/local/anacondaz/lib/python2.7/multiprocessing/queues.py", line 268, in _feed
#     send(obj)
#   IOError: [Errno 32] Broken pipe
#   Tests run OK even when error message appears

class gui_mimic_energy_minimization():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'energy_minimization'

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
    run_gui = gui_mimic_energy_minimization(test, paths)
    print "time used: ", time.time() - start
