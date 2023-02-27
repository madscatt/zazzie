'''
Driver method to run the merge_utilities module
'''

import sys
import os
import shutil
import time
import sassie.tools.merge_utilities.merge_utilities as merge_utilities
#import merge_utilities as merge_utilities
import sassie.interface.input_filter as input_filter
import sassie.interface.merge_utilities.merge_utilities_filter as merge_utilities_filter
#import merge_utilities_filter as merge_utilities_filter
import multiprocessing


def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname = 'run_0'
    self.pdb_file = os.path.join('./', 'hiv1_gag.pdb')
    self.trajectory_names = os.path.join('./','run_m1.dcd')+','+os.path.join('./','run_m2.dcd')
#    self.trajectory_names = os.path.join('./','run_m1.pdb')+','+os.path.join('./','run_m2.pdb')    
    self.output_filename = 'merged_run_m1_m2.dcd'
#    self.output_filename = 'merged_run_m1_m2.pdb'    
    self.number_of_runs = '2'
    self.local_value = '0'
#    self.local_value = os.path.join('./','weights_file_m1.txt')+','+os.path.join('./','weights_file_m2.txt')
#    self.local_value = '2'
    self.merge_option = '0'
    self.merge_type_option = '0'
    self.sas_type = '0'
#    self.sas_paths = os.path.join('./', 'merge_files_0','run_0','sascalc', 'neutron_D2Op_80')+','+os.path.join('./', 'merge_files_0','run_1', 'sascalc', 'neutron_D2Op_80') 
    self.sas_paths = os.path.join('./', 'merge_files_0','run_0','sascalc')+','+os.path.join('./', 'merge_files_0','run_1', 'sascalc')    
#    self.sas_paths = os.path.join('./', 'merge_files_0','run_0','xtal2sas')+','+os.path.join('./', 'merge_files_0','run_1','xtal2sas')
#    self.sas_paths = os.path.join('./', 'merge_files_0','run_0','cryson')+','+os.path.join('./', 'merge_files_0','run_1','cryson')
#    self.sas_paths = os.path.join('./', 'merge_files_0','run_0','crysol')+','+os.path.join('./', 'merge_files_0','run_1','crysol')

    self.testflag = False

    ### END USER INPUT ###
    ### END USER INPUT ###
    ### END USER INPUT ###


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
    self.pdb_file = os.path.join(pdb_data_path, 'hiv1_gag.pdb')
    self.trajectory_names = os.path.join(dcd_data_path,'run_m1.dcd')+','+os.path.join(dcd_data_path,'run_m2.dcd')
    self.output_filename = 'all.dcd'
    self.number_of_runs = '2'
    self.local_value = '0'
    self.merge_option = '0'
    self.merge_type_option = '0'
    self.sas_type = '1'
    self.sas_paths = os.path.join(other_data_path, 'merge_files_0','run_0','xtal2sas')+','+os.path.join(other_data_path, 'merge_files_0','run_1','xtal2sas')
    
    self.testflag = True
    self.precision = 3


def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['runname'] = (self.runname, 'string')
    svariables['pdb_file'] = (self.pdb_file, 'string')
    svariables['trajectory_names'] = (self.trajectory_names, 'string')
    svariables['output_filename'] = (self.output_filename, 'string')
    svariables['number_of_runs'] = (self.number_of_runs, 'int')
    svariables['local_value'] = (self.local_value, 'string')
    svariables['merge_option'] = (self.merge_option, 'int')
    svariables['merge_type_option'] = (self.merge_type_option, 'int')
    svariables['sas_type'] = (self.sas_type, 'int')
    svariables['sas_paths'] = (self.sas_paths, 'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)

    if(len(error) > 0):
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = merge_utilities_filter.check_merge_utilities(
                self.variables)
    except:
        error = merge_utilities_filter.check_merge_utilities(
            self.variables, no_file_check="true")

    if(len(error) > 0):
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

    path = os.path.join(runname, self.module)

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_merge_utilities = merge_utilities.merge_utilities()
    this_merge_utilities.main(self.variables, txtQueue)
    

class gui_mimic_merge_utilities():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'merge_utilities'

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
            __file__)), '..', '..', 'data', 'interface', 'merge_utilities') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_merge_utilities(test, paths)
    print "time used: ", time.time() - start
