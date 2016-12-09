'''
Driver method to run the extract_utilities module
'''

import sys
import os
import shutil
import time
import sassie.tools.extract_utilities.extract_utilities as extract_utilities
#import extract_utilities as extract_utilities
import sassie.interface.input_filter as input_filter
import sassie.interface.extract_utilities.extract_utilities_filter as extract_utilities_filter
#import extract_utilities_filter as extract_utilities_filter
import multiprocessing


def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname = 'run_0'
    self.path = './'
    self.pdb_filename = os.path.join(self.path, 'hiv1_gag.pdb')
    self.trajectory_filename = os.path.join(self.path, 'hiv1_gag_20_frames.dcd')
    self.option = 'weight_file'
#    self.option = 'text_file'
#    self.option ='single_frame'
#    self.option ='sampling_frequency'
#    self.option ='all'
#    self.option ='range'
#    self.local_value = '1'
#    self.local_value = '2-5'
    self.local_value = os.path.join(self.path, 'hiv1_gag_weight_file.txt')
#    self.local_value = os.path.join(self.path, 'hiv1_gag_text_file.txt')
    self.output_filename = 'chosen_weights.dcd'
#    self.output_filename = 'chosen_text.dcd'
#    self.output_filename = 'chosen_range.dcd'
#    self.output_filename = 'single_frame.dcd' 
#    self.output_filename = 'periodic.dcd'    
#    self.output_filename = 'all.dcd'
    self.extract_trajectory = True
    self.extract_sas = True
    self.sas_type = '0'
#    self.sas_type = '1'
#    self.sas_type = '2'
#    self.sas_type = '3'
#    self.sas_paths = os.path.join(self.path, 'hiv1_gag_0', 'sascalc', 'neutron_D2Op_100')
    self.sas_paths = os.path.join(self.path, 'hiv1_gag_0', 'sascalc', 'neutron_D2Op_100')+','+os.path.join(self.path, 'hiv1_gag_0', 'sascalc', 'neutron_D2Op_0') 
#    self.sas_paths = os.path.join(self.path,'hiv1_gag_0','xtal2sas')
#    self.sas_paths = os.path.join(self.path,'hiv1_gag_0','cryson')
#    self.sas_paths = os.path.join(self.path,'hiv1_gag_0','crysol')


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
    self.path = ''
    self.pdb_filename = os.path.join(pdb_data_path, 'hiv1_gag.pdb')
    self.trajectory_filename = os.path.join(dcd_data_path, 'hiv1_gag_20_frames.dcd')
    self.option = 'weight_file'
    self.local_value = os.path.join(other_data_path, 'hiv1_gag_weight_file.txt')
    self.output_filename = 'chosen_weights.dcd'
    self.extract_trajectory = True
    self.extract_sas = True
    self.sas_type = '0'
    self.sas_paths = os.path.join(other_data_path, 'hiv1_gag_0', 'sascalc', 'neutron_D2Op_100')+','+os.path.join(other_data_path, 'hiv1_gag_0', 'sascalc', 'neutron_D2Op_0')


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
    svariables['pdb_filename'] = (self.pdb_filename, 'string')
    svariables['trajectory_filename'] = (self.trajectory_filename, 'string')
    svariables['option'] = (self.option, 'string')
    svariables['local_value'] = (self.local_value, 'string')
    svariables['output_filename'] = (self.output_filename, 'string')
    svariables['extract_trajectory'] = (self.extract_trajectory, 'boolean')
    svariables['extract_sas'] = (self.extract_sas, 'boolean')
    svariables['path'] = (self.path, 'string')
    svariables['sas_type'] = (self.sas_type, 'int')
    svariables['sas_paths'] = (self.sas_paths, 'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)

    if(len(error) > 0):
        print 'error = ', error
#        sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = extract_utilities_filter.check_extract_utilities(
                self.variables)
    except:
        error = extract_utilities_filter.check_extract_utilities(
            self.variables, no_file_check="true")

    if(len(error) > 0):
        print 'error = ', error
#        sys.exit()
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
    this_extract_utilities = extract_utilities.extract_utilities()
    this_extract_utilities.main(self.variables, txtQueue)

#    txtQueue = multiprocessing.JoinableQueue()
#    extract_utilities.extract_data(self.variables, txtQueue)


class gui_mimic_extract_utilities():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'extract_utilities'

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
            __file__)), '..', '..', 'data', 'interface', 'extract_utilities') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_extract_utilities(test, paths)
    print "time used: ", time.time() - start
