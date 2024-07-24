'''
Driver method to run the asaxs module
'''

import sys
import os
import shutil
import time

import sassie.analyze.asaxs.asaxs as asaxs
import sassie.interface.input_filter as input_filter
import sassie.interface.asaxs.asaxs_filter as asaxs_filter
import multiprocessing


def user_variables(self, **kwargs):

    #### user input ####
    #### user input ####
    #### user input ####

    self.run_name = 'run_0'
    
    self.pdb_file_name = os.path.join("test_pdb_inputs", "1E3B_label.pdb")

    self.dcd_file_flag = False

    if self.dcd_file_flag:
        self.dcd_name = "run_0.dcd"

    self.gaussian_noise_fraction = "0.001"

    self.number_of_q_values = "31"
    self.maximum_q_value = "1.0"

    self.number_of_d_values = "31"
    self.maximum_d_value = "1.0"

    self.nanocluster_flag = True

    if self.nanocluster_flag:
        # choice 1
        self.nanocluster_atom_name = "Au"

        self.nanocluster_radius = "7.0"

    self.central_energy_value = "12000"
    self.energy_range = "800"
    self.number_of_energy_values = "5"

    self.number_of_d_values_for_pair_distribution = "60"
    self.maximum_d_value_for_pair_distribution = "60"

    self.alternative_scattering_calculator_flag = False

    if self.alternative_scattering_calculator_flag:
        # choice 1
        self.crysol_flag = True
        self.sascalc_flag = False
        self.experimental_data_flag = False

        if self.crysol_file_flag:
            self.crysol_file_name = os.path.join("test_crysol_inputs", "DNA10Au2.int")


    self.plot_flag = False
    self.test_flag = False

#

    #### end user input ####
    #### end user input ####
    #### end user input ####


def test_variables(self, paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the gui_mimic_asaxs class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    module_data_path = paths['module_data_path']
    other_data_path = paths['other_data_path']

    self.run_name = 'run_0'



    self.plot_flag = False

    self.test_flag = True


def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['run_name'] = (self.run_name, 'string')
    svariables['pdb_file_name'] = (self.pdb_file_name, 'string')
    svariables['dcd_file_flag'] = (self.dcd_file_flag, 'boolean')

    if(self.dcd_file_flag):
        svariables['dcd_file_name'] = (self.dcd_file_name, 'string')

    svariables['gaussian_noise_fraction'] = (self.gaussian_noise_fraction, 'float')

    svariables['number_of_q_values'] = (self.number_of_q_values, 'integer')
    svariables['maximum_q_value'] = (self.maximum_q_value, 'float')

    svariables['number_of_d_values'] = (self.number_of_d_values, 'integer')
    svariables['maximum_d_value'] = (self.maximum_d_value, 'float')

    svariables['nanocluster_flag'] = (self.nanocluster_flag, 'boolean')

    if self.nanocluster_flag:
        # choice 1
        svariables['nanocluster_atom_name'] = (self.nanocluster_atom_name, 'string')
        svariables['nanocluster_radius'] = (self.nanocluster_radius, 'float')

    svariables['central_energy_value'] = (self.central_energy_value, 'float')
    svariables['energy_range'] = (self.energy_range, 'float')
    svariables['number_of_energy_values'] = (self.number_of_energy_values, 'integer')

    svariables['number_of_d_values_for_pair_distribution'] = (self.number_of_d_values_for_pair_distribution, 'integer')
    svariables['maximum_d_value_for_pair_distribution'] = (self.maximum_d_value_for_pair_distribution, 'float')

    svariables['alternative_scattering_calculator_flag'] = (self.alternative_scattering_calculator_flag, 'boolean')

    if self.alternative_scattering_calculator_flag:
        # choice 1
        svariables['crysol_file_flag'] = (self.crysol_file_flag, 'boolean')
        svariables['sascalc_file_flag'] = (self.sascalc_file_flag, 'boolean')
        svariables['experimental_data_flag'] = (self.experimental_data_flag, 'boolean')

        if self.crysol_file_flag:
            svariables['crysol_file_name'] = (self.crysol_file_name, 'string')
        elif self.sascalc_file_flag:
            svariables['experimenal_data_file_name'] = (self.experimenal_data_file_name, 'string')
        elif self.sascalc_file_flag:
            svariables['sascalc_file_name'] = (self.sascalc_file_name, 'string')

    svariables['plot_flag'] = (self.plot_flag, 'boolean')

    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print('error = ', error)
        if not(self.test_flag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = asaxs_filter.check_asaxs(self.variables)
    except:
        error = asaxs_filter.check_asaxs(
            self.variables, no_file_check="true")

    if(len(error) > 0):
        print('error = ', error)
        if not(self.test_flag):
            sys.exit()
        return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass

    run_name = self.variables['run_name'][0]

    if os.path.exists(os.path.join(run_name, self.module)):
        shutil.rmtree(os.path.join(run_name, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_asaxs = asaxs.asaxs()


#HERE
    import sys ; sys.exit()

    this_asaxs.main(self.variables, txtQueue)

    return

class gui_mimic_asaxs():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'asaxs'

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
        module_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'interface', 'align') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_asaxs(test, paths)
    print("time used: ", time.time() - start)
