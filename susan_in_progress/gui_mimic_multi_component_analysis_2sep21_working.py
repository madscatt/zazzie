# -*- coding: utf-8 -*-
'''
Driver method to run the multi_component_analysis module
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import shutil
import time

#import sassie.contrast.multi_component_analysis.multi_component_analysis as multi_component_analysis
import multi_component_analysis_2sep21_working as multi_component_analysis
import sassie.interface.input_filter_sasmol as input_filter
#import sassie.interface.multi_component_analysis.multi_component_analysis_filter as multi_component_analysis_filter
import multi_component_analysis_filter_2sep21_working as multi_component_analysis_filter
import multiprocessing

def user_variables(self, **kwargs):

    #### user input ####
    #### user input ####
    #### user input ####

    self.run_name = 'run_0'
    self.path = './'
    self.number_of_contrast_points = '3' # must be >= number of components (2 in this case)
    self.number_of_components = '2'
    self.read_from_file = False
    self.output_file_name = os.path.join(self.path,'99_12_41.out')
    self.input_file_name = os.path.join(self.path,'input_contrast.txt') #this only matters if read_from_file = True
    self.match_point_flag = False
    self.stuhrmann_parallel_axis_flag = False
    self.decomposition_flag = False
    self.stoichiometry_flag = True
    
#at some point we will need a way to know which variables will be used for the chosen method.  Otherwise, a blank string may give an error?  Check how the input filter handles a blank string.
    self.fraction_d2o = '0.99, 0.12, 0.41'  # 1 value for each contrast
    self.izero = '11.8, 0.6, 0.17'  # 1 value for each contrast
    self.concentration = '3.7, 3.6, 3.1'  # 1 value for each contrast
#if read_from_file = True, then partial_specific_volume needs to be input IN THE SAME ORDER as the delta_rho values will be read from the file.
    self.partial_specific_volume = '0.745, 0.903'  # 1 value for each component

#still making changes to input filter to accommodate a string of delta rho values
    self.delta_rho = [[u'-3.2', u'-5.7'], [u'1.6', u'0.26'], [u'0.031', u'-1.74']]  # 2 values for each contrast since there are 2 components.  This will be superceded by the new values read from an input file if read_from_file = True.

    self.test_flag = False
    
    #### end user input ####
    #### end user input ####
    #### end user input ####

def test_variables(self,paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the gui_mimic_align class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    module_data_path = paths['module_data_path']
    other_data_path = paths['other_data_path']

    self.run_name = 'run_0'
    self.path = ''
    self.number_of_contrast_points = '3'
    self.number_of_components = '2'
    self.read_from_file = True
    self.output_file_name = os.path.join(other_data_path,'99_12_41.out')
    self.input_file_name = os.path.join(other_data_path,'input_contrast.txt')
    self.match_point_flag = False
    self.stuhrmann_parallel_axis_flag = False
    self.decomposition_flag = False
    self.stoichiometry_flag = True
    self.fraction_d2o = '0.99, 0.12, 0.41' 
    self.izero = '11.8, 0.6, 0.17' 
    self.concentration = '3.7, 3.6, 3.1'  
    self.partial_specific_volume = '0.745, 0.903' 

    self.delta_rho = [[u'-3.2', u'-5.7'], [u'1.6', u'0.26'], [u'0.031', u'-1.74']]

    self.precision = 3 #not needed? What is the default?
    self.test_flag = True

def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables={}

    svariables['run_name'] = (self.run_name,'string')
    svariables['path'] = (self.path,'string')
    svariables['output_file_name'] = (self.output_file_name,'string')
    svariables['input_file_name'] = (self.input_file_name,'string')
    svariables['stoichiometry_flag'] = (self.stoichiometry_flag,'boolean')
    svariables['match_point_flag'] = (self.match_point_flag,'boolean')
    svariables['stuhrmann_parallel_axis_flag'] = (self.stuhrmann_parallel_axis_flag,'boolean')
    svariables['decomposition_flag'] = (self.decomposition_flag,'boolean')
    svariables['read_from_file'] = (self.read_from_file,'boolean')
    svariables['number_of_contrast_points'] = (self.number_of_contrast_points,'int')
    svariables['number_of_components'] = (self.number_of_components,'int')
    svariables['fraction_d2o'] = (self.fraction_d2o,'float_array')
    svariables['izero'] = (self.izero, 'float_array')
    svariables['concentration'] = (self.concentration, 'float_array')
    svariables['partial_specific_volume'] = (self.partial_specific_volume, 'float_array')


    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print('error = ', error)
        if not(self.test_flag):
            sys.exit()
        return error

    print('after input filter')
    print(self.variables)

#only stoichiometry variables will be checked for now; others to be added later
    try:
        if kwargs['file_check']:
            error = multi_component_analysis_filter.check_multi_component_analysis(self.variables, self.delta_rho)
    except:
            error = multi_component_analysis_filter.check_multi_component_analysis(self.variables, self.delta_rho, no_file_check="true")

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

#I don't want the directory tree to be deleted since I want to be able to put more than one output file in the directory for scenarios involving different contrasts.  But, then the sassie_json and sassie_log files don't get removed even if an output file is overwritten.
    if os.path.exists(os.path.join(run_name, self.module)):
        shutil.rmtree(os.path.join(run_name, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_multi_component_analysis = multi_component_analysis.multi_component_analysis()
    this_multi_component_analysis.main(self.variables, self.delta_rho, txtQueue)


class gui_mimic_multi_component_analysis():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'multi_component_analysis'

    def __init__(self, test, paths):

        if not test:
            user_variables(self)
        else:
            test_variables(self, paths)

        run_module(self)


if __name__ == '__main__':

    test = False  # option to run with test variables not implemented
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
    run_gui = gui_mimic_multi_component_analysis(test, paths)
    print('time used: ', time.time() - start)


