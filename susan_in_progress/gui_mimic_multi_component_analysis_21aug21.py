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
import multi_component_analysis_19aug21 as multi_component_analysis
import sassie.interface.input_filter_sasmol as input_filter
#import sassie.interface.multi_component_analysis.multi_component_analysis_filter as multi_component_analysis_filter
#no module input filer for now
import multiprocessing

def user_variables(self, **kwargs):

    #### user input ####
    #### user input ####
    #### user input ####

#NOTE: 'runname' is hardwired in run.utils.general_setup, so it can't be changed to 'run_name' here!!
    self.runname = 'run_0'
    self.path = './'
    self.number_of_contrast_points = '3' # must be >= number of components (2 in this case)
    self.number_of_components = '2'
    self.read_from_file = '0' #0=don't read from contrast calculator output file, 1=read from contrast calculator output file    
    self.output_file_name = os.path.join(self.path,'99_12_41.out')
    self.input_file_name = os.path.join(self.path,'input_contrast.txt') #this only matters if read_from_file = 1
    self.fraction_d2o = [u'0.99', u'0.12', u'0.41']  # 1 value for each contrast
    self.partial_specific_volume = [u'0.745', u'0.903']  # 1 value for each component
    self.izero = [u'11.8', u'0.6', u'0.17']  # 1 value for each contrast
    self.concentration = [u'3.7', u'3.6', u'3.1']  # 1 value for each contrast
    self.delta_rho = [[u'-3.2', u'-5.7'], [u'1.6', u'0.26'], [u'0.031', u'-1.74']]  # 2 values for each contrast since there are 2 components.  This will be superceded by the new values read from an input file if read_from_file = 1.
    self.match_point_flag = False

    self.stuhrman_parallel_axis_flag = False

    self.decomposition_flag = False

    self.stoichiometry_flag = True
    self.stoichiometry_variables = [self.fraction_d2o,self.partial_specific_volume,self.izero,self.concentration,self.delta_rho]  

    self.testflag = False
    
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

    self.runname = 'run_0'
    self.path = ''
    self.number_of_contrast_points = '3'
    self.number_of_components = '2'
    self.read_from_file = '0'
    self.output_file_name = os.path.join(other_data_path,'99_12_41.out')
    self.input_file_name = os.path.join(other_data_path,'input_contrast.txt')    
    self.fraction_d2o = [u'0.99', u'0.12', u'0.41']
    self.partial_specific_volume = [u'0.745', u'0.903']
    self.izero = [u'11.8', u'0.6', u'0.17']
    self.concentration = [u'3.7', u'3.6', u'3.1']
    self.delta_rho = [[u'-3.2', u'-5.7'], [u'1.6', u'0.26'], [u'0.031', u'-1.74']]
    self.match_point_flag = False
    self.stuhrman_parallel_axis_flag = False
    self.decomposition_flag = False
    self.stoichiometry_flag = True
    self.stoichiometry_variables = [self.fraction_d2o,self.partial_specific_volume,self.izero,self.concentration,self.delta_rho]  

    self.precision = 3 #not needed? What is the default?
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
    svariables['path'] = (self.path,'string')
    svariables['output_file_name'] = (self.output_file_name,'string')
    svariables['input_file_name'] = (self.input_file_name,'string')
    svariables['stoichiometry_flag'] = (self.stoichiometry_flag,'boolean')
    svariables['match_point_flag'] = (self.match_point_flag,'boolean')
    svariables['stuhrman_parallel_axis_flag'] = (self.stuhrman_parallel_axis_flag,'boolean')
    svariables['decomposition_flag'] = (self.decomposition_flag,'boolean')
    svariables['read_from_file'] = (self.read_from_file,'int')
    svariables['number_of_contrast_points'] = (self.number_of_contrast_points,'int')
    svariables['number_of_components'] = (self.number_of_components,'int')


    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print('error = ', error)
        if not(self.testflag):
            sys.exit()
        return error

#only stoichiometry variables are checked for now; others to be added later
#    try:
#        if kwargs['file_check']:
#            error = multi_component_analysis_filter.check_multi_component_analysis(self.variables, self.stoichiometry_variables)
#    except:
#            error = multi_component_analysis_filter.check_multi_component_analysis(self.variables, #self.stoichiometry_variables, no_file_check="true")

#    if(len(error) > 0):
#        print('error = ', error)
#        if not(self.testflag):
#            sys.exit()
#        return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass

    runname = self.variables['runname'][0]

#I don't want the directory tree to be deleted since I want to be able to put more than one output file in the directory for scenarios involving different contrasts.  But, then the sassie_json and sassie_log files don't get removed even if an output file is overwritten.
#    if os.path.exists(os.path.join(runname, self.module)):
#        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_multi_component_analysis = multi_component_analysis.multi_component_analysis()
    this_multi_component_analysis.main(self.variables, self.stoichiometry_variables, txtQueue)


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


