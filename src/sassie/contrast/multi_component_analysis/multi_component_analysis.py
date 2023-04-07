# -*- coding: utf-8 -*-

#    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#       MULTI-COMPONENT ANALYSIS
#
#       08/09/2021       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Multi-component Analysis** is the module that is used to analyze data from SANS contrast variation experiments. It can be used on model data prior to an experiment or on measured data.
    
    Included are 4 independent methods:
     
        1. Match Point Analysis
        2. Stuhrmann and Parallel Axis Theorem Analyses
        3. Decomposition Analysis
        4. Stoichiometry Analysis

    **Inputs:**
    
        The inputs depend on which method is used. The inputs are noted in the independent methods.
     
        The method is chosen based on which of the following flags is True:
    
        - match_point_flag: Match Point Analysis
        - stuhrmann_parallel_axis_flag: Sturhmann and Parallel Axis Theorem Analyses
        - decomposition_flag: Decomposition Analysis
        - stoichiometry_flag: Stoichiometry Analysis
        
        Only one flag at a time can be True.
    
    **Outputs:**
    
        The outputs depend on which method is used. The names of output files and their contents are noted in the independent methods.


    TODO: Need to encapsulate the flow chart, use cases, etc., flow that is hardcoded with if statements, for example.
    
    TODO: Need boilerplate line that shows the flow, i.e.,  module utilities, setup logging, unpack variables, run main, etc. This will be in all modules.

'''

import os
import io
import sys
import string
import locale
import time

import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.contrast.multi_component_analysis.read_contrast_output_files as read_contrast_output_files
import sassie.contrast.multi_component_analysis.match_point as match_point
import sassie.contrast.multi_component_analysis.stoichiometry as stoichiometry
import sassie.contrast.multi_component_analysis.stuhrmann_parallel_axis as stuhrmann_parallel_axis


if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'multi_component_analysis'


class module_variables():
    ''' Module variables class'''

    def __init__(self, parent=None):
        self.app = app


class multi_component_analysis_variables():
    ''' Multi-component analysis variables class'''

    def __init__(self, parent=None):
        pass


class multi_component_analysis():
    ''' Base class containing methods that control the program execution, unpacking of the module variables, intialization of the multi-component analysis variables and clean up after execution.'''

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):
        '''The main method that calls all of the other methods.'''

        self.module_variables = module_variables()

        self.multi_component_analysis_variables = multi_component_analysis_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.multi_component_analysis_initialization()

        if self.module_variables.match_point_flag:

            match_point.get_match_point(self)

        elif self.module_variables.stoichiometry_flag:

            stoichiometry.get_molecular_weights(self)

        elif self.module_variables.stuhrmann_parallel_axis_flag:

            stuhrmann_parallel_axis.parallel_axis(self)
            stuhrmann_parallel_axis.stuhrmann(self)

        self.epilogue()

        return

    def unpack_variables(self, variables):
        '''Method to unpack variables passed from the GUI.'''

        log = self.log
        mvars = self.module_variables
        log.debug('in unpack_variables')

# Unpack the variables that are used for all methods.
        mvars.run_name = variables['run_name'][0]
        mvars.output_file_name = variables['output_file_name'][0]
        mvars.number_of_contrast_points = variables['number_of_contrast_points'][0]
        mvars.stoichiometry_flag = variables['stoichiometry_flag'][0]
        mvars.match_point_flag = variables['match_point_flag'][0]
        mvars.stuhrmann_parallel_axis_flag = variables['stuhrmann_parallel_axis_flag'][0]
        mvars.decomposition_flag = variables['decomposition_flag'][0]
        mvars.fraction_d2o = variables['fraction_d2o'][0]

# Unpack the additional variables that are unique to each method. The read_from_file if statement may not be needed if the GUI is going to be populated with these values prior to running the main program.
        if mvars.stoichiometry_flag == True:
            #mvars.input_file_name = variables['input_file_name'][0]
            #mvars.read_from_file = variables['read_from_file'][0]
            mvars.izero = variables['izero'][0]
            mvars.concentration = variables['concentration'][0]
            mvars.partial_specific_volume = variables['partial_specific_volume'][0]
            mvars.number_of_components = variables['number_of_components'][0]

            # if mvars.read_from_file == True:
            # read_contrast_output_files.read_contrast_file(self)
#               # print('delta rho after reading from file: ', mvars.delta_rho)
            #log.debug('delta rho after reading from file: %s' % (str(mvars.delta_rho)))
            # pass
            # elif mvars.read_from_file == False:
            #mvars.delta_rho = variables['delta_rho'][0]
            mvars.delta_rho = variables['delta_rho'][0]

        elif mvars.match_point_flag == True:
            mvars.izero = variables['izero'][0]
            mvars.izero_error = variables['izero_error'][0]
            mvars.concentration = variables['concentration'][0]
            mvars.concentration_error = variables['concentration_error'][0]
            mvars.initial_match_point_guess = variables['initial_match_point_guess'][0]

        elif mvars.stuhrmann_parallel_axis_flag == True:
            #mvars.input_file_name = variables['input_file_name'][0]
            #mvars.read_from_file = variables['read_from_file'][0]
            mvars.partial_specific_volume = variables['partial_specific_volume'][0]
            mvars.molecular_weight = variables['molecular_weight'][0]
            mvars.number_of_components = variables['number_of_components'][0]
            mvars.radius_of_gyration = variables['radius_of_gyration'][0]
            mvars.radius_of_gyration_error = variables['radius_of_gyration_error'][0]

            # if mvars.read_from_file == True:
            #    #read_contrast_output_files.read_contrast_file(self)
#           #    # print('delta rho after reading from file: ', mvars.delta_rho)
            #    #log.debug('delta rho after reading from file: %s' % (str(mvars.delta_rho)))
            #    pass
            # elif mvars.read_from_file == False:
            #    mvars.delta_rho = variables['delta_rho'][0]
            mvars.delta_rho = variables['delta_rho'][0]

#        print(vars(mvars))

        log.debug(vars(mvars))

        return

    def multi_component_analysis_initialization(self):
        '''
        Method to initialize Multi-component analysis variables and write initial information to an output file.

        Parameters
        ----------

            run_name: string
                run name
            output_file_name: string
                user-specified output file name
            stoichiometry_flag: boolean
                flag to determine if stoichiometry analysis is being used
            match_point_flag: boolean
                flag to determine if match point analysis is being used 
            stuhrmann_parallel_axis_flag: boolean
                flag to determine if Stuhrmann and Parallel Axis methods are being used
            decomposition_flag: boolean
                flag to determine if decomposition analysis is being used
            partial_specific_volume: float array (dimension = number of components)
                partial specific volume of each component (used only if stuhrmann_parallel_axis_flag = True)
            molecular_weight: float array (dimension = number of components)
                molecular_weight of each component (used only if stuhrmann_parallel_axis_flag = True)

        Returns
        -------

            multi_component_analysis_path: string
                sub-path where output file will be written: run_name + \'multi_component_analysis\' + method-dependent sub-path
            outfile: string
                output file name (with full path): path + output_file_name
            volume_fraction: float array (dimension = number of components)
                volume fraction of each component (returned only if stuhrmann_parallel_axis_flag = True)

        '''

        log = self.log
        log.debug('in multi_component_analysis_initialization')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        mcavars = self.multi_component_analysis_variables


# Need to ask which method is being initialized to put a sub-path for the method used.
        if mvars.stoichiometry_flag == True:
            if (mvars.run_name[-1] == '/'):
                log.debug('run_name(1) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    'multi_component_analysis/stoichiometry/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))
            else:
                log.debug('run_name(2) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    '/multi_component_analysis/stoichiometry/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))
        elif mvars.match_point_flag == True:
            if (mvars.run_name[-1] == '/'):
                log.debug('run_name(1) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    'multi_component_analysis/match_point/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))
            else:
                log.debug('run_name(2) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    '/multi_component_analysis/match_point/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))
        elif mvars.stuhrmann_parallel_axis_flag == True:
            if (mvars.run_name[-1] == '/'):
                log.debug('run_name(1) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    'multi_component_analysis/stuhrmann_parallel_axis/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))
            else:
                log.debug('run_name(2) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    '/multi_component_analysis/stuhrmann_parallel_axis/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))
# calculate the volume fraction of each component from the partial specific volume and molecular weight
            total_volume = 0.0
            mcavars.volume_fraction = []
            for i in range(mvars.number_of_components):
                total_volume = total_volume + \
                    mvars.molecular_weight[i]*mvars.partial_specific_volume[i]
#            print('total volume: ', total_volume)
            for i in range(mvars.number_of_components):
                volume_fraction = mvars.molecular_weight[i] * \
                    mvars.partial_specific_volume[i]/total_volume
                mcavars.volume_fraction.append(volume_fraction)
#            print('volume fraction: ', mcavars.volume_fraction)
        direxist = os.path.exists(mcavars.multi_component_analysis_path)
        if(direxist == 0):
            os.system('mkdir -p ' + mcavars.multi_component_analysis_path)

# this is just to write the contrast calculator filename into the output file so we know this option was used
        mcavars.outfile = io.open(
            mcavars.multi_component_analysis_path+mvars.output_file_name, 'w')

        # if mvars.stoichiometry_flag == True or mvars.stuhrmann_parallel_axis_flag == True:
        #    if mvars.read_from_file == True:
        #        #pgui('\ninput file: %s' % (mvars.input_file_name))
        #        #mcavars.outfile.write('input file: ' + mvars.input_file_name +'\n')
        #        pass
        #    else:
        #        pgui('\ninput file: None')
        #        mcavars.outfile.write('input file: None\n')

# TODO: write out the individual input variables on separate lines and only write out the ones that are relevant for the method being used.  For stuhrmann_parallel_axis, include volume fraction.
        mcavars.outfile.write('input variables: ' + repr(vars(mvars)) + '\n')

        log.debug(vars(mvars))
        log.debug(vars(mcavars))

    def epilogue(self):
        '''Method to print out results and to move results to appropriate places.'''

        log = self.log
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        st = ''.join(['=' for x in range(60)])
        pgui('\n%s \n\n' % (st))
        pgui('MULTI-COMPONENT ANALYSIS IS DONE')

        time.sleep(1.0)

        return
