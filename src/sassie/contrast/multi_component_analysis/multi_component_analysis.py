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
#       08/09/2021       --      initial coding               :  Susan Krueger
#       04/03/2023       --      python 3 coding              :  Joseph E. Curtis
#       06/27/2023       --      added decomposition variables:  Susan Krueger
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
import numpy

import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
#import sassie.contrast.multi_component_analysis.read_contrast_output_files as read_contrast_output_files
import sassie.contrast.multi_component_analysis.match_point as match_point
#import match_point as match_point
import sassie.contrast.multi_component_analysis.stoichiometry as stoichiometry
import sassie.contrast.multi_component_analysis.stuhrmann_parallel_axis as stuhrmann_parallel_axis
import sassie.contrast.multi_component_analysis.decomposition as decomposition
#import decomposition as decomposition


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

        elif self.module_variables.decomposition_flag:

            decomposition.get_composite_scattering_intensities(self)

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
# The read_from_contrast_calculator_output_file variables won't be needed since this is going to be done at the GUI level prior to running the main program.
        #mvars.read_from_contrast_calculator_output_file = variables['read_from_contrast_calculator_output_file'][0]
        # if mvars.read_from_contrast_calculator_output_file:
        #    mvars.read_from_contrast_calculator_output_file = variables['read_from_contrast_calculator_output_file'][0]

# Unpack the additional variables that are unique to each method.
        if mvars.stoichiometry_flag == True:
            mvars.izero = variables['izero'][0]
            mvars.concentration = variables['concentration'][0]
            mvars.partial_specific_volume = variables['partial_specific_volume'][0]
            mvars.number_of_components = variables['number_of_components'][0]

            # if mvars.read_from_contrast_calculator_output_file == True:
            # read_contrast_output_files.read_contrast_file(self)
#               # print('delta rho after reading from file: ', mvars.delta_rho)
            #log.debug('delta rho after reading from file: %s' % (str(mvars.delta_rho)))
            # pass
            # elif mvars.read_from_contrast_calculator_output_file == False:
            #mvars.delta_rho = variables['delta_rho'][0]
            mvars.delta_rho = variables['delta_rho'][0]

        elif mvars.match_point_flag == True:
            mvars.izero = variables['izero'][0]
            mvars.izero_error = variables['izero_error'][0]
            mvars.concentration = variables['concentration'][0]
            mvars.concentration_error = variables['concentration_error'][0]

        elif mvars.stuhrmann_parallel_axis_flag == True:
            mvars.partial_specific_volume = variables['partial_specific_volume'][0]
            mvars.molecular_weight = variables['molecular_weight'][0]
            mvars.number_of_components = variables['number_of_components'][0]
            mvars.radius_of_gyration = variables['radius_of_gyration'][0]
            mvars.radius_of_gyration_error = variables['radius_of_gyration_error'][0]

            # if mvars.read_from_contrast_calculator_output_file == True:
            #    #read_contrast_output_files.read_contrast_file(self)
#           #    # print('delta rho after reading from file: ', mvars.delta_rho)
            #    #log.debug('delta rho after reading from file: %s' % (str(mvars.delta_rho)))
            #    pass
            # elif mvars.read_from_file == False:
            #    mvars.delta_rho = variables['delta_rho'][0]
            mvars.delta_rho = variables['delta_rho'][0]

        elif mvars.decomposition_flag == True:
            mvars.partial_specific_volume = variables['partial_specific_volume'][0]
            mvars.molecular_weight = variables['molecular_weight'][0]
            mvars.number_of_components = variables['number_of_components'][0]
            mvars.delta_rho = variables['delta_rho'][0]
            mvars.concentration = variables["concentration"][0]
            # mvars.concentration_error is not currently used
            #mvars.concentration_error = variables["concentration_error"][0]
            mvars.data_file_name = variables["data_file_name"][0]
            mvars.q_rg_limit_guinier = variables["q_rg_limit_guinier"][0]
            mvars.starting_data_point_guinier = variables["starting_data_point_guinier"][0]
            mvars.initial_points_to_use_guinier = variables["initial_points_to_use_guinier"][0]
            mvars.refine_scale_factor_flag = variables["refine_scale_factor_flag"][0]

#        print(vars(mvars))

        log.debug(vars(mvars))

        return

    def multi_component_analysis_initialization(self):
        r'''
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
        partial_specific_volume: float array (dimension = number_of_components)
            partial specific volume of each component (used only if stuhrmann_parallel_axis_flag or decomposition_flag = True)
        molecular_weight: float array (dimension = number_of_components)
            molecular_weight of each component (used only if stuhrmann_parallel_axis_flag or decomposition_flag = True)
        concentration: float array (dimension = number_of_contrast_points)
            total concentration of the complex at each fraction D\ :sub:`2`\ O  (used only if decomposition_flag = True)
        data_file_name: string array (dimension = number_of_contrast_points)
            contrast variation data file name at each fraction D\ :sub:`2`\ O  (used only if decomposition_flag = True)

        Returns
        -------

        multi_component_analysis_path: string
            sub-path where output file will be written: run_name + \'multi_component_analysis\' + method-dependent sub-path
        outfile: string
            output file name (with full path): path + output_file_name
        volume_fraction: float array (dimension = number_of_components)
            volume fraction of each component (returned only if stuhrmann_parallel_axis_flag or decomposition_flag = True)
        initial_scale_factor:  float array (dimension = number_of_contrast_points)
            initial scale factor for the data at each fraction D\ :sub:`2`\ O  (returned only if decomposition_flag = True)
        scale_factor:  float array (dimension = number_of_contrast_points)
            scale factor for the data at each fraction D\ :sub:`2`\ O  that is the same as the initial scale factor before the Guinier analysis is performed (returned only if decomposition_flag = True)
        delta_rho_v:  float array (dimension = number_of_contrast_points)
            :math:`\Delta \rho V` at each fraction D\ :sub:`2`\ O  as defined in the Guinier analysis helper program (returned only if decomposition_flag = True)
        composite_intensity_file_name:  string array (dimension = 3)
            names of the composite scattering intensity output files (returned only if decomposition_flag = True)
        rescaled_data_file_name:  string array (dimension = number_of_contrast_points)
            names of the rescaled data output files (returned only if decomposition_flag = True)
        calculated_data_file_name:  string array (dimension = number_of_contrast_points)
            names of the calculated data output files (returned only if decomposition_flag = True)

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
        elif mvars.decomposition_flag == True:
            if (mvars.run_name[-1] == '/'):
                log.debug('run_name(1) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    'multi_component_analysis/decomposition/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))
            else:
                log.debug('run_name(2) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + \
                    '/multi_component_analysis/decomposition/'
                log.debug('multi_component_analysis_path = %s' %
                          (mcavars.multi_component_analysis_path))

# Calculate the volume fraction of each component from the partial specific volume and molecular weight if stuhrmann_parallel_axis_flag or decomposition_flag = True
# Note:  Avogadro's number cancels out when calculating the volume fraction, i.e., V1/(V1+V2)
        if mvars.stuhrmann_parallel_axis_flag == True or mvars.decomposition_flag == True:
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
# TODO:  Add a check to make sure that the volume fractions add up to 1?

# Additional initialization of mcavars if decomposition_flag = True
        if mvars.decomposition_flag == True:

            # initial and refined scale factor for each contrast based on the ratio of the concentrations. The scale factors are normalized to that of the first data set, i.e., scale_factor[i] = c[0]/c[i]. The initial scale factor is based on the input concentrations. The refined scale factor starts out the same as the initial scale factor but is refined later if refine_scale_factor_flag = True. We want to keep track of the initial scale factor to write to output file.
            # TODO: Now the normalization is to the first data set by default. Allow the user to specify which data set?
            mcavars.initial_scale_factor = []
            mcavars.scale_factor = []
            for i in range(mvars.number_of_contrast_points):
                ratio = mvars.concentration[0]/mvars.concentration[i]
                mcavars.initial_scale_factor.append(ratio)
                mcavars.scale_factor.append(ratio)
            #print('initial scale factor from concentrations: ', mcavars.initial_scale_factor)

            # calculate delta_rho_v = drho1*vf1 + drho2*vf2 for the scale factor refinement
            # I(0) = n*(delta_rho*volume)**2 is found from the Guinier analysis and then the scale factor is found using the ratio I(0)[0]/I(0)[i]; volume fraction is used here instead of volume since vf1 = v1/(v1+v2) but the denominator cancels out when calculating the ratio.

            mcavars.delta_rho_v = numpy.zeros(
                mvars.number_of_contrast_points, float)
#           print('delta_rho_v: ', mcavars.delta_rho_v)
#           delta_rho_v = vf1*delta_rho1 + vf2*delta_rho2
            for i in range(mvars.number_of_contrast_points):
                mcavars.delta_rho_v[i] = mcavars.volume_fraction[0]*mvars.delta_rho[i][0] + \
                    mcavars.volume_fraction[1]*mvars.delta_rho[i][1]
#               print(i, mcavars.delta_rho_v[i])
#           print('delta_rho_v: ', mcavars.delta_rho_v)

            # define composite intensity file names
            # we want the first few characters to be i11, i12 and i22; the run name has been added to give the user some flexibility in file naming.
# NOTE:  run name is maybe not needed because the output files will be in the output path, which includes the run name.  But, we include it in sascalc output files, monte carlo output files, etc., so it is included here.
            mcavars.composite_intensity_file_name = [
                'i11_'+mvars.run_name+'.dat', 'i12_'+mvars.run_name+'.dat', 'i22_'+mvars.run_name+'.dat']
            #print('composite intensity file name: ', mcavars.composite_intensity_file_name)

            # define rescaled and calculated data file names based on the data file names
            # first, strip the input path from the file name since it may be different from the output file path
            stripped_data_file_name = []
            for i in range(len(mvars.data_file_name)):
                numslash = mvars.data_file_name[i].count('/')
                if numslash == 0:
                    stripped_file = mvars.data_file_name[i]
                else:
                    groups = mvars.data_file_name[i].split('/')
#                   print('groups: ', groups)
                    stripped_file = ('.'.join(groups[numslash:]))
                stripped_data_file_name.append(stripped_file)
            #print('stripped file name: ', stripped_data_file_name)

            # then, get new output file names from the stripped data file name (output path is added later)
            mcavars.rescaled_data_file_name = []
            mcavars.calculated_data_file_name = []
# find the number of '.' in each stripped file name and locate the last one so that more characters can be added to the filename at that location
            for i in range(mvars.number_of_contrast_points):
                #print('data file name: ', data_file_name[i])
                numdots = stripped_data_file_name[i].count('.')
#               print('numdots: ', numdots)
                if numdots == 0:
                    rescaled_file = stripped_data_file_name[i] + '_rescaled'
                    calculated_file = stripped_data_file_name[i] + \
                        '_calculated'
#                   print('rescaled file, calculated file: ', rescaled_file, calculated_file)
                else:
                    groups = stripped_data_file_name[i].split('.')
#                   print('groups: ', groups)
                    rescaled_file = ('.'.join(
                        groups[:numdots]) + '_rescaled.' + '.'.join(groups[numdots:]))
                    calculated_file = ('.'.join(
                        groups[:numdots]) + '_calculated.' + '.'.join(groups[numdots:]))
#                   print('rescaled file, calculated file: ', rescaled_file, calculated_file)
                mcavars.rescaled_data_file_name.append(rescaled_file)
                mcavars.calculated_data_file_name.append(calculated_file)
            #print('data file, rescaled data file: ', mvars.data_file_name, mcavars.rescaled_data_file_name)

        # check for existence of output file path and create if necessary

        direxist = os.path.exists(mcavars.multi_component_analysis_path)
        if(direxist == 0):
            os.system('mkdir -p ' + mcavars.multi_component_analysis_path)

# open the general output file for writing
        mcavars.outfile = io.open(
            mcavars.multi_component_analysis_path+mvars.output_file_name, 'w')

        # this is just to write the contrast calculator filename into the output file so we know this option was used.  Commented out since we are no longer unpacking the contrast_calculator_output_file variables.
        # if mvars.stoichiometry_flag == True or mvars.stuhrmann_parallel_axis_flag == True:
        #    if mvars.read_from_contrast_calculator_output_file == True:
        #        #pgui('\ncontrast calculator output file used: %s' % (mvars.contrast_calculator_output_file_name))
        #        #mcavars.outfile.write('contrast calculator output file used: ' + mvars.contrast_calculator_output_file_name +'\n')
        #        pass
        #    else:
        #        pgui('\ncontrast calculator output file used: None')
        #        mcavars.outfile.write('contrast calculator output file used: None\n')

        # write the input variables in the output file
# TODO: write out the individual input variables on separate lines to make it easier to read. Perhaps make a separate method for this?
        mcavars.outfile.write('input variables: ' + repr(vars(mvars)) + '\n')

        log.debug(vars(mvars))
        log.debug(vars(mcavars))

        return

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
