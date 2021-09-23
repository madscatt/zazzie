# -*- coding: utf-8 -*-
'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import io
import sys
import string
import locale
import time

import numpy
import scipy.optimize


import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.interface.input_filter_sasmol as input_filter
#import sassie.contrast.muli_component_analysis.read_contrast_output_files as read_contrast_output_files
import read_contrast_output_files_2sep21_working as read_contrast_output_files #this is older version of the method that isn't in the 2.0 format

#       MULTI-COMPONENT ANALYSIS
#
#       08/09/2021       --      initial coding         :   sk
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    MULTI-COMPONENT ANALYSIS is the module that contains the methods
    that are used to analyze data from contrast variation experiments.

'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'multi_component_analysis'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class multi_component_analysis_variables():

    def __init__(self, parent=None):
        pass

class multi_component_analysis():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, contrast_variables, txtOutput):

        self.module_variables = module_variables()

        self.multi_component_analysis_variables = multi_component_analysis_variables() 

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)  

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.process_input_variables(contrast_variables)

        self.run_utils.general_setup(self)

        self.multi_component_analysis_initialization()

        if self.module_variables.stoichiometry_flag:

            self.get_molecular_weights()

        self.epilogue()

        return



    def unpack_variables(self,variables):

        log = self.log
        mvars = self.module_variables
        print('in unpack_variables')

        mvars.run_name = variables['run_name'][0]
        mvars.output_file_name = variables['output_file_name'][0]
        mvars.input_file_name = variables['input_file_name'][0]
        mvars.read_from_file = variables['read_from_file'][0]
        mvars.number_of_contrast_points = variables['number_of_contrast_points'][0]
        mvars.number_of_components = variables['number_of_components'][0]
        mvars.stoichiometry_flag = variables['stoichiometry_flag'][0]
        mvars.match_point_flag = variables['match_point_flag'][0] 
        mvars.stuhrmann_parallel_axis_flag = variables['stuhrmann_parallel_axis_flag'][0]
        mvars.decomposition_flag = variables['decomposition_flag'][0]
        mvars.fraction_d2o = variables['fraction_d2o'][0]
        mvars.izero = variables['izero'][0]
        mvars.concentration = variables['concentration'][0]
        mvars.partial_specific_volume = variables['partial_specific_volume'][0]
        
        print(vars(mvars))

        log.debug(vars(mvars))

        return


    def process_input_variables(self,contrast_variables):

        log = self.log
        mvars = self.module_variables
        pgui = self.run_utils.print_gui

        print('in process_input_variables')

        print('contrast_variables: %s' %(str(contrast_variables)))
        print('len(contrast_variables): %i' % len(contrast_variables))
        log.debug('number of contrast points: %i' % (mvars.number_of_contrast_points))
        log.debug('number of components: %i' % (mvars.number_of_components))
        log.debug('read from file: %s' % (mvars.read_from_file))

#this will eventually be an if, elif depending on the chosen method
        if mvars.stoichiometry_flag == True:
            if mvars.read_from_file == True:
#               delta_rho_read_from_file = read_contrast_output_files.read_contrast_file(mvars.input_file_name, mvars.fraction_d2o, mvars.number_of_contrast_points, mvars.number_of_components)
#               delta_rho_read_from_file = read_contrast_output_files.read_contrast_file(self)
#remove this print statement
#               print('delta rho after reading from file: ', delta_rho_read_from_file)
#               log.debug('delta rho after reading from file: %s' % (str(delta_rho_read_from_file)))
                read_contrast_output_files.read_contrast_file(self)
#                print('delta rho after reading from file: ', mvars.delta_rho)
                log.debug('delta rho after reading from file: %s' % (str(mvars.delta_rho)))
#               mvars.delta_rho = delta_rho_read_from_file
            elif mvars.read_from_file == False:
                mvars.delta_rho = contrast_variables

        print(vars(mvars))
        log.debug(vars(mvars))

        return




#Will there be multiple initialization methods depending on which analysis is being done?

    def multi_component_analysis_initialization(self):
        '''
        method to prepare for get_molecular_weights
        '''

        log = self.log
        log.debug('in multi_component_analysis_initialization')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        mcavars = self.multi_component_analysis_variables

        pgui('setting up output file path')
        
# Need to ask which method is being initialized to put a sub-path for the method used, i.e., multi_component_analysis/method
        if mvars.stoichiometry_flag == True:
            if (mvars.run_name[-1] == '/'):
                log.debug('run_name(1) = %s' % (mvars.run_name))
                multi_component_analysis_path = mvars.run_name + 'multi_component_analysis/stoichiometry/'
                log.debug('multi_component_analysis_path = %s' % (multi_component_analysis_path))
            else:
                log.debug('run_name(2) = %s' % (mvars.run_name))
                multi_component_analysis_path = mvars.run_name + '/multi_component_analysis/stoichiometry/'
                log.debug('multi_component_analysis_path = %s' % (multi_component_analysis_path))

        direxist = os.path.exists(multi_component_analysis_path)
        if(direxist == 0):
            os.system('mkdir -p ' + multi_component_analysis_path)

        # this is just to write the contrast calculator filename into the output file so we know this option was used
        mcavars.outfile = io.open(multi_component_analysis_path+mvars.output_file_name, 'w')

        if mvars.read_from_file == True:
#            pgui('input file: %s' % (mvars.input_file_name))
            mcavars.outfile.write('input file: ' + mvars.input_file_name +'\n')
        else:
#            pgui('input file: None')
            mcavars.outfile.write('input file: None\n')
            
#TODO: write the other relevant input variables to the output file for future reference. 

        log.debug(vars(mvars))
        log.debug(vars(mcavars))


    def get_molecular_weights(self):
        '''
        method to calculate the molecular weights
        
        
        INPUT:  variable descriptions
        
        OUTPUT:
        '''
#mvars used:  fraction_d2o, number_of_contrast_points, fraction_d2o, partial_specific_volume, izero, concentration, delta_rho 
#avars used:  outfile

    # TODO: this function could be outside of this method. If keeping it inside, put it right above the curve fit. NOTE:  I tried the latter and it didn't work.  Error said that molecular_weight_1 was referenced before assignment.
    # This function is the I(0) equation for two components rearranged such that the x array represents the known coefficients of the equation and the y array is == 0. This function is satisfied at each contrast.
    #TODO: put the actual equation in the doc string

        def izero_function(x, molecular_weight_1, molecular_weight_2):
            return x[0] * molecular_weight_1**2 + x[1] * molecular_weight_1*molecular_weight_2 + x[2] * molecular_weight_2**2 - x[3] * molecular_weight_1 - x[3] * molecular_weight_2

        log = self.log
        log.debug('in get_molecular_weights')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        mcavars = self.multi_component_analysis_variables
        
        log.debug(vars(mvars))

        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        Na = 6.023  #10**23

        partial_specific_volume_1 = mvars.partial_specific_volume[0]
        partial_specific_volume_2 = mvars.partial_specific_volume[1]
        print ('partial_specific_volume_1,partial_specific_volume_2: ', partial_specific_volume_1, partial_specific_volume_2)

        # print ("number of contrast points: ", number_of_contrast_points)

        pgui('setting up I(0) equation coefficients')
        
        x = numpy.zeros(shape=(mvars.number_of_contrast_points, 4))

        # Here we are defining the x coeficients at each contrast for the set of simultaneous I(0) equations
        # NOTE that range is being used here since xrange in Python 2 == range in Python 3. So, we aren't doing exactly the same thing here as will be done in Python 3, as range will behave like Python 2's xrange in Python 3.  But, for this case, it shouldn't matter (if I understand correctly).
        for i in range(mvars.number_of_contrast_points):
            delta_rho_1 = locale.atof(mvars.delta_rho[i][0])
            delta_rho_2 = locale.atof(mvars.delta_rho[i][1])
            izero_1 = mvars.izero[i]
            concentration_1 = mvars.concentration[i]
        # print('delta_rho_1, delta_rho_2, izero, concentration: ', delta_rho_1, delta_rho_2, izero_1, concentration_1)
            x0 = delta_rho_1**2*partial_specific_volume_1**2
            x1 = 2*delta_rho_1*delta_rho_2*partial_specific_volume_1*partial_specific_volume_2
            x2 = delta_rho_2**2*partial_specific_volume_2**2
            x3 = izero_1*Na/concentration_1
            x[i] = (x0, x1, x2, x3)

        # print ('x: ', x)
        # print ('length of x: ', len(x))

        y = numpy.zeros(len(x))
        # print ('y: ', y)

        # transpose x to have correct inputs for curve_fit
        x = x.T  

        pgui('calculating molecular weights')
        

        '''
scipy.optimize.curve_fit is using scipy version 1.2.1 in Python 2.7.  There is no known change of syntax or usage if using scipy 1.7 under Python 3.

        '''

        # optimized_molecular_weights are the optimal values for the molecular weights so that the sum of the squared error of izero_function(x, optimized_molecular_weights) - y is minimized
        # covariance, correlation coefficient and standard deviations are not reported to the user since the option to define errors on the coefficients isn't used since errors on partial specific volume and delta_rho are difficult to estimate. The calculations are left here for reference.

        optimized_molecular_weights, covariance = scipy.optimize.curve_fit(izero_function, x, y)

#TODO: decide what results should be in the output file and what should be printed out to the GUI screen
#        pgui('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
#        avars.outfile.write('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
#        pgui('calculated covariance: '+str(covariance)+'\n')
#        avars.outfile.write('calculated covariance: '+str(covariance)+'\n')

        pgui('-------------------------------')
        mcavars.outfile.write('--------------------------------\n')
        pgui('Final Results\n')
        mcavars.outfile.write('Final Results\n')

        # convert molecular weights to kDa
        molecular_weight_1 = optimized_molecular_weights[0]*10**3
        molecular_weight_2 = optimized_molecular_weights[1]*10**3
        pgui('molecular_weight_1, molecular_weight_2 (kDa): '+str(molecular_weight_1)+'\t'+str(molecular_weight_2)+'\n')
        mcavars.outfile.write('molecular_weight_1, molecular_weight_2 (kDa): '+str(molecular_weight_1)+'\t'+str(molecular_weight_2)+'\n')

        # The covariance matrix will be infinite if there are only two contrasts, so we don't want to calculate standard deviations in this case; numpy.all tests whether all array elements along a given axis evaluate to True.
#        if numpy.all(covariance != numpy.inf):

            # rescale values in covariance matrix since molecular weights were converted to kDa
#            rescaled_covariance = covariance*10**6
#            pgui('covariance matrix: ' + str(rescaled_covariance)+'\n')
#            avars.outfile.write('covariance matrix: ' + str(rescaled_covariance)+'\n')
#            # compute one standard deviation of the molecular weights
#            standard_deviation_1 = numpy.sqrt(rescaled_covariance[0][0])
#            standard_deviation_2 = numpy.sqrt(rescaled_covariance[1][1])
#            pgui('standard_deviation_1, standard_deviation_2 (kDa): ' +
#                      str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
#            avars.outfile.write('standard_deviation_1, standard_deviation_2 (kDa): ' +
#                      str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
#            correlation_coefficient = rescaled_covariance[0][1]/(standard_deviation_1*standard_deviation_2)
#            pgui('correlation coefficient: '+str(correlation_coefficient)+'\n')
#            avars.outfile.write('correlation coefficient: '+str(correlation_coefficient)+'\n')

        total_molecular_weight = molecular_weight_1 + molecular_weight_2
        weight_fraction_1 = molecular_weight_1/total_molecular_weight
        weight_fraction_2 = molecular_weight_2/total_molecular_weight

        # Do we want to calculate volume fractions here as well? Will they be needed?

        pgui('Total Mw (kDa), weight_fraction_1, weight_fraction_2: '+str(total_molecular_weight)+'\t'+str(weight_fraction_1)+'\t'+str(weight_fraction_2))
        mcavars.outfile.write('Total Mw (kDa), weight_fraction_1, weight_fraction_2: '+str(total_molecular_weight)+'\t'+str(weight_fraction_1)+'\t'+str(weight_fraction_2))
        mcavars.outfile.close()

        return


    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        st = ''.join(['=' for x in xrange(60)])
        pgui("\n%s \n\n" % (st))        
        pgui('MULTI-COMPONENT ANALYSIS IS DONE')

        time.sleep(1.0)

        return

