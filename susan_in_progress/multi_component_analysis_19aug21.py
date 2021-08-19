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

#import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
#import sassie.multi_component_analysis.stoichiometry as stoichiometry
#import stoichiometry as stoichiometry
import sassie.interface.input_filter_sasmol as input_filter
#import sassie.contrast.muli_component_analysis.read_contrast_output_files as read_contrast_output_files
import read_contrast_output_files as read_contrast_output_files

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

class multi_component_analysis_input_variables():

    def __init__(self, parent=None):
        pass

class multi_component_analysis():

    def __init__(self, parent=None):
        pass

#how is this going to work?  The variables depend on which method is chosen.  Will we be listing all of them and then only using the ones that we need?  In any case, will delta_rho need default values even if reading from a file?  How is this handled currently when we don't use advanced options, for instance?
    def main(self, input_variables, stoichiometry_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = multi_component_analysis_input_variables() #follows SasCalc convention; not sure why they are called "input" variables since they are defined within the method and not input from the user

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

#I'm envisioning one process_input_variables method that will handle the variables from all methods inside of it.
        self.process_input_variables(stoichiometry_variables)

        self.run_utils.general_setup(self)

#The variables processed from here will depend on the method being used, i.e, the chosen flag. How are we going to do this?  With an if, elif ...?

        if self.mvars.stoichiometry_flag:

            self.stoichiometry_initialization()

            self.get_molecular_weights()
#OR
#Do we need one initialization and a call to one executable that will call multiple methods?
        self.epilogue()

        return



    def unpack_variables(self,variables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.output_file_name = variables['output_file_name'][0]
        mvars.input_file_name = variables['input_file_name'][0]
        mvars.read_from_file = variables['read_from_file'][0]
        mvars.number_of_contrast_points = variables['number_of_contrast_points'][0]
        mvars.number_of_components = variables['number_of_components'][0]
        mvars.stoichiometry_flag = variables['stoichiometry_flag'][0]
        mvars.match_point_flag = variables['match_point_flag'][0] 
        mvars.stuhrman_parallel_axis_flag = variables['stuhrman_parallel_axis_flag'][0]
        mvars.decomposition_flag = variables['decomposition_flag'][0]

        log.debug(vars(mvars))

        return


#Will this be a separate method only for stoichiometry variables? Or will ariables for other methods be added later?
    def process_input_variables(self,stoichiometry_variables):

        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui

        log.debug('in process_input_variables')

        log.debug('stoichiometry_variables: %s' %(str(stoichiometry_variables)))
        log.debug('len(stoichiometry_variables): %i' % len(stoichiometry_variables))
        log.debug('number of contrast points: %i' % (mvars.number_of_contrast_points))
        log.debug('number of components: %i' % (mvars.number_of_components))
        log.debug('read from file: %i' % (mvars.read_from_file))        

        mvars.fraction_d2o = stoichiometry_variables[0]
        mvars.partial_specific_volume = stoichiometry_variables[1]
        mvars.izero = stoichiometry_variables[2]
        mvars.concentration = stoichiometry_variables[3]
        if mvars.read_from_file == 1:
#where is this input file?  The gui mimic input assumes it is in the current path.
            delta_rho_read_from_file = read_contrast_output_files.read_contrast_file(mvars.input_file_name, mvars.fraction_d2o, mvars.number_of_contrast_points, mvars.number_of_components)
#            delta_rho_read_from_file = read_contrast_output_files.read_contrast_file()
#remove this print statement
            print('delta rho after reading from file: ', delta_rho_read_from_file)
            log.debug('delta rho after reading from file: %s' % (str(delta_rho_read_from_file)))
            mvars.delta_rho = delta_rho_read_from_file
        elif mvars.read_from_file == 0:
            mvars.delta_rho = stoichiometry_variables[4]
        else:
            message = 'read_from_file is not defined'
            message += ': stopping here'
            pgui(message)
#TODO:  test read_from_file variable in the input filter rather than here.  This is just for diagnostics

        log.debug(vars(mvars))

        return




#Will there be multiple initialization methods depending on which analysis is being done?

    def stoichiometry_initialization(self):
        '''
        method to prepare for get_molecular_weights
        '''

        log = self.log
        log.debug('in stoichiometry_initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        pgui('setting up output file path')
        
#if doing one "big" intialization, we need to ask which method is being initialized to put a sub-path for the method used, i.e., multi_component_analysis/stoichiometry
        if (mvars.runname[-1] == '/'):
            log.debug('runname(1) = %s' % (mvars.runname))
            multi_component_analysis_path = mvars.runname + 'multi_component_analysis/stoichiometry/'
            log.debug('multi_component_analysis_path = %s' % (multi_component_analysis_path))
        else:
            log.debug('runname(2) = %s' % (mvars.runname))
            multi_component_analysis_path = mvars.runname + '/multi_component_analysis/stoichiometry/'
            log.debug('multi_component_analysis_path = %s' % (multi_component_analysis_path))

        direxist = os.path.exists(multi_component_analysis_path)
        if(direxist == 0):
            os.system('mkdir -p ' + multi_component_analysis_path)


        # this is just to write the contrast calculator filename into the output file so we know this option was used
        avars.outfile = io.open(multi_component_analysis_path+mvars.output_file_name, 'w')
        if mvars.read_from_file == 1:
            pgui('input file: %s' % (mvars.input_file_name))
            avars.outfile.write('input file: ' + mvars.input_file_name +'\n')
        else:
            pgui('input file: None')
            avars.outfile.write('input file: None\n')

        log.debug(vars(mvars))
        log.debug(vars(avars))


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

        mvars = self.mvars
        avars = self.avars

        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        Na = 6.023  #10**23

        partial_specific_volume_1 = locale.atof(mvars.partial_specific_volume[0])
        partial_specific_volume_2 = locale.atof(mvars.partial_specific_volume[1])
        print ('partial_specific_volume_1,partial_specific_volume_2: ', partial_specific_volume_1, partial_specific_volume_2)

        # print ("number of contrast points: ", number_of_contrast_points)

        pgui('setting up I(0) equation coefficients')
        
        x = numpy.zeros(shape=(mvars.number_of_contrast_points, 4))

    # Here we are defining the x coeficients at each contrast for the set of simultaneous I(0) equations
    # NOTE that range is being used here since xrange in Python 2 == range in Python 3. So, we aren't doing exactly the same thing here as will be done in Python 3, as range will behave like Python 2's xrange in Python 3.  But, for this case, it shouldn't matter (if I understand correctly).
        for i in range(mvars.number_of_contrast_points):
            delta_rho_1 = locale.atof(mvars.delta_rho[i][0])
            delta_rho_2 = locale.atof(mvars.delta_rho[i][1])
            izero_1 = locale.atof(mvars.izero[i])
            concentration_1 = locale.atof(mvars.concentration[i])
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
        
        # optimized_molecular_weights are the optimal values for the molecular weights so that the sum of the squared error of izero_function(x, *optimized_molecular_weights) - y is minimized
        # covariance is the estimated covariance of molecular weights. The diagonals are used to obtained the standard deviations (below).
        optimized_molecular_weights, covariance = scipy.optimize.curve_fit(izero_function, x, y)

        print ('calculated molecular_weight_1,molecular_weight_2: ', optimized_molecular_weights)
        avars.outfile.write('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
        print ('calculated covariance: ', covariance)
        avars.outfile.write('calculated covariance: '+str(covariance)+'\n')

        print ('-------------------------------')
        avars.outfile.write('--------------------------------\n')
        print ('Final Results:')
        avars.outfile.write('Final Results\n')

        # convert molecular weights to kDa
        molecular_weight_1 = optimized_molecular_weights[0]*10**3
        molecular_weight_2 = optimized_molecular_weights[1]*10**3
        print ('molecular_weight_1, molecular_weight_2 (kDa): ', molecular_weight_1, molecular_weight_2)
        avars.outfile.write('molecular_weight_1, molecular_weight_2 (kDa): '+str(molecular_weight_1)+'\t'+str(molecular_weight_2)+'\n')

        # The covariance matrix will be infinite if there are only two contrasts, so we don't want to calculate it in this case; numpy.all tests whether all array elements along a given axis evaluate to True.
        if numpy.all(covariance != numpy.inf):

            # rescale values in covariance matrix since molecular weights were converted to kDa
            rescaled_covariance = covariance*10**6
            print ('covariance matrix: ', rescaled_covariance)
            avars.outfile.write('covariance matrix: ' + str(rescaled_covariance)+'\n')
            # compute one standard deviation of the molecular weights
            standard_deviation_1 = numpy.sqrt(rescaled_covariance[0][0])
            standard_deviation_2 = numpy.sqrt(rescaled_covariance[1][1])
            print ('standard deviations, standard_deviation_1, standard_deviation_2: ', standard_deviation_1, standard_deviation_2)
            avars.outfile.write('standard deviations, standard_deviation_1, standard_deviation_2: ' +
                      str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
            correlation = rescaled_covariance[0][1]/(standard_deviation_1*standard_deviation_2)
            print ('correlation: ', correlation)
            avars.outfile.write('correlation: '+str(correlation)+'\n')

        total_molecular_weight = molecular_weight_1 + molecular_weight_2
        weight_fraction_1 = molecular_weight_1/total_molecular_weight
        weight_fraction_2 = molecular_weight_2/total_molecular_weight

        # Do we want to output the covariance matrix? How much meaning do the standard deviations have in this case since we are performing an "unweighted" fit? We aren't specifying errors on x, as they are hard to determine. While errors on I(0) and concentration can be estimated, it is harder to estimate the errors on delta rho and partial specific volume. 
        # Do we want to calculate volume fractions here as well? Will they be needed?
        # TODO: We will want to write all of the input values to the output file so that we know what values were used in the calculation.

        print ('Mw (kDa), weight_fraction_1, weight_fraction_2: ', total_molecular_weight, weight_fraction_1, weight_fraction_2)
        avars.outfile.write('Mw (kDa), weight_fraction_1, weight_fraction_2: '+str(total_molecular_weight)+'\t'+str(weight_fraction_1)+'\t'+str(weight_fraction_2))
        avars.outfile.close()

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

