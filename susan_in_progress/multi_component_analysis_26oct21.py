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
    MULTI-COMPONENT ANALYSIS is the module that contains the methods
    that are used to analyze data from contrast variation experiments.
    It can be used on model data prior to an experiment or on measured data.
    
    Included are 4 independent methods (give file names if they end up in other files): 
        * Match Point Analysis
        * Stuhrmann and Parallel Axis Theorem Analyses
        * Decomposition Analysis
        * Stoichiometry Analysis

    This method requires read_contrast_output_files (and any more files that will eventually be called; perhaps put them under the individual method names unless a required file is common to all).

    INPUTS:
        The inputs depend on which method is used. The inputs are noted in the independent methods. 
        The method is chosen based on which of the following flags is True:
            match_point_flag
            stuhrmann_parallel_axis_flag
            decomposition_flag
            stoichiometry_flag
        Only one flag at a time can be True.
    
    OUTPUTS:
            The outputs depend on which method is used. The names of output files and their contents are noted in the independent methods.

    (Need to encapsulate the flow chart, use cases, etc., flow that is hardcoded with if statements, for example.)
    
    (Need boilerplate line that shows the flow, i.e.,  module utilities, setup logging, unpack variables, run main, etc. This will be in all modules.)

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
#import sassie.contrast.multi_component_analysis.read_contrast_output_files as read_contrast_output_files
import read_contrast_output_files as read_contrast_output_files
import linear_fit as linear_fit

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

    def main(self, input_variables,txtOutput):

        self.module_variables = module_variables()

        self.multi_component_analysis_variables = multi_component_analysis_variables() 

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)  

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.multi_component_analysis_initialization()

        if self.module_variables.match_point_flag:

            self.get_match_point()

        elif self.module_variables.stoichiometry_flag:

            self.get_molecular_weights()

        self.epilogue()

        return


    def unpack_variables(self,variables):
        '''
        method to unpack variables passed from the GUI
        '''

        log = self.log
        mvars = self.module_variables
        log.debug('in unpack_variables')

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

#unpack the variables that are unique to each method. If some variables turn out to be needed for all methods, then they can be unpacked above before checking the flags. The read_from_file if statement may not be needed if the GUI is going to be populated with these values prior to running the main program.
        if mvars.stoichiometry_flag == True:
            mvars.izero = variables['izero'][0]
            mvars.concentration = variables['concentration'][0]
            mvars.partial_specific_volume = variables['partial_specific_volume'][0]

            if mvars.read_from_file == True:
                read_contrast_output_files.read_contrast_file(self)
#                print('delta rho after reading from file: ', mvars.delta_rho)
                log.debug('delta rho after reading from file: %s' % (str(mvars.delta_rho)))
#               mvars.delta_rho = delta_rho_read_from_file
            elif mvars.read_from_file == False:
                mvars.delta_rho = variables['delta_rho'][0]
                
        elif mvars.match_point_flag == True:
            mvars.izero = variables['izero'][0]
            mvars.izero_error = variables['izero_error'][0]
            mvars.concentration = variables['concentration'][0]
            mvars.concentration_error = variables['concentration_error'][0]
            mvars.initial_match_point_guess = variables['initial_match_point_guess'][0]

#        print(vars(mvars))

        log.debug(vars(mvars))

        return


#Will there be multiple initialization methods depending on which analysis is being done?  If so, the initialization and main method for each analysis should be in a separate file (formatted like read_contrast_output_files).

    def multi_component_analysis_initialization(self):
        '''
        Method to initialize multi_component_analysis variables and 
        write initial information to output files.
        
        INPUTS:
            module variables
                stoichiometry_flag
                match_point_flag
                stuhrmann_parallel_axis_flag 
                decomposition_flag
                run_name 
                read_from_file   
        
        OUTPUTS:
            variables added to multi_component analysis variables
                multi_component_analysis_path
                outfile        
        '''

        log = self.log
        log.debug('in multi_component_analysis_initialization')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        mcavars = self.multi_component_analysis_variables

        
# Need to ask which method is being initialized to put a sub-path for the method used, i.e., multi_component_analysis/method.
        if mvars.stoichiometry_flag == True:
            if (mvars.run_name[-1] == '/'):
                log.debug('run_name(1) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + 'multi_component_analysis/stoichiometry/'
                log.debug('multi_component_analysis_path = %s' % (mcavars.multi_component_analysis_path))
            else:
                log.debug('run_name(2) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + '/multi_component_analysis/stoichiometry/'
                log.debug('multi_component_analysis_path = %s' % (mcavars.multi_component_analysis_path))
        elif mvars.match_point_flag == True:
            if (mvars.run_name[-1] == '/'):
                log.debug('run_name(1) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + 'multi_component_analysis/match_point/'
                log.debug('multi_component_analysis_path = %s' % (mcavars.multi_component_analysis_path))
            else:
                log.debug('run_name(2) = %s' % (mvars.run_name))
                mcavars.multi_component_analysis_path = mvars.run_name + '/multi_component_analysis/match_point/'
                log.debug('multi_component_analysis_path = %s' % (mcavars.multi_component_analysis_path))

        direxist = os.path.exists(mcavars.multi_component_analysis_path)
        if(direxist == 0):
            os.system('mkdir -p ' + mcavars.multi_component_analysis_path)

        # this is just to write the contrast calculator filename into the output file so we know this option was used
        mcavars.outfile = io.open(mcavars.multi_component_analysis_path+mvars.output_file_name, 'w')

        if mvars.read_from_file == True:
#            pgui('input file: %s' % (mvars.input_file_name))
            mcavars.outfile.write('input file: ' + mvars.input_file_name +'\n')
        else:
#            pgui('input file: None')
            mcavars.outfile.write('input file: None\n')
            
#TODO: write out the individual input variables on separate lines and only write out the ones that are relevant for the method being used
        mcavars.outfile.write('input variables: ' + repr(vars(mvars)) + '\n')
 
        log.debug(vars(mvars))
        log.debug(vars(mcavars))


    def get_match_point(self):
        '''
        GET MATCH POINT is the Match Point Analysis method that calculates the match point from SANS contrast variation data. I(0) is divided by the concentration (c) to normalize for differences in concentration at each contrast. To make the relationship between I(0)/c and fraction D2O linear, sqrt[I(0)/c] is calculated and the sign in changed for the fraction D2O values that are less than the match point where I(0) ~ 0. Thus, an initial guess for the match point must be provided. This can usually be guessed from the I(0) vs fraction D2O data.
        This method is expected to be called when analyzing I(0) vs fraction D2O data. A weighted linear fit is performed to sqrt[I(0)/c] vs fraction D2O in order to obtain the match point from the slope and intercept of the fitted curve using the relation:  match point = -intercept/slope.
                
        INPUTS:
            module variables
                run_name:                           name of directory containing the outputs
                number_of_contrast_points:          number of contrasts being used in the molecular weights calculation
                fraction_d2o array:                 the fraction of D2O in the solvent for each contrast
                izero array:                        I(0) value at each contrast
                izero error array:                  I(0) error value at each contrast
                concentration array:                concentration at each contrast
                concentration error array:          concentration error at each contrast
                initial match point guess           fraction of D2O to be used as initial match point guess so that 

            multi_component_analysis_variables
                output_file_name:                   name of output file
        
        OUTPUTS:
            files stored in output file directory (need boilerplate language to add variable name/path for the directory):

                output file with chosen output_file_name containing:
                    module variable input parameters
                    number of points used in the linear fit
                    calculated match point and error 
                    r value and reduced chi-squared for the linear fit
                    table containing:
                        fraction D2O
                        sqrt[I(0)/c]
                        sqrt[I(0)/c] error
                        calculated sqrt[I(0)/c]
                        sqrt[I(0)/c] - calculated sqrt[I(0)/c]

        Requires line_fit
        '''

#mvars used:  fraction_d2o, number_of_contrast_points, izero, izero_error, concentration, concentration_error, initial_match_point_guess 
#mcavars used:  outfile

       
        log = self.log
        log.debug('in get_match_point')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        mcavars = self.multi_component_analysis_variables
        
        log.debug(vars(mvars))
        log.debug(vars(mcavars))

        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui('\n%s \n' % (st))
        pgui('DATA FROM RUN: %s \n\n' % (ttxt))

#   In most cases, this routine will only be called when analyzing data (since planning is done with the contrast calculator), so it is assumed that the user has input I(0), c and an error on both I(0) and c, and a weighted fit will be used.  Values of "0" for the errors aren't allowed, so this is checked in the module filter.

# calculate sqrt[I(0)/c] and propagate the error on I(0)/c to make the relationship between I(0) and fraction D2O linear.  The user needs to to input a frac_d2o value for which the sign changes, i.e, an initial guess for the match point. Several values may need to be tried to get an optimal fit.  An initial guess should be able to be made from the data (perhaps with help from contrast calculator output).

        mode = 1 #for weighted fit
        normalized_izero = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
        normalized_izero_error = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
        square_root_izero = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
        square_root_izero_error = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
 
#       print('fraction d2o, izero, error: ', fraction_d2o, izero, izero_error)

        pgui('Match Point Method\n')
        pgui('calculating sqrt[I(0)/c] and propagating errors')        

        for i in range(mvars.number_of_contrast_points):
            normalized_izero[i] = mvars.izero[i]/mvars.concentration[i]
            normalized_izero_error[i] = numpy.sqrt(mvars.izero_error[i]**2/mvars.concentration[i]**2 + mvars.izero[i]**2*mvars.concentration_error[i]**2/mvars.concentration[i]**4)
#To obtain a straight line, set sqrt(I(0)/c) as negative if fraction D2O is <= to the initial matchpoint guess
            if (mvars.fraction_d2o[i] <= mvars.initial_match_point_guess):
                square_root_izero[i] = numpy.sqrt(normalized_izero[i])
            else:
                square_root_izero[i] = -numpy.sqrt(normalized_izero[i])
            square_root_izero_error[i] = normalized_izero_error[i]/(2.0*numpy.sqrt(normalized_izero[i]))

        log.debug('I(0)/c, error: ' + str(normalized_izero) + ',' + str(normalized_izero_error) + '\n')
        log.debug('sqrt[I(0)/c], error: ' + str(square_root_izero) + ',' + str(square_root_izero_error) + '\n')


# calculate the linear fit of sqrt[I(0)] vs fraction D2O

        slope, slope_error, intercept, intercept_error, r_value, reduced_chi_squared, ycalc, diff = linear_fit.line_fit(mvars.fraction_d2o,square_root_izero, square_root_izero_error,mvars.number_of_contrast_points, mode)

        log.debug('slope, error, intercept, error, r_value, reduced_chi_squared: ' + str(slope) +','+ str(slope_error) +','+ str(intercept) +','+ str(intercept_error) +','+ str(r_value) +','+ str(reduced_chi_squared) +'\n')

# calculate the match point from the slope and intercept

        match_point = -intercept/slope
        match_point_error = numpy.sqrt((intercept_error/slope)**2 + ((slope_error*intercept)/(slope**2))**2)


        pgui('results written to output file: %s' % (mcavars.multi_component_analysis_path+mvars.output_file_name))
        pgui('-------------------------------')
        mcavars.outfile.write('--------------------------------\n')
        pgui('Final Results\n')
        mcavars.outfile.write('Final Results\n')
        pgui('number of points fit: ' + str(mvars.number_of_contrast_points) + '\n')
        mcavars.outfile.write('number of points fit: ' + str(mvars.number_of_contrast_points) + '\n')
        pgui('match point: ' + str(round(match_point,4)) + ' +/- ' + str(round(match_point_error,4)) + '\n')
        mcavars.outfile.write('match point: ' + str(round(match_point,4)) + ' +/- ' + str(round(match_point_error,4)) + '\n')
        pgui('r value: ' + str(round(r_value,4)) + '\n')
        mcavars.outfile.write('r value: ' + str(round(r_value,4)) + '\n')
        pgui('reduced chi-squared: ' + str(round(reduced_chi_squared,4)) + '\n')    
        mcavars.outfile.write('reduced chi-squared: ' + str(round(reduced_chi_squared,4)) + '\n')
        mcavars.outfile.write('fraction_d2o  sqrt[I(0)/c]  sqrt[I(0)/c]_error  sqrt[I(0)/c]_calc  sqrt[I(0)/c]-sqrt[I(0)/c]_calc\n')
        for i in range(mvars.number_of_contrast_points):
            mcavars.outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\t\t%9.4f\n' % (mvars.fraction_d2o[i], square_root_izero[i], square_root_izero[i], ycalc[i], diff[i]))   
        mcavars.outfile.close()

#TODO: plot sqrt[I(0)/c] with errorbars vs fraction D2O and the residuals, i.e., sqrt[I(0)/c]-sqrt[I(0)/c]_calc, with residuals on a separate plot

        return


    def izero_function(self, x, molecular_weight_1, molecular_weight_2):
    #TODO: generalize this function any number of components? Or offer hardwired options for 2 and 3 components?
        r'''
        
        IZERO FUNCTION is the I(0) equation 
            
        n(\sum_{i} \Delta  \rho_{i}V_{i})^{2} - I(0) = 0, where 
            
        I(0) is the intensity at q=0
        \Delta \rho_{i} is the contrast of i^{th} component
        V_{i} is the volume of i^{th} component
        n = c\frac{N_{A}}{\sum_{i} (M_{w})_{i}} is the number density of particles where
        c is the concentration
        N_{A} is Avogadro's number and
        \sum_{i} (M_{w})_{i} is the total molecular weight of the complex.
            
        If the volume is rewritten in terms of the partial specific volume, \overline{v},
            
        V = \overline{v} \frac {M_{w}}{N_{A}},
            
        and the sum over M_{w} is expanded, then the equation can be rewritten in terms
        of the individual M_{w} values such that they can be obtained by solving a set
        of simultaneous equations at each contrast.
        
        This function is called from GET MOLECULAR WEIGHTS.
            
        INPUTS:
            module variables
                x array:                array of x coefficients at each contrast (currently hardwired for 2 components)
                    x0:                 \Delta \rho_{1}^{2} \overline{v}_{1}^{2}
                    x1:                 2 \Delta \rho_{1} \Delta \rho_{2} \overline{v}_{1} \overline{v}_{2}
                    x2:                 \Delta \rho_{2}^{2} \overline{v}_{2}^{2}
                    x3:                 \frac(I(0) N_{A}/c)
        
        OUTPUTS: 
            optimized molecular weights, (M_{w})_{1} and (M_{w})_{2}, for the components
                        
        '''
        f = x[0] * molecular_weight_1**2 + x[1] * molecular_weight_1*molecular_weight_2 + x[2] * molecular_weight_2**2 - x[3] * molecular_weight_1 - x[3] * molecular_weight_2
        return f


    def get_molecular_weights(self):
        '''
        GET MOLECULAR WEIGHTS is the Stoichiometry Analysis method that calculates the molecular weights of the components in a complex containing multiple copies of two or more components.
        Executed if stoichiometry_flag == True.
        
        This method calls IZERO FUNCTION.
                       
        INPUTS:

            module variables
                run_name:                           name of directory containing the outputs
                input_file_name:                    name of file containing contrast values calculated from the contrast calculator module
                read_from_file:                     True if contrast values are being read from an input file with chosen input_file_name; False if contrast values are being input by the user
                number_of_contrast_points:          number of contrasts being used in the molecular weights calculation
                number_of_components:               number of components with different contrasts in the complex
                fraction_d2o array:                 the fraction of D2O in the solvent for each contrast
                izero array:                        I(0) value at each contrast
                concentration array:                concentration at each contrast
                partial_specific_volume array:      partial specific volume of each component
                delta_rho matrix:                   contrast value for each component at each fraction of D2O

            multi_component_analysis_variables
                output_file_name:                   name of output file
        
        OUTPUTS:
            files stored in output file directory (need boilerplate language to add variable name/path for the directory):

                output file with chosen output_file_name containing:
                    module variable input parameters
                    calculated molecular weights of the components in kDa
                    total molecular weight of the complex in kDa
                    weight fractions of each component

        '''
#mvars used:  fraction_d2o, number_of_contrast_points, partial_specific_volume, izero, concentration, delta_rho 
#mcavars used:  outfile


       
        log = self.log
        log.debug('in get_molecular_weights')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        mcavars = self.multi_component_analysis_variables
        
        log.debug(vars(mvars))
        log.debug(vars(mcavars))

        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui('\n%s \n' % (st))
        pgui('DATA FROM RUN: %s \n\n' % (ttxt))

        Na = 6.023  #10**23

        partial_specific_volume_1 = mvars.partial_specific_volume[0]
        partial_specific_volume_2 = mvars.partial_specific_volume[1]
#        print ('partial_specific_volume_1,partial_specific_volume_2: ', partial_specific_volume_1, partial_specific_volume_2)

        # print ("number of contrast points: ", number_of_contrast_points)

        pgui('Stoichiometry method\n')
        pgui('setting up I(0) equation coefficients')
        
        if mvars.read_from_file == True:
            #pgui('\n')
            pgui('contrast values read from input file: %s' % (mvars.input_file_name))
        
        x = numpy.zeros(shape=(mvars.number_of_contrast_points, 4))

        # Here we are defining the x coeficients at each contrast for the set of simultaneous I(0) equations.
        # NOTE that range is being used here since xrange in Python 2 == range in Python 3. So, we aren't doing exactly the same thing here as will be done in Python 3, as range will behave like Python 2's xrange in Python 3.  But, for this case, it shouldn't matter (if I understand correctly).
        for i in range(mvars.number_of_contrast_points):
            delta_rho_1 = mvars.delta_rho[i][0]
            delta_rho_2 = mvars.delta_rho[i][1]
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

        # transpose x to have correct inputs for curve_fit, which accepts a 1 column matrix with number of rows = number of contrast points.  We have a 1 row matrix of length = number of contrast points.
        x = x.T  

        pgui('calculating molecular weights')
        

        '''
        scipy.optimize.curve_fit is being used to solve the set of simultaneous equations
        using non-linear least squares fitting.  It uses scipy version 1.2.1 in Python 2.7.
        There is no known change of syntax or usage if using scipy 1.7 under Python 3.
        
        optimized_molecular_weights are the optimal values for the molecular weights such that
        the sum of the squared error of izero_function(x, optimized_molecular_weights) - y
        is minimized. 
        
        The y array = 0 since izero_function is defined to be equal to zero
        The x array = the values of the coefficients that consist of the known parameters in izero_function
        once it is written in terms of the individual M_{w} values.
        
        '''

        optimized_molecular_weights, covariance = scipy.optimize.curve_fit(self.izero_function, x, y)

#TODO: decide what results should be in the output file and what should be printed out to the GUI screen.  Also, Test some errors on the coefficients (hardwired) here to see if we should be allowing users to enter errors. Look at Watson et al power law paper for ranges of values for molecular volume. Would we ask for errors on each parameter, i.e., izero, delta_rho, concentration, etc. and then propagate them here when we calculate the coefficients? If we decide to include errors, then we can report the covariance matrix and standard deviations. 

#        pgui('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
#        mcavars.outfile.write('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
#        pgui('calculated covariance: '+str(covariance)+'\n')
#        mcavars.outfile.write('calculated covariance: '+str(covariance)+'\n')

        pgui('results written to output file: %s' % (mcavars.multi_component_analysis_path+mvars.output_file_name))
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
#            mcavars.outfile.write('covariance matrix: ' + str(rescaled_covariance)+'\n')
#            # compute one standard deviation of the molecular weights
#            standard_deviation_1 = numpy.sqrt(rescaled_covariance[0][0])
#            standard_deviation_2 = numpy.sqrt(rescaled_covariance[1][1])
#            pgui('standard_deviation_1, standard_deviation_2 (kDa): ' +
#                      str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
#            mcavars.outfile.write('standard_deviation_1, standard_deviation_2 (kDa): ' +
#                      str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
#            correlation_coefficient = rescaled_covariance[0][1]/(standard_deviation_1*standard_deviation_2)
#            pgui('correlation coefficient: '+str(correlation_coefficient)+'\n')
#            mcavars.outfile.write('correlation coefficient: '+str(correlation_coefficient)+'\n')

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
        pgui('\n%s \n\n' % (st))        
        pgui('MULTI-COMPONENT ANALYSIS IS DONE')

        time.sleep(1.0)

        return

