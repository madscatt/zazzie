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
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#       MULTI COMPONENT ANALYSIS: GET MATCH POINT
#
#       10/22/2021       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''

    MATCH POINT contains the GET MATCH POINT method that calculates the match point from SANS contrast variation data.
    
    INPUTS: 
        I(0), concentration, fraction D2O at each contrast
    
    OUTPUTS:
        contrast match point
    
    Called by the MULTI-COMPONENT ANALYSIS module.
    
'''


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import io
import os
import sys
import string
import locale
import time
import numpy
import linear_fit as linear_fit


def get_match_point(other_self):
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

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    mvars = other_self.module_variables

    log.debug('in get match point')    

#mvars used:  fraction_d2o, number_of_contrast_points, izero, izero_error, concentration, concentration_error, initial_match_point_guess 
#mcavars used:  outfile


    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables
        
    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in xrange(60)])

    pgui('\n%s \n' % (st))
    pgui('DATA FROM RUN: %s \n\n' % (ttxt))

#   In most cases, this routine will only be called when analyzing data (since planning is done with the contrast calculator), so it is assumed that the user has input I(0), c and an error on both I(0) and c, and a weighted fit will be used.  Values of "0" for the errors aren't allowed, so this is checked in the module filter.

# calculate sqrt[I(0)/c] and propagate the error on I(0)/c to make the relationship between I(0) and fraction D2O linear.  The user needs to to input a frac_d2o value for which the sign changes, i.e, an initial guess for the match point. Several values may need to be tried to get an optimal fit.  An initial guess should be able to be made from the data (perhaps with help from contrast calculator output).
#TODO:  consider doing a 2nd order polynomial fit to I(0)/c vs fraction_d2o to get an estimate for the match point instead of having the user provide one. (This is done in MulCh, for instance.)

    mode = 1 #for weighted fit
    normalized_izero = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
    normalized_izero_error = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
    square_root_izero = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
    square_root_izero_error = numpy.zeros(mvars.number_of_contrast_points, numpy.float)
 
#    print('fraction d2o, izero, error: ', fraction_d2o, izero, izero_error)

    pgui('Match Point Method\n')
    pgui('calculating sqrt[I(0)/c] and propagating errors')        

#TODO:  should concentration be converted to g/cm^3?  See what is done in MulCh and contrast calculator. 

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

