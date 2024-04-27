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
#       5/9/2023         --      Python 3               :   Susan Krueger
#       7/26/2023        --      added polynomial fit   :   Susan Krueger
#       3/28/2024        --      revised polynomial fit :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''

    **Match Point** contains the **Get Match Point** method that calculates the match point from SANS contrast variation data.

    **Inputs:**

    I(0), concentration, fraction D\ :sub:`2`\ O at each contrast

    **Outputs:**

    contrast match point

    Called by the **Multi-component Analysis** module.

'''

import io
import os
import time
import numpy
import scipy.optimize
import sassie.contrast.multi_component_analysis.polynomial_fit as polynomial_fit
import sassie.contrast.multi_component_analysis.chi_squared_correlation as chi_squared_correlation
import json

#import polynomial_fit as polynomial_fit
#import chi_squared_correlation as chi_squared_correlation

def save_data_to_plot_as_json(other_self, square_root_izero, square_root_izero_error, square_root_izero_calculated, diff):

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    data_dict = {
        'fraction_d2o': [],
        'sqrt[I(0)/c]': [],
        'sqrt[I(0)/c]_error': [],
        'sqrt[I(0)/c]_calc': [],
        'sqrt[I(0)/c]-sqrt[I(0)/c]_calc': []
    }

    for i in range(mvars.number_of_contrast_points):
    # Append new values to the lists in the dictionary
        data_dict['fraction_d2o'].append(mvars.fraction_d2o[i])
        data_dict['sqrt[I(0)/c]'].append(square_root_izero[i])
        data_dict['sqrt[I(0)/c]_error'].append(square_root_izero_error[i])
        data_dict['sqrt[I(0)/c]_calc'].append(square_root_izero_calculated[i])
        data_dict['sqrt[I(0)/c]-sqrt[I(0)/c]_calc'].append(diff[i])
        
    json_data = json.dumps(data_dict)

    mcavars.json_outfile.write(json_data)
    mcavars.json_outfile.close()

    return


def get_match_point(other_self):
    '''
    **Get Match Point** is the Match Point Analysis method that calculates the match point from SANS contrast variation data. I(0) is divided by the concentration (c) to normalize for differences in concentration at each contrast. To make the relationship between I(0)/c and fraction D2O linear, sqrt[I(0)/c] is calculated and the sign in changed for the D\ :sub:`2`\ O values that are less than the match point where I(0) ~ 0. Thus, an initial guess for the match point is determined by first performing a 2nd order polynomial fit to I(0)/c vs fraction D\ :sub:`2`\ O and then calculating the fraction D\ :sub:`2`\ O value where the derivative of the polynomial function is zero.  If this fails to provide a reasonable initial match point guess, the user can input a value as an advanced option.

    This method is expected to be called when analyzing I(0) vs fraction D\ :sub:`2`\ O data. A weighted linear fit is performed to sqrt[I(0)/c] vs fraction D\ :sub:`2`\ O in order to obtain the match point from the slope and intercept of the fitted curve using the relation:  match point = -intercept/slope.

    Notes:

        Output file with chosen output_file_name is stored in the output file directory.
        TODO: Need boilerplate language to add variable name/path for the directory.

        Output file contains:
            - module variable input parameters
            - number of points used in the linear fit
            - calculated match point and error
            - R value and reduced chi-squared for the linear fit
            - table containing:
                - fraction D\ :sub:`2`\ O
                - sqrt[I(0)/c]
                - sqrt[I(0)/c] error
                - calculated sqrt[I(0)/c]
                - sqrt[I(0)/c] - calculated sqrt[I(0)/c]

    Requires **Polynomial Function**, **Get Reduced Chi Squared**, **Correlation from Covariance** and **scipy.optimize.curve_fit**.

    Parameters
    ----------
    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    fraction_d2o:   float array (dimension = number_of_contrast_points)
        The fraction D\ :sub:`2`\ O values that define the contrasts
    izero:  float array (dimension = number_of_contrast_points)
        I(0) value at each contrast in cm\ :sup:`-1`\
    izero_error:  float array (dimension = number_of_contrast_points)
        I(0) error value at each contrast
    concentration:  float array (dimension = number_of_contrast_points)
        concentration at each contrast in mg/mL
    concentration_error:  float array (dimension = number_of_contrast_points)
        concentration error at each contrast
    initial_match_point_guess_flag:  boolean (default = False)
        advanced option flag to allow the initial match point guess to be input directly rather than calculated from an 2nd order polynomial fit to I(0)/c vs fraction D\ :sub:`2`\ O 
    initial_match_point_guess:  float (default = 0.5)
        initial match point guess that is read if intial_match_point_guess_flag = True
    outfile: string
        output file name (with full path): path + output_file_name

    Returns
    -------
    match_point: float
        The fraction D\ :sub:`2`\ O value corresponding to the contrast match point
    match_point_error: float
        The error on the contrast match point

    '''

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    mvars = other_self.module_variables

    log.debug('in get match point')

# mvars used:  fraction_d2o, number_of_contrast_points, izero, izero_error, concentration, concentration_error, initial_matchpoint_guess_flag, initial_matchpoint_guess
# mcavars used:  outfile

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in range(60)])

    pgui('\n%s \n' % (st))
    pgui('DATA FROM RUN: %s \n\n' % (ttxt))

#   In most cases, this routine will only be called when analyzing data (since planning is done with the contrast calculator), so it is assumed that the user has input I(0), c and an error on both I(0) and c, and a weighted fit will be used.  Values of "0" for the errors aren't allowed, so this is handled in the module filter.

# calculate the normalized I(0)/c and propagate the errors to get I(0)/c error
# calculate sqrt[I(0)/c] and propagate the errors to get sqrt[I(0)/c] error. To make the relationship between I(0) and fraction D2O linear the sign of sqrt[I(0)/c] is changed for all frac_d2o values that are greater than the match point. An initial guess for the match point is made by fitting a 2nd order polynomial to I(0)/c vs frac_d2o.

    normalized_izero = numpy.zeros(
        mvars.number_of_contrast_points, numpy.float)
    normalized_izero_error = numpy.zeros(
        mvars.number_of_contrast_points, numpy.float)
    square_root_izero = numpy.zeros(
        mvars.number_of_contrast_points, numpy.float)
    square_root_izero_error = numpy.zeros(
        mvars.number_of_contrast_points, numpy.float)

#    print('fraction d2o, izero, error: ', fraction_d2o, izero, izero_error)

    pgui('Match Point Method\n')
    pgui('calculating I(0)/c and propagating errors')

# TODO:  should concentration be converted to g/cm^3?  See what is done in MulCh and contrast calculator.

    for i in range(mvars.number_of_contrast_points):
        normalized_izero[i] = mvars.izero[i]/mvars.concentration[i]
        normalized_izero_error[i] = numpy.sqrt(mvars.izero_error[i]**2/mvars.concentration[i]
                                               ** 2 + mvars.izero[i]**2*mvars.concentration_error[i]**2/mvars.concentration[i]**4)

    log.debug('I(0)/c, error: ' + str(normalized_izero) +
              ',' + str(normalized_izero_error) + '\n')

# polynomial fit to izero = a2 + a1*fraction_d2o + a0*fraction_d2o**2; note that a0 is the coefficient for the quadratic term and a2 is the constant term; this seems backwards to me, but that is the way the polynomial fit function was written

    log.debug('initial match point guess flag: ' +
              str(mvars.initial_match_point_guess_flag) + '\n')
    if (mvars.initial_match_point_guess_flag == False):
        pgui('performing polynomial fit to determine initial match point guess')
# polynomial fit to izero = a0 + a1*fraction_d2o + a2*fraction_d2o**2
        order_polynomial_fit = 2
        initial_guess_polynomial_fit = [1.]*(order_polynomial_fit + 1)
        log.debug('initial guess for polynomial coefficients: ' +
                  str(initial_guess_polynomial_fit) + '\n')
        polynomial_fit_coefficients, polynomial_fit_covariance = scipy.optimize.curve_fit(
            polynomial_fit.polynomial_function, mvars.fraction_d2o, normalized_izero, sigma=normalized_izero_error, absolute_sigma=True, p0=initial_guess_polynomial_fit)
# for debugging
        log.debug('polynomial fit coefficients: ' +
                  str(polynomial_fit_coefficients) + '\n')
        log.debug('polynomial fit covariance matrix: ' +
                  str(polynomial_fit_covariance) + '\n')
        normalized_izero_calculated = polynomial_fit.polynomial_function(
            mvars.fraction_d2o, *polynomial_fit_coefficients)
        log.debug('I(0)/c calculated: ' +
                  str(normalized_izero_calculated) + '\n')
        polynomial_fit_reduced_chi_squared = chi_squared_correlation.get_reduced_chi_squared(
            normalized_izero, normalized_izero_error, normalized_izero_calculated, mvars.number_of_contrast_points, 3)
        log.debug('polynomial fit reduced chi-squared: ' +
                  str(polynomial_fit_reduced_chi_squared) + '\n')
        # these numbers agree with a weighted polynomial fit in IGOR
        polynomial_fit_error = numpy.sqrt(
            numpy.diag(polynomial_fit_covariance))
        log.debug('polynomial fit standard deviation: ' +
                  str(polynomial_fit_error) + '\n')
# end debugging

# find where derivative = 0 and set to initial match point guess
# 0 = a1 + 2*a2*match_point_guess; match_point_guess = -a1/(2*a2)

        match_point_guess = - \
            polynomial_fit_coefficients[1]/(2*polynomial_fit_coefficients[2])
        pgui('initial match point guess from polynomial fit: ' +
             str(round(match_point_guess, 4)) + '\n')
    else:
        match_point_guess = mvars.initial_match_point_guess
        pgui('initial match point guess: ' +
             str(round(match_point_guess, 4)) + '\n')

    pgui('calculating sqrt[I(0)/c] and propagating errors')

# To obtain a straight line, set sqrt(I(0)/c) as negative if fraction D2O is >= to the initial matchpoint guess
    for i in range(mvars.number_of_contrast_points):
        if (mvars.fraction_d2o[i] <= match_point_guess):
            square_root_izero[i] = numpy.sqrt(normalized_izero[i])
        else:
            square_root_izero[i] = -numpy.sqrt(normalized_izero[i])
        square_root_izero_error[i] = normalized_izero_error[i] / \
            (2.0*numpy.sqrt(normalized_izero[i]))

    log.debug('sqrt[I(0)/c], error: ' + str(square_root_izero) +
              ',' + str(square_root_izero_error) + '\n')

# calculate the linear fit of sqrt[I(0)] vs fraction D2O

    order_line_fit = 1
    initial_guess_line_fit = [1.]*(order_line_fit + 1)
    line_fit, line_fit_covariance = scipy.optimize.curve_fit(
        polynomial_fit.polynomial_function, mvars.fraction_d2o, square_root_izero, sigma=square_root_izero_error, absolute_sigma=True, p0=initial_guess_line_fit)

    pgui('performing linear fit to calculate match point')
    log.debug('line fit coefficients: ' + str(line_fit) + '\n')
    log.debug('line fit covariance matrix: ' + str(line_fit_covariance) + '\n')
    square_root_izero_calculated = polynomial_fit.polynomial_function(
        mvars.fraction_d2o, *line_fit)
    log.debug('sqrt[I(0)/c calculated]: ' +
              str(square_root_izero_calculated) + '\n')
    # for writing to file and plotting
    diff = square_root_izero - square_root_izero_calculated
    log.debug('diff: ' + str(diff) + '\n')

    line_fit_reduced_chi_squared = chi_squared_correlation.get_reduced_chi_squared(
        square_root_izero, square_root_izero_error, square_root_izero_calculated, mvars.number_of_contrast_points, 2)

    log.debug('line fit reduced chi-squared: ' +
              str(line_fit_reduced_chi_squared) + '\n')

# NOTE: To get the correlation matrix and the std deviation (line_fit_error):
    line_fit_error, line_fit_correlation = chi_squared_correlation.correlation_from_covariance(
        line_fit_covariance)
    # these numbers agree with a weighted polynomial fit in IGOR
    log.debug('line fit std deviation: ' + str(line_fit_error) + '\n')
    log.debug('line fit correlation matrix: ' +
              str(line_fit_correlation) + '\n')
    correlation_coefficient = line_fit_covariance[0][1] / (
        line_fit_error[0]*line_fit_error[1])
    log.debug('line fit correlation coefficient: ' +
              str(line_fit_correlation[0][1]) + '\n')

# calculate the match point from the slope and intercept
    intercept = line_fit[0]
    slope = line_fit[1]
    intercept_error = line_fit_error[0]
    slope_error = line_fit_error[1]
    match_point = -intercept/slope
    match_point_error = numpy.sqrt(
        (intercept_error/slope)**2 + ((slope_error*intercept)/(slope**2))**2)
    log.debug('match point, error: ' + str(match_point) +
              ',' + str(match_point_error) + '\n')

    pgui('results written to output file: %s' %
         (mcavars.multi_component_analysis_path+mvars.output_file_name))
    pgui('-------------------------------')
    mcavars.outfile.write('--------------------------------\n')
    pgui('Final Results\n')
    mcavars.outfile.write('Final Results\n')
    pgui('number of points fit: ' + str(mvars.number_of_contrast_points) + '\n')
    mcavars.outfile.write('number of points fit: ' +
                          str(mvars.number_of_contrast_points) + '\n')
    if (mvars.initial_match_point_guess_flag == False):
        pgui('initial match point guess from polynomial fit: ' +
             str(round(match_point_guess, 4)) + '\n')
        mcavars.outfile.write('initial match point guess from polynomial fit: ' +
                              str(round(match_point_guess, 4)) + '\n')
    else:
        pgui('initial match point guess: ' +
             str(round(match_point_guess, 4)) + '\n')
        mcavars.outfile.write('initial match point guess: ' +
                              str(round(match_point_guess, 4)) + '\n')
    pgui('match point: ' + str(round(match_point, 4)) +
         ' +/- ' + str(round(match_point_error, 4)) + '\n')
    mcavars.outfile.write('match point: ' + str(round(match_point, 4)) +
                          ' +/- ' + str(round(match_point_error, 4)) + '\n')
# TODO:  check for reduced chi-squared = -1 and print out N/A message instead of value; could occur if user only has 2 contrast points, i.e., H2O and D2O; this is not typical for a contrast variation study and the polynomial fit might fail, so user would want to enter an initial match point guess
    mcavars.outfile.write('reduced chi-squared for linear fit: ' +
                          str(round(line_fit_reduced_chi_squared, 4)) + '\n')
    mcavars.outfile.write(
        'fraction_d2o  sqrt[I(0)/c]  sqrt[I(0)/c]_error  sqrt[I(0)/c]_calc  sqrt[I(0)/c]-sqrt[I(0)/c]_calc\n')
    for i in range(mvars.number_of_contrast_points):
        mcavars.outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\t\t%9.4f\n' % (
            mvars.fraction_d2o[i], square_root_izero[i], square_root_izero_error[i], square_root_izero_calculated[i], diff[i]))
    mcavars.outfile.close()

# TODO: plot sqrt[I(0)/c] with errorbars and sqrt[I(0)/c]_calc vs fraction D2O and the residuals, i.e., diff = sqrt[I(0)/c]-sqrt[I(0)/c]_calc, with residuals on a separate plot.  The values will be read from mcavars.outfile.
    
    save_data_to_plot_as_json(other_self, square_root_izero, square_root_izero_error, square_root_izero_calculated, diff)


    return
