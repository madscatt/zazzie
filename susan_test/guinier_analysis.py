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
#       GUINIER ANALYSIS
#
#       4/23/2023       --  Initial coding               :  Kathryn Sarachan
#       4/25/2023       --  Converted to Python 3        :  Susan Krueger
#       5/23/2023       --  SASSIE 3.0 test program      :  Susan Krueger
#       6/28/2023       --  SASSIE 3.0 coding            :  Susan Krueger
#       3/2024  --  modified to use scipy for Guinier fit:  Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **

'''
    **Guinier Analysis** performs a Guinier fit to the data in the range qR:sub:`gmin`\  to qR\:sub:`gmax`\ , where qR:sub:`gmin`\  is defined by the first data point and qR\:sub:`gmax`\  is based on the qR\:sub:`g`\  limit, both as specified by the user.

    Called by **Get Composite Scattering Intensities**
    
    Requires **Polynomial Function** and **scipy.optimize.curve_fit**
    
'''

import math
import scipy.optimize
import numpy
# import sassie.contrast.multi_component_analysis.polynomial_fit as polynomial_fit
# import sassie.contrast.multi_component_analysis.chi_squared_correlation as chi_squared_correlation
import polynomial_fit as polynomial_fit
import chi_squared_correlation as chi_squared_correlation


def guinier_fit(other_self):
    r'''

    **Guinier Fit** fits a first order polynomial to the Guinier equation. The original code was written by Andrew Whitten (2/2006) as part of the MULCh program. Rewritten in Python by Kathryn Sarachan (4/2023). 

    **Reference:** Whitten, A.E., et al. (2007), "MULCh: modules for the analysis of small-angle neutron contrast variation data from biomolecular assemblies", _Journal of Applied Crystallography_ **41**, 222-226. 


    Note
    ----

    The Guinier equation is given by:

    :math:`ln{I(q)} = ln{I_{0}} - \frac{q^2R_{g}^2}{3}`

    Once :math:`I_{0}` and :math:`R_{g}` are found at each contrast, scale factors can be calculated, if the user chooses to do so, using the equation:

    :math:`I_{0} = \frac{cM_{w}}{N_{A}}(\Delta \rho V)^2`

    where 
    :math:`M_{w}` is the total molecular weight of the complex  
    :math:`c` is the total concentration of the complex  
    :math:`N_{A}` is Avogadro's number  
    and 
    :math:`\Delta \rho V = \Delta \rho_{1} V_{1} + \Delta \rho_{2} V_{2}`

    Taking the ratio of the concentrations between two different contrasts, :math:`i` and :math:`j`,

    :math:`\frac{c_{i}}{c_{j}} = \frac{I_{0i}(\Delta \rho V)^2_{j}}{I_{0j}(\Delta \rho V)^2_{i}}`

    Currently, the scale factors are calculated at each contrast with respect to the first data set in the scattering data array. The scale factors at each contrast are returned so that the data can be rescaled before undergoing further processing.  If the measured concentrations are not accurate, rescaling the data in this way can produce a more accurate result for the composite scattering intensities.


    Parameters
    ----------
    number_of_contrast_points: int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    number_of_data_points: int
        The number of points in the contrast variation data files
    scattering data: 3D float array (dimension = number_of_contrast_points x number_of_data_points x 3)
        The scattering data (q, I(q), I\ :sub: `err`\ (q)) at each fraction D\ :sub:`2`\ O
    q_rg_limit_guinier: float array (dimension = number_of_contrast_points)
        qR\ :sub: `g`\  limit for the Guinier analysis at each fraction D\ :sub:`2`\ O
    starting_data_point_guinier: int array (dimension = number_of_contrast_points)
        The index of the starting data point for the Guinier fit at each fraction D\ :sub:`2`\ O  (index of the first data point = 1)
    initial_points_to_use_guinier: int array (dimension = number_of_contrast_points)
        The number of data points to use initially for the Guinier fit at each fraction D\ :sub:`2`\ O  (the final number of points used depends on the qR\ :sub: `g`\  limit)
    initial_guess_guinier: float array (dimension = 2 for a line fit)
        The initial guess for the Guinier fit parameters (default = [1., 1.])
    scale_factor: float array (dimension = number_of_contrast_points)
        The initial scale factor for the data at each fraction D\ :sub:`2`\ O 
    delta_rho_v: float array (dimension = number_of_contrast_points)
        :math:`\Delta \rho V` as defined above at each fraction D\ :sub:`2`\ O 
    refine_scale_factor_flag: boolean
        Indicates whether the scale factor at each fraction D\ :sub:`2`\ O  will be adjusted based on the I(0) values calculated from the Guinier fit and delta_rho_v    

    Returns
    -------
    rg_guinier: float array( dimension = number_of_contrast_points)
        The radius of gyration from the Guinier fit at each fraction D\ :sub:`2`\ O
    rg_guinier_error: float array (dimension = number_of_contrast_points)
        The error in the Gunier radius of gyration at each fraction D\ :sub:`2`\ O
    izero_guinier: float array 
        The I(0) value from the Guinier fit at each fraction D\ :sub:`2`\ O 
    izero_error_guinier: float array
        The error in the Guinier I(0) value at each fraction D\ :sub:`2`\ O 
    points_used_guinier: int array (dimension = number_of_contrast_points)
        The final number of points used in the Guinier fit at each fraction D\ :sub:`2`\ O 
    chi_squared_guinier: float array (dimension = number_of_contrast_points)
        The chi-square value for the Guinier fit at each fraction D\ :sub:`2`\ O 
    q_min_guinier: float array (dimension = number_of_contrast_points)
        The minimum q value for the Guinier fit at each fraction D\ :sub:`2`\ O 
    q_max_guinier: float array (dimension = number_of_contrast_points)
        The maximum q value for the Guinier fit at each fraction D\ :sub:`2`\ O 
    q_rg_min_guinier: float array (dimension = number_of_contrast_points)
        The minimum qR\ :sub: `g`\  value for the Guinier fit at each fraction D\ :sub:`2`\ O 
    q_rg_max_guinier: float array (dimension = number_of_contrast_points)
        The maximum qR\ :sub: `g`\  value for the Guinier fit at each fraction D\ :sub:`2`\ O 
    scale_factor: float array (dimension = number_of_contrast_points)
        The final scale factor at each fraction D\ :sub:`2`\ O  (rescaling is only preformed if refine_scale_factor_flag is True) 

    '''
# mvars used:  number_of_contrast_points, number_of_data_points, scattering_data, q_rg_limit_guinier, starting_data_point_guinier, initial_points_to_use_guinier, refine_scale_factor_flag, initial_guess_guinier
# mcavars used:  number_of_data_points, scale_factor, delta_rho_v
# mcavars returned:  points_used_guinier, chi_squared_guinier, rg_guinier, rg_error_guinier, izero_guinier, izero_error_guinier, q_min_guinier, q_max_guinier, q_rg_min_guinier, q_rg_max_guinier

    log = other_self.log
    log.debug('\nin Guinier fit\n')
    pgui = other_self.run_utils.print_gui

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    log.debug('variables before Guinier fit\n')
    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    mcavars.points_used_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.chi_squared_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.rg_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.rg_error_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.izero_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.izero_error_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.q_rg_min_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.q_rg_max_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.q_min_guinier = numpy.zeros(mvars.number_of_contrast_points)
    mcavars.q_max_guinier = numpy.zeros(mvars.number_of_contrast_points)

    pgui('performing Guinier analysis\n')
#    print('qRg limit: ', mvars.q_rg_limit_guinier)
#    print('number of data points: ', mcavars.number_of_data_points)
    log.debug('scale factor passed to Guinier: ' +
              str(mcavars.scale_factor) + '\n')

    for j in range(mvars.number_of_contrast_points):

        # if the starting point is not the first point, adjust the q**2, ln(I) and ln(I) error arrays accordingly
        start = mvars.starting_data_point_guinier[j] - 1
        number_of_points = mcavars.number_of_data_points - start
        q_squared = numpy.zeros(number_of_points)
        ln_i = numpy.zeros(number_of_points)
        ln_i_error = numpy.zeros(number_of_points)

        for i in range(number_of_points):
            index = start+i
#            print('i, index, q: ', i, index, scattering_data[j][index][0])
            q_squared[i] = mcavars.scattering_data[j][index][0] * \
                mcavars.scattering_data[j][index][0]
#            print('i, q: ', i, math.sqrt(q_squared[i]))
# Katie's note: in C, ln(something negative) throws a floating point error that just goes into the array. Here I needed a place holder to serve the same purpose (for now): 999.999
# NOTE: numpy.log returns a RuntimeWarning and not a ValueError as math.log does and "except RuntimeWarning" doesn't set the new value as desired.
            try:
                ln_i[i] = math.log(mcavars.scattering_data[j][index][1])
            except ValueError as e:
                ln_i[i] = 999.999
            ln_i_error[i] = mcavars.scattering_data[j][index][2] / \
                mcavars.scattering_data[j][index][1]

        status = 0
# points to use is somewhat arbitrary.  It is input by the user since it depends on delta_q as to how many points will be in the Guinier region.
        points_to_use = mvars.initial_points_to_use_guinier[j]
#        print('initial points to use: ', points_to_use)

#        print('initial_guess_guinier: ', mvars.initial_guess_guinier)
        while (status == 0):
         # print('points to use: ', points_to_use)
            # print('q_squared, ln_i: ', q_squared[0:points_to_use], ln_i[0:points_to_use])
            line_fit, line_fit_covariance = scipy.optimize.curve_fit(
                polynomial_fit.polynomial_function, q_squared[0:points_to_use], ln_i[0:points_to_use], sigma=ln_i_error[0:points_to_use], absolute_sigma=True, p0=mvars.initial_guess_guinier)
            log.debug('line fit coefficients: ' + str(line_fit) + '\n')
            log.debug('line fit covariance matrix: ' +
                      str(line_fit_covariance) + '\n')
            ln_i_calculated = polynomial_fit.polynomial_function(
                q_squared[0:points_to_use], *line_fit)
#            print('ln(I) calculated: ', ln_i_calculated)
#            diff = ln_i[0:points_to_use] - ln_i_calculated
#            print('diff: ', diff)

            line_fit_reduced_chi_squared = chi_squared_correlation.get_reduced_chi_squared(
                ln_i[0:points_to_use], ln_i_error[0:points_to_use], ln_i_calculated, points_to_use, 2)

            log.debug('line fit reduced chi-squared: ' +
                      str(line_fit_reduced_chi_squared) + '\n')

# NOTE: To get the correlation matrix and the std deviation (line_fit_error):
            line_fit_error, line_fit_correlation = chi_squared_correlation.correlation_from_covariance(
                line_fit_covariance)
            log.debug('line fit std deviation: ' + str(line_fit_error) + '\n')
#            print('line fit correlation matrix: ', line_fit_correlation)
#            print('line fit correlation coefficient: ', line_fit_correlation[0][1])

            intercept = line_fit[0]
            slope = line_fit[1]
            intercept_error = line_fit_error[0]
            slope_error = line_fit_error[1]
            log.debug('j, slope, error, intercept, error: ' + str(j) + ', ' +
                      str(slope) + ', ' + str(slope_error) + ', ' + str(intercept) + ', ' + str(intercept_error) + '\n')

# Katie's note: another spot where C exception handling is dopey - negative sqrts just populate the answer with a code, rather than throwing an exception. Here, we're going to call any negative sqrt 0.00, which should always be safely < q_rg_limit_guinier
# NOTE: numpy.sqrt returns a RuntimeWarning and not a ValueError as math.sqrt does and "except RuntimeWarning" doesn't set the new value as desired.
#            print('qmax: ', mcavars.scattering_data[j][start + points_to_use - 1][0])
            try:
                limit_test = math.sqrt(-3.0*slope) * \
                    mcavars.scattering_data[j][start + points_to_use - 1][0]
            except ValueError as e:
                limit_test = 0.00
#            print('j, limit_test: ', j, limit_test)

# WHY are we OK with slope > 0?  slope should be negative, otherwise not in a valid Guinier region!
#            if (limit_test < q_rg_limit_guinier) or slope > 0.0:
            if (limit_test < mvars.q_rg_limit_guinier):
                mcavars.points_used_guinier[j] = points_to_use
                mcavars.chi_squared_guinier[j] = line_fit_reduced_chi_squared
#                print('j, points used, chi squared: ', j, mcavars.points_used_guinier[j], mcavars.chi_squared_guinier[j])

                # Rg, where slope = -Rg**2/3; Rg = sqrt(-3.0*slope)
                try:
                    mcavars.rg_guinier[j] = math.sqrt(-3.0*slope)
                except ValueError as e:
                    mcavars.rg_guinier[j] = 0.00
                # RgErr
                try:
                    mcavars.rg_error_guinier[j] = 3.0*math.sqrt(
                        slope_error)/2.0/math.sqrt(-3.0*slope)
                except ValueError as e:
                    mcavars.rg_error_guinier[j] = 0.01
                log.debug('j, Rg, error: ' + str(j) + ', ' +
                          str(mcavars.rg_guinier[j]) + ', ' + str(mcavars.rg_error_guinier[j]) + '\n')

                # I0, where intercept = ln[I(0)]
                mcavars.izero_guinier[j] = math.exp(intercept)
                # I0Err
# TODO: need try: except: below because of square root
                mcavars.izero_error_guinier[j] = math.exp(
                    intercept)*math.sqrt(intercept_error)

                mcavars.q_min_guinier[j] = mcavars.scattering_data[j][start][0]
                mcavars.q_max_guinier[j] = mcavars.scattering_data[j][start +
                                                                      points_to_use - 1][0]
                mcavars.q_rg_min_guinier[j] = mcavars.rg_guinier[j] * \
                    mcavars.scattering_data[j][start][0]
                mcavars.q_rg_max_guinier[j] = mcavars.rg_guinier[j] * \
                    mcavars.scattering_data[j][start + points_to_use - 1][0]
                log.debug('j, qmin, qmax, rg, qrgmin, qrgmax: ' + str(j) + ', ' +
                          str(mcavars.q_min_guinier[j]) + ', ' + str(mcavars.q_max_guinier[j]) + ', ' + str(mcavars.rg_guinier[j]) + ', ' + str(mcavars.q_rg_min_guinier[j]) + ', ' + str(mcavars.q_rg_max_guinier[j]) + '\n')

                points_to_use += 1

            else:
                status = 1

        if (mvars.refine_scale_factor_flag == True):
            mcavars.scale_factor[j] = mcavars.izero_guinier[0]/mcavars.delta_rho_v[0] / \
                mcavars.delta_rho_v[0]/(mcavars.izero_guinier[j] /
                                        mcavars.delta_rho_v[j]/mcavars.delta_rho_v[j])
            log.debug('j, scale factor: ' + str(j) + ', ' +
                      str(mcavars.scale_factor[j]) + '\n')
    log.debug('refined scale factor: ' + str(mcavars.scale_factor) + '\n')
# the refined scale factor gets passed back to the main program

    log.debug('mcavars after Guinier fit\n')
    log.debug(vars(mcavars))

    return
