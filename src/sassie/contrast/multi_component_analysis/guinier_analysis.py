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
#       4/23/2023       --  Initial coding              :  Kathryn Sarachan
#       4/25/2023       --  Converted to Python 3       :  Susan Krueger
#       5/23/2023       --  SASSIE 3.0 test program     :  Susan Krueger
#       6/28/2023       --  SASSIE 3.0 coding           :  Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **

'''
    **Guinier Analysis** performs a Guinier fit to the data in the range qR\ :sub:`gmin`\  to qR\ :sub:`gmax`\ , where qR\ :sub:`gmin`\  is defined by the first data point and qR\ :sub:`gmax`\  is based on the qR\ :sub:`g`\  limit, both as specified by the user.

    Called by **Get Composite Scattering Intensities**
    
    Requires **Polynomial Fit**
    
'''
import io
import os
import math
import sys
import sassie.contrast.multi_component_analysis.polynomial_function_fit as polynomial_function_fit
#import polynomial_function_fit as polynomial_function_fit
import numpy


def guinier_fit(other_self):
    r'''

    **Guinier Fit** fits a first order polynomial to the Guinier equation. The original code was written by Andrew Whitten (2/2006) as part of the MULCh program. Rewritten in Python by Kathryn Sarachan (4/2023). 

    **Reference:** Whitten, A.E., Cai, S, Trewhella, J. (2008). "MULCh: modules for the analysis of small-angle neutron contrast variation data from biomolecular assemblies", *J. Appl. Cryst.* **41**, 222 - 226. 


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
    scattering_data: 3D float array (dimension = number_of_contrast_points x number_of_data_points x 3)
        The scattering data, q, I(q), I\ :sub:`err`\ (q) , at each fraction D\ :sub:`2`\ O 
    q_rg_limit_guinier: float
        qR\ :sub:`g`\  limit for the Guinier analysis; a single value applies for all contrasts
    starting_data_point_guinier: int array (dimension = number_of_contrast_points)
        The index of the starting data point for the Guinier fit at each fraction D\ :sub:`2`\ O  (index of the first data point = 1)
    initial_points_to_use_guinier: int array (dimension = number_of_contrast_points)
        The number of data points to use initially for the Guinier fit at each fraction D\ :sub:`2`\ O  (the final number of points used depends on the qR\ :sub:`g`\  limit) 
    scale_factor: float array (dimension = number_of_contrast_points)
        The initial scale factor for the data at each fraction D\ :sub:`2`\ O 
    delta_rho_v: float array (dimension = number_of_contrast_points)
        :math:`\Delta \rho V` as defined above at each fraction D\ :sub:`2`\ O 
    refine_scale_factor_flag: boolean
        Indicates whether the scale factor at each fraction D\ :sub:`2`\ O  will be adjusted based on the I(0) values calculated from the Guinier fit and delta_rho_v    

    Returns
    -------
    rg_guinier: float array (dimension = number_of_contrast_points)
        The radius of gyration from the Guinier fit at each fraction D\ :sub:`2`\ O
    rg_guinier_error: float array (dimension = number_of_contrast_points)
        The error in the Gunier radius of gyration at each fraction D\ :sub:`2`\ O
    izero_guinier: float array (dimension = number_of_contrast_points)
        The I(0) value from the Guinier fit at each fraction D\ :sub:`2`\ O 
    izero_error_guinier: float array (dimension = number_of_contrast_points)
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
        The minimum qR\ :sub:`g`\  value for the Guinier fit at each fraction D\ :sub:`2`\ O 
    q_rg_max_guinier: float array (dimension = number_of_contrast_points)
        The maximum qR\ :sub:`g`\  value for the Guinier fit at each fraction D\ :sub:`2`\ O 
    scale_factor: float array (dimension = number_of_contrast_points)
        The final scale factor at each fraction D\ :sub:`2`\ O  (rescaling is only preformed if refine_scale_factor_flag is True) 

    '''

    log = other_self.log
    log.debug('in Guinier fit')
    pgui = other_self.run_utils.print_gui

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    fit = numpy.zeros((2, 1))
    correlation = numpy.zeros((2, 2))
    points_used_guinier = numpy.zeros(mvars.number_of_contrast_points)
    chi_squared_guinier = numpy.zeros(mvars.number_of_contrast_points)
    rg_guinier = numpy.zeros(mvars.number_of_contrast_points)
    rg_error_guinier = numpy.zeros(mvars.number_of_contrast_points)
    izero_guinier = numpy.zeros(mvars.number_of_contrast_points)
    izero_error_guinier = numpy.zeros(mvars.number_of_contrast_points)
    q_rg_min_guinier = numpy.zeros(mvars.number_of_contrast_points)
    q_rg_max_guinier = numpy.zeros(mvars.number_of_contrast_points)
    q_min_guinier = numpy.zeros(mvars.number_of_contrast_points)
    q_max_guinier = numpy.zeros(mvars.number_of_contrast_points)

#    print('qRg limit: ', mvars.q_rg_limit_guinier)
#    print('scale factor passed to Guinier: ', mcavars.scale_factor)

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
# Katie's NOTE: In C, ln(something negative) throws a floating point error that just goes into the array. Here I needed a place holder to serve the same purpose (for now): 999.999
# NOTE: numpy.log returns a RuntimeWarning and not a ValueError as math.log does and "except RuntimeWarning" doesn't set the new value as desired. So, we are using math rather than numpy here to capture the ValueError
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

        while (status == 0):
            reduced_chi_squared, fit, correlation = polynomial_function_fit.polynomial_fit(
                1, q_squared, ln_i, ln_i_error, points_to_use)
            reduced_chi_squared = reduced_chi_squared[0, 0]
            correlation00 = correlation[0, 0]
            correlation11 = correlation[1, 1]
            slope = fit[0, 0]
            intercept = fit[1, 0]
#            print('j, slope, intercept: ', j, slope, intercept)

# another spot where C exception handling is dopey - negative sqrts just populate the answer with a code, rather than throwing an exception. Here, we're going to call any negative sqrt 0.00, which should always be safely  < q_rg_limit_guinier
# NOTE: numpy.sqrt returns a RuntimeWarning and not a ValueError as math.sqrt does and "except RuntimeWarning" doesn't set the new value as desired.
#            print('qmax: ', mcavars.scattering_data[j][start + points_to_use - 1][0])
            try:
                limit_test = math.sqrt(-3.0*slope) * \
                    mcavars.scattering_data[j][start + points_to_use - 1][0]
            except ValueError as e:
                limit_test = 0.00
#            print('j, limit_test: ', j, limit_test)

# WHY check for slope > 0? Slope should always be negative, otherwise not in a valid Guinier region?  Testing for qRg limit only.
#            if (limit_test < q_rg_limit_guinier) or slope > 0.0:
            if (limit_test < mvars.q_rg_limit_guinier):
                points_used_guinier[j] = points_to_use
                chi_squared_guinier[j] = reduced_chi_squared
#                print('j, points used, chi squared: ', j, points_used_guinier[j], chi_squared_guinier[j])

                # Rg, where slope = -Rg**2/3; Rg = sqrt(-3.0*slope)
                try:
                    rg_guinier[j] = math.sqrt(-3.0*slope)
                except ValueError as e:
                    rg_guinier[j] = 0.00
                # RgErr
                try:
                    rg_error_guinier[j] = 3.0*math.sqrt(
                        reduced_chi_squared*correlation00)/2.0/math.sqrt(-3.0*slope)
                except ValueError as e:
                    rg_error_guinier[j] = 0.01
#                print('j, Rg, error: ', j, rg_guinier[j], rg_error_guinier[j])

                # I0, where intercept = ln[I(0)]
                izero_guinier[j] = math.exp(intercept)
                # I0Err
                izero_error_guinier[j] = math.exp(
                    intercept)*math.sqrt(reduced_chi_squared*correlation11)

                q_min_guinier[j] = mcavars.scattering_data[j][start][0]
                q_max_guinier[j] = mcavars.scattering_data[j][start +
                                                              points_to_use - 1][0]
                q_rg_min_guinier[j] = rg_guinier[j] * \
                    mcavars.scattering_data[j][start][0]
                q_rg_max_guinier[j] = rg_guinier[j] * \
                    mcavars.scattering_data[j][start + points_to_use - 1][0]
#                print('j, qmin, qmax, rg, qrgmin, qrgmax: ', j, q_min_guinier[j], q_max_guinier[j], rg_guinier[j], q_rg_min_guinier[j], q_rg_max_guinier[j])

                points_to_use += 1

            else:
                status = 1

        if(mvars.refine_scale_factor_flag == True):

            # scale factor = c[0]/C[j] = I0[0]*delta_rho_v[j]**2/(I0[j]/delta_rho_v[0]**2)
            mcavars.scale_factor[j] = izero_guinier[0]/mcavars.delta_rho_v[0] / \
                mcavars.delta_rho_v[0]/(izero_guinier[j] /
                                        mcavars.delta_rho_v[j]/mcavars.delta_rho_v[j])
#            print(j, mcavars.scale_factor[j])
#    print('refined scale factor: ', mcavars.scale_factor)
# the refined scale factor gets passed back to the main program, no need to explicitely return it

    log.debug(vars(mcavars))

    return points_used_guinier, chi_squared_guinier, rg_guinier, rg_error_guinier, izero_guinier, izero_error_guinier, q_min_guinier, q_max_guinier, q_rg_min_guinier, q_rg_max_guinier
