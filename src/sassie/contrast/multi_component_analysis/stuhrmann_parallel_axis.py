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
#    STUHRMANN PARALLEL AXIS is the method that performs both a Sturhmann
#    and parallel axis theorem analysis to obtain the radii of gyration of
#    the individual components in a two component complex and the distance
#    between their center of mass.  The results from the two methods should
#    be similar and can be compared.#
#
#       STUHRMANN is the method that fits Rg**2 vs 1/mvars.delta_rho data to a
#       2nd order polynomial and then calculates the radii of gyration for
#       each of two components and the distance between their center of mass
#
#       PARALLEL AXIS is the method that calculates the radii of gyration for
#       each of two components and the distance between their center of mass
#       using the parallel axis theorem.
#
#       7/18/2021       --  Parallel Axis initial coding:   Kathryn Sarachan
#       2/10/2022       --  Revised for SASSIE 3.0      :   Susan Krueger
#       6/3/2021        --  Stuhrmann initial coding    :   Kathryn Sarachan
#       2/14/2022       --  Revised for SASSIE 3.0      :   Susan Krueger
#       2/2024   --  Rewrote to solve eqs using scipy   :   Susan Krueger
#       3/2024          --  Revised error calculations  :   Susan Krueger
#
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''

    **Stuhrmann Parallel Axis** contains the **Sturhmann** and **Parallel Axis** methods that calculate the radii of gyration of the individual components in a two component complex and the distance between their center of mass. The original code was written by Andrew Whitten (2/2006) as part of the MULCh program. Rewritten in Python by Kathryn Sarachan (6/2021). 
    
    **Reference:** Whitten, A.E., Cai, S, Trewhella, J. (2008). "MULCh: modules for the analysis of small-angle neutron contrast variation data from biomolecular assemblies", *J. Appl. Cryst.* **41**, 222 - 226. 
    
    **Inputs:**
    
    radius of gyration, fraction D\ :sub:`2`\ O for each contrast
    
    contrast values (delta_rho) for each component at each contrast
    
    **Outputs:**
    
    radii of gyration for the components, distance between their centers of mass
    
    Called by the **Multi-component Analysis** module.

    Requires **Polynomial Function**, **Parallel Axis Function**, **Do MC Error Analysis**, **Get Reduced Chi Squared**, **numpy.linalg.solve**, **numpy.random.default_rng().uniform** and **scipy.optimize.curve_fit**.
        
'''

import time
import numpy
import scipy.optimize
import sassie.contrast.multi_component_analysis.polynomial_fit as polynomial_fit
import sassie.contrast.multi_component_analysis.chi_squared_correlation as chi_squared_correlation
# import polynomial_fit as polynomial_fit
# import chi_squared_correlation as chi_squared_correlation
import json


def save_data_to_plot_as_json(other_self, delta_rho_inverse, rg_squared, rg_squared_error, rg_squared_calculated, diff):
    '''
    **Save Data to Plot as JSON** is the method that saves the data to be plotted to a JSON file.
    '''

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    data_dict = {
        '1/contrast': [],
        'Rg^2': [],
        'Rg^2_error': [],
        'calculated_Rg^2': [],
        'Rg^2-calculated_Rg^2': []
    }

    for i in range(mvars.number_of_contrast_points):
        # Append new values to the lists in the dictionary
        data_dict['1/contrast'].append(delta_rho_inverse[i])
        data_dict['Rg^2'].append(rg_squared[i])
        data_dict['Rg^2_error'].append(rg_squared_error[i])
        data_dict['calculated_Rg^2'].append(rg_squared_calculated[i])
        data_dict['Rg^2-calculated_Rg^2'].append(diff[i])

    json_data = json.dumps(data_dict)

    mcavars.json_outfile.write(json_data)
    mcavars.json_outfile.close()

    return


def parallel_axis_function(x, r1_squared, r2_squared, cm_distance_squared):
    r'''The **Parallel Axis Function**. Called from **Parallel Axis**. 

    Note
    ----

    The radii of gyration for the two components, :math:`R_{1}` and :math:`R_{2}`, and the distance between their centers of mass, :math:`D`, are obtained by simultaneously solving the equation for the parallel axis theorem:

    :math:`R_{g}^2 = \frac{{\Delta \rho}_{1} V_{1}}{{\Delta \rho}_{1} V_{1} + {\Delta \rho}_{2} V_{2}} R_{1}^2 + \frac{{\Delta \rho}_{2} V_{2}}{{\Delta \rho}_{1} V_{1} + {\Delta \rho}_{2} V_{2}} R_{2}^2 + \frac{{\Delta \rho}_{1} {\Delta \rho}_{2} V_{1} V_{2}}{({\Delta \rho}_{1} V_{1} + {\Delta \rho}_{2} V_{2})^2} D^2`

    Parameters
    ----------
    x:  float array (dimension = 3 x number_of_contrast_points)
        :math:`x` coefficients at each contrast 

    Returns
    -------
    r1_squared:  float
        The square of the radius of gyration for component 1, :math:`R_{1}^2`
    r2_squared:  float
        The square of the radius of gyration for component 2, :math:`R_{2}^2`
    cm_distance_squared:  float
        The square of the distance between the center of mass of the two components, :math:`D^2`

    '''
    f = x[0] * r1_squared + x[1] * r2_squared + x[2] * cm_distance_squared
    return f


def parallel_axis(other_self):
    r'''

    The **Parallel Axis** method uses the parallel axis therorem to solve for the radii of gyration for the two components, :math:`R_{1}` and :math:`R_{2}`, and the distance between their centers of mass, :math:`D`.

    Executed if stuhrmann_parallel_axis_flag == True.

    Notes:

        Output file with chosen output_file_name is stored in the output file directory.
        TODO: Need boilerplate language to add variable name/path for the directory. 

        Output file contains:
            - module variable input parameters
            - radii of gyration for each of two components and the distance between their centers of mass 
                - :math:`R_{1}`, error in :math:`R_{1}`, :math:`R_{2}`, error in :math:`R_{2}`, :math:`D`, the distance between the centers of mass of the components, error in :math:`D`

    This method calls **Parallel Axis Function**, **Get Reduced Chi Squared** and **scipy.optimize.curve_fit**

    Parameters
    ----------
    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    volume_fraction: float array (dimension = number_of_components)
        The volume fraction of each component        
    radius_of_gyration: float array( dimension = number_of_contrast_points)
        The radius of gyration at each fraction D\ :sub:`2`\ O
    radius_of_gyration_error: float array (dimension = number_of_contrast_points)
        The error in radius of gyration at each fraction D\ :sub:`2`\ O
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
    outfile: string
        output file name (with full path): path + output_file_name

    Returns
    -------
    r1_parallel:  float
        The radius of gyration for component 1, :math:`R_{1}`, obtained from the parallel axis theorem analysis
    r1_parallel_error: float
        The error on :math:`R_{1}`
    r2_parallel:  float
        The radius of gyration for component 2, :math:`R_{2}`, obtained from the parallel axis theorem analysis
    r2_parallel_error:  float
        The error on :math:`R_{2}`
    cm_distance_parallel:  float
        The distance between the center of mass of the two components, :math:`D`, obtained from the parallel axis theorem analysis
    cm_distance_parallel_error:  float
        The error on :math:`D`
    chi_squared_parallel:  float
        The reduced chi-squared for fit to the parallel axis equation, i.e, chi-squared divided by: the number of equations (contrast points) - the number of unknowns (3); if there are only 3 contrast points, the reduced chi-squared is not reported and a message is written to the output file

    '''

# mvars used:  number_of_contrast_points, delta_rho, radius_of_gyration, radius_of_gyration_error
# mcavars used:  outfile, volume_fraction

#    print('in parallel axis')

    log = other_self.log
    log.debug('in parallel axis')
    pgui = other_self.run_utils.print_gui

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    log.debug('volume fraction in parallel axis: ' +
              str(mcavars.volume_fraction) + '\n')

    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in range(60)])

    pgui('\n%s \n' % (st))
    pgui('DATA FROM RUN: %s \n\n' % (ttxt))

    pgui('results written to output file: %s \n' %
         (mcavars.multi_component_analysis_path+mvars.output_file_name))
    pgui('-------------------------------\n')
    pgui('\nNumber of points fit: ' + str(mvars.number_of_contrast_points) + '\n')

# Parallel axis analysis

    pgui('Parallel Axis Method\n')
    pgui('Calculating individual Rg values and distance between centers of mass \n')

    rg_squared = numpy.zeros(mvars.number_of_contrast_points)
    rg_squared_error = numpy.zeros(mvars.number_of_contrast_points)
    x = numpy.zeros((mvars.number_of_contrast_points, 3))
    for i in range(mvars.number_of_contrast_points):
        rg_squared[i] = mvars.radius_of_gyration[i]*mvars.radius_of_gyration[i]
        rg_squared_error[i] = 2.0*mvars.radius_of_gyration[i] * \
            mvars.radius_of_gyration_error[i]

    log.debug('rg_squared, rg_squared_error: ' + str(rg_squared) +
              ' +/- ' + str(rg_squared_error) + '\n')

    # define the coefficients for parallel axis equation
    # Rg**2 = f1*Rg1**2 + f2*Rg2**2 + f1*f2*D**2, where f1 = drho1*vf1/(drho1*vf1 + drho2*vf2), f2 = 1 - f1
    # vf1 = V1/(V1+V2); vf2 = V2/(V1+V2) = 1 - vf1
    # Note: The equation used here is equivalent to that given in the docstring above since the total volume (V1+V2) cancels out.

    for k in range(mvars.number_of_contrast_points):
        # print('k, mvars.delta_rho[k][0], mvars.delta_rho[k][1]: ', k, mvars.delta_rho[k][0], mvars.delta_rho[k][1])
        # f1 = drho1*vf1/(vf1*drho1 + vf2*drho2)
        x0 = mvars.delta_rho[k][0]*mcavars.volume_fraction[0] / \
            (mvars.delta_rho[k][0]*mcavars.volume_fraction[0] +
             mvars.delta_rho[k][1]*mcavars.volume_fraction[1])
        x1 = 1.0 - x0  # f2 = 1 - f1
        x2 = x0*x1  # f1*f2
        x[k] = (x0, x1, x2)

#    print ('x: ', x)
#    print ('length of x: ', len(x))
#    print ('rg_squared: ', rg_squared)
#    print ('length of rg_squared: ', len(rg_squared))

    # transpose x to have correct inputs for curve_fit, which accepts a 1 column matrix with number of rows = number of contrast points.  We have a 1 row matrix of length = number of contrast points.
    x = x.T

#    print('transposed x: ', x)
#    print ('length of transposed x: ', len(x))

    fit, covariance = scipy.optimize.curve_fit(
        parallel_axis_function, x, rg_squared, sigma=rg_squared_error, absolute_sigma=True)

    log.debug('fit: ' + str(fit) + '\n')
    log.debug('covariance: ' + str(covariance) + '\n')

    rg_squared_calculated = parallel_axis_function(x, *fit)
    log.debug('rg_squared_calculated: ' + str(rg_squared_calculated) + '\n')

    reduced_chi_squared = chi_squared_correlation.get_reduced_chi_squared(
        rg_squared, rg_squared_error, rg_squared_calculated, mvars.number_of_contrast_points, 3)

    chi_squared_parallel = reduced_chi_squared
    rg1_squared = fit[0]
    rg1_squared_error = numpy.sqrt(covariance[0][0])
    rg2_squared = fit[1]
    rg2_squared_error = numpy.sqrt(covariance[1][1])
    cm_distance_squared = fit[2]
    cm_distance_squared_error = numpy.sqrt(covariance[2][2])
#    print(rg1_squared, rg1_squared_error)
#    print(rg2_squared, rg2_squared_error)
#    print(cm_distance_squared, cm_distance_squared_error)
    rg1_parallel = numpy.sqrt(rg1_squared)
    rg1_error_parallel = 0.5/rg1_parallel*rg1_squared_error
    rg2_parallel = numpy.sqrt(rg2_squared)
    rg2_error_parallel = 0.5/rg2_parallel*rg2_squared_error
    cm_distance_parallel = numpy.sqrt(cm_distance_squared)
    cm_distance_error_parallel = 0.5/cm_distance_parallel*cm_distance_squared_error
#    print(chi_squared_parallel)
#    print(rg1_parallel, rg1_error_parallel)
#    print(rg2_parallel, rg2_error_parallel)
#    print(cm_distance_parallel, cm_distance_error_parallel)

    pgui('Results from Parallel Axis Theorem Analysis: \n')
    pgui('Errors are one standard deviation\n')
    pgui('R1: ' + str(round(rg1_parallel, 4)) + ' +/- ' +
         str(round(rg1_error_parallel, 4)) + '\n')
    pgui('R2: ' + str(round(rg2_parallel, 4)) + ' +/- ' +
         str(round(rg2_error_parallel, 4)) + '\n')
    pgui('D: ' + str(round(cm_distance_parallel, 4)) + ' +/- ' +
         str(round(cm_distance_error_parallel, 4)) + '\n')
    if (chi_squared_parallel >= 0):
        pgui('reduced chi-squared parallel axis: ' +
             str(round(chi_squared_parallel, 4)) + '\n\n')
    else:
        pgui('reduced chi-squared parallel axis: N/A (number of contrast points = number of unknowns for the parallel axis equation)\n\n')

    mcavars.outfile.write('Results from Parallel Axis Theorem Analysis: \n')
    mcavars.outfile.write('Errors are one standard deviation\n')
    mcavars.outfile.write('R1: ' + str(round(rg1_parallel, 4)) +
                          ' +/- ' + str(round(rg1_error_parallel, 4)) + '\n')
    mcavars.outfile.write('R2: ' + str(round(rg2_parallel, 4)) +
                          ' +/- ' + str(round(rg2_error_parallel, 4)) + '\n')
    mcavars.outfile.write('D: ' + str(round(cm_distance_parallel, 4)) +
                          ' +/- ' + str(round(cm_distance_error_parallel, 4)) + '\n')
    if (chi_squared_parallel >= 0):
        mcavars.outfile.write('reduced chi-squared parallel axis: ' +
                              str(round(chi_squared_parallel, 4)) + '\n\n')
    else:
        mcavars.outfile.write(
            'reduced chi-squared parallel axis: N/A (number of contrast points = number of unknowns for the parallel axis equation)\n\n')

    pgui('\n%s \n' % (st))

    time.sleep(0.5)

    return


def stuhrmann(other_self):
    r'''

    The **Stuhrmann** method fits :math:`R_{g}^2` vs :math:`\frac{1}{\Delta \rho}` data to a 2nd order polynomial and then uses the results to calculates the radii of gyration for each of two components and the distance between their centers of mass.  

    Executed if stuhrmann_parallel_axis_flag == True.

    Notes:

        Output file with chosen output_file_name is stored in the output file directory.
        TODO: Need boilerplate language to add variable name/path for the directory. 

        Output file contains:
            - module variable input parameters
            - data for Stuhrmann plot
                - :math:`\frac{1}{\Delta \rho}`, :math:`R_{g}^2`, error in :math:`R_{g}^2`, calculated :math:`R_{g}^2`, :math:`R_{g}^2` - calculated :math:`R_{g}^2`
            - the coefficients of the Stuhrmann equation
                - :math:`\alpha`, error in :math:`\alpha`, :math:`\beta`, error in :math:`\beta`, :math:`R_{m}`, the :math:`R_{g}` of an equivalent homogeneous particle, error in :math:`R_{m}`
            - radii of gyration for each of two components and the distance between their centers of mass 
                - :math:`R_{1}`, error in :math:`R_{1}`, :math:`R_{2}`, error in :math:`R_{2}`, :math:`D`, the distance between the centers of mass of the components, error in :math:`D`

    Requires **Polynomial Function** **Do MC Error Analysis**, **Get Reduced Chi Squared**, **numpy.linalg.solve** and **scipy.optimize.curve_fit**.


    Note
    ----

    :math:`R_{g}` is the measured radius of gyration at each contrast (fraction D\ :sub:`2`\ O in the solvent)

    :math:`\Delta \rho` is the calculated contrast of the complex under the same solvent conditions.

    :math:`\Delta \rho` is calculated from the contrasts of the individual components, :math:`\Delta \rho_{1}` and :math:`\Delta \rho_{2}`.  

    :math:`\Delta \rho = f_{1} \Delta \rho_{1} + f_{2} \Delta \rho_{2}`

    where the volume fractions of the components, :math:`f_{1}` and :math:`f_{2}`, are defined as

    :math:`f_{1} = \frac{V_{1}}{(V_{1} + V_{2})}`

    :math:`f_{2} = 1 - f_{1}`

    The volumes of the components are calculated from their molecular weights and partial specific volumes:

    :math:`V = \overline{v} \frac {M_{w}}{N_{A}}` 

    The Stuhrmann equation:

    :math:`R_{g}^2 = R_{m}^2 + \frac{\alpha}{\Delta \rho} - \frac{\beta}{(\Delta \rho)^2}`

    is used to obtain :math:`\alpha`, :math:`\beta` and :math:`R_{m}^2`, where :math:`R_{m}` is the radius of gyration of the complex at infinite contrast , i.e., the radius of gyration of an equivalent homogeneous particle.

    Then, the radii of gyration for the two components, :math:`R_{1}` and :math:`R_{2}`, and the distance between their centers of mass, :math:`D`, are obtained by simultaneously solving:

    :math:`R_{m}^2 = f_{1} R_{1}^2 + f_{2} R_{2}^2 + f_{1} f_{2} D^2`

    :math:`\alpha = ({\Delta \rho}_{1} - {\Delta \rho}_{2}) f_{1} f_{2} [R_{1}^2 - R_{2}^2 + (f_{2}^2 - f_{1}^2) D^2]`

    :math:`\beta = ({\Delta \rho}_{1} - {\Delta \rho}_{2})^2 f_{1}^2 f_{2}^2 D^2`

    **Reference:** Olah, G.A., et al. (1994). "Troponin I Encompasses an Extended Troponin C in the Ca\ :sup:`2+`\ -Bound Complex: A Small-Angle X-ray and Neutron Scattering Study", *Biochemistry* **33**, 8233 - 8239.


    Parameters
    ----------
    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    volume_fraction: float array (dimension = number_of_components)
        The volume fraction of each component
    radius_of_gyration: float array( dimension = number_of_contrast_points)
        The radius of gyration at each fraction D\ :sub:`2`\ O
    radius_of_gyration_error: float array (dimension = number_of_contrast_points)
        The error in radius of gyration at each fraction D\ :sub:`2`\ O
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
    outfile: string
        output file name (with full path): path + output_file_name
    initial_guess_stuhrmann: float array (dimension = order of polynomial + 1 = 3)
        The initial guess for the coefficients of the 2nd order polynomial (default = [1.0, 1.0, 1.0]; can be changed as an advanced option)

    Returns
    -------
    alpha:  float
        The Stuhrmann equation coefficient, :math:`\alpha`
    alpha_error:  float
        The error on :math:`\alpha`
    beta:  float
        The Stuhrmann equation coefficient, :math:`\beta`
    beta_error:  float
        The error on :math:`\beta`
    rg_infinite_contrast:  float
        The Stuhrmann equation coefficient, :math:`R_{m}`
    rg_infinite_contrast_error:  float
        The error on :math:`R_{m}`
    r1_stuhrmann:  float
        The radius of gyration for component 1, :math:`R_{1}`, obtained from the Stuhrmann analysis
    r1_stuhrmann_error: float
        The error on :math:`R_{1}`
    r2_stuhrmann:  float
        The radius of gyration for component 2, :math:`R_{2}`, obtained from the Stuhrmann analysis
    r2_stuhrmann_error:  float
        The error on :math:`R_{2}`
    cm_distance_stuhrmann:  float
        The distance between the center of mass of the two components, :math:`D`, obtained from the Stuhrmann analysis
    cm_distance_stuhrmann_error:  float
        The error on :math:`D`
    chi_squared_stuhrmann: float
        The reduced chi-squared for the polynomial fit to the Stuhrmann equation, i.e, chi-squared divided by: the number of equations (contrast points) - the number of unknowns (3); if there are only 3 contrast points, the reduced chi-squared is not reported and a message is written to the output file

    '''
# mvars used:  number_of_contrast_points, delta_rho, radius_of_gyration, radius_of_gyration_error, initial_guess_stuhrmann
# mcavars used:  outfile, volume_fraction

    log = other_self.log
    log.debug('in sturhmann')
    pgui = other_self.run_utils.print_gui

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    pgui('Stuhrmann Analysis Method\n')
    pgui('Calculating alpha, beta and Rm')

    delta_rho_inverse = numpy.zeros(mvars.number_of_contrast_points)
    rg_squared = numpy.zeros(mvars.number_of_contrast_points)
    rg_squared_error = numpy.zeros(mvars.number_of_contrast_points)
    rg_squared_calculated = numpy.zeros(mvars.number_of_contrast_points)
    diff = numpy.zeros(mvars.number_of_contrast_points)

    log.debug('volume fraction in stuhrmann: ' +
              str(mcavars.volume_fraction) + '\n')

    # Rg**2 = Rm**2 + alpha/delta_rho - beta/delta_rho**2
    # delta_rho_inverse are the x-coordinates (1/delta_rho) obtained from calculated delta_rho values
    # rg_squared are the y-coordinates (Rg**2) obtained from the experimental Rg values, and
    # rg_squared_error are the errors in the y-coordinates propagated from the experimental error in Rg: err(Rg**2) = 2*Rg*err(Rg)

    for i in range(mvars.number_of_contrast_points):
        delta_rho_inverse[i] = 1.0/(mvars.delta_rho[i][0] *
                                    mcavars.volume_fraction[0] + mvars.delta_rho[i][1]*mcavars.volume_fraction[1])
        rg_squared[i] = mvars.radius_of_gyration[i]*mvars.radius_of_gyration[i]
        rg_squared_error[i] = 2.0*mvars.radius_of_gyration[i] * \
            mvars.radius_of_gyration_error[i]
    log.debug('delta_rho_inverse, rg_squared, rg_squared_error: ' + str(delta_rho_inverse) +
              ', ' + str(rg_squared) + ' +/- ' + str(rg_squared_error) + '\n')

#    print('initial guess for polynomial coefficients: ', mvars.initial_guess_stuhrmann)
    fit, covariance = scipy.optimize.curve_fit(
        polynomial_fit.polynomial_function, delta_rho_inverse, rg_squared, sigma=rg_squared_error, absolute_sigma=True, p0=mvars.initial_guess_stuhrmann)

    log.debug('fit: ' + str(fit) + '\n')
    log.debug('covariance: ' + str(covariance) + '\n')

    rg_squared_calculated = polynomial_fit.polynomial_function(
        delta_rho_inverse, *fit)
    log.debug('rg_squared_calculated: ' + str(rg_squared_calculated) + '\n')

    diff = rg_squared - rg_squared_calculated

    reduced_chi_squared = chi_squared_correlation.get_reduced_chi_squared(
        rg_squared, rg_squared_error, rg_squared_calculated, mvars.number_of_contrast_points, 3)

    chi_squared_stuhrmann = reduced_chi_squared
    beta = - fit[2]
    beta_error = numpy.sqrt(covariance[2][2])
    alpha = fit[1]
    alpha_error = numpy.sqrt(covariance[1][1])
    rg_infinite_contrast_squared = fit[0]
    rg_infinite_contrast_squared_error = numpy.sqrt(covariance[0][0])
    rg_infinite_contrast = numpy.sqrt(rg_infinite_contrast_squared)
    rg_infinite_contrast_error = 0.5/rg_infinite_contrast * \
        rg_infinite_contrast_squared_error

    log.debug('\nbeta, beta_error: ' + str(beta) +
              ' +/- ' + str(beta_error) + '\n')
    log.debug('alpha, alpha_error: ' + str(alpha) +
              ' +/- ' + str(alpha_error) + '\n')
    log.debug('rg_infinite_contrast_squared, rg_infinite_contrast_squared_error: ' +
              str(rg_infinite_contrast_squared) + ' +/- ' + str(rg_infinite_contrast_squared_error) + '\n')
    log.debug('rg_infinite_contrast, rg_infinite_contrast_error: ' +
              str(rg_infinite_contrast) + ' +/- ' + str(rg_infinite_contrast_error) + '\n')

# screen output
    pgui('Results from Stuhrmann Analysis: \n')
    pgui('Errors are one standard deviation\n')
    pgui('beta: ' + str(round(beta, 4)) +
         ' +/- ' + str(round(beta_error, 4)) + '\n')
    pgui('alpha: ' + str(round(alpha, 4)) +
         ' +/- ' + str(round(alpha_error, 4)) + '\n')
    if (alpha > 0):
        pgui(
            'alpha is positive, indicating that the component with a higher scattering length density lies toward the periphery of the complex.\n')
    elif (alpha < 0):
        pgui(
            'alpha is negative, indicating that the component with a higher scattering length density lies toward the interior of the complex.\n')
    pgui('Rm: ' + str(round(rg_infinite_contrast, 4)) +
         ' +/- ' + str(round(rg_infinite_contrast_error, 4)) + '\n')
    if (chi_squared_stuhrmann >= 0):
        pgui('reduced chi-squared Stuhrmann: ' +
             str(round(chi_squared_stuhrmann, 4)) + '\n\n')
    else:
        pgui('reduced chi-squared Stuhrmann: N/A (number of contrast points = number of unknowns for the Stuhrmann equation)\n')
# do we want to output the following to the screen?
    pgui('1/delta_rho\t Rg^2 exp\t Rg^2err\t Rg^2 calc\t diff\n')
    for i in range(mvars.number_of_contrast_points):
        pgui('%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\n' % (
            delta_rho_inverse[i], rg_squared[i], rg_squared_error[i], rg_squared_calculated[i], diff[i]))

# Now solve for R1, R2 and D using alpha, beta and Rm**2
# B is the RHS of the equations 5a-c given by Olah 1994; we are solving 3 equations in 3 unknowns using numpy.linalg.solve.  Since this is an exact solution (3 equations and 3 unknowns), the errors on alpha, beta and Rm will be propagated using a MC error analysis.
# NOTE: In the original code, only the delta_rho values from the first contrast point are used, i.e. delta_rho[0]. The results should be the same no matter which contrast values are used. This was tested and found to be the case, so we will use the mean R1, R2 and D values obtained using all of the contrasts.  These mean values for all MC trials will be averaged to obtain the final R1, R2, D values and the standard errors of the mean will be calculated for each and reported as the errors on R1, R2 and D.

    pgui('Calculating individual Rg values and distance between centers of mass from alpha, beta and Rm**2 \n')

# initial calculation for R1, R2, D using delta_rho at all contrasts
    B = numpy.zeros((3, 3))
    solution = []

    for i in range(mvars.number_of_contrast_points):
        # print('i, mvars.delta_rho[i][0], mvars.delta_rho[i][1]: ', i, mvars.delta_rho[i][0], mvars.delta_rho[i][1])
        B[0][0] = mcavars.volume_fraction[0]
        B[0][1] = mcavars.volume_fraction[1]
        B[0][2] = mcavars.volume_fraction[0]*mcavars.volume_fraction[1]
        B[1][0] = (mvars.delta_rho[i][0] - mvars.delta_rho[i][1]) * \
            mcavars.volume_fraction[0]*mcavars.volume_fraction[1]
        B[1][1] = -(mvars.delta_rho[i][0] - mvars.delta_rho[i][1]) * \
            mcavars.volume_fraction[0]*mcavars.volume_fraction[1]
        B[1][2] = (mvars.delta_rho[i][0] - mvars.delta_rho[i][1])*mcavars.volume_fraction[0]*mcavars.volume_fraction[1] * \
            (mcavars.volume_fraction[1]*mcavars.volume_fraction[1] -
             mcavars.volume_fraction[0]*mcavars.volume_fraction[0])
        B[2][0] = 0.0
        B[2][1] = 0.0
# Note the minus sign in the equation below.  This is because fit[2] (beta) is negative.  However, beta is positive according the eq 5c in Olah 1994.
        B[2][2] = -(mvars.delta_rho[i][0] - mvars.delta_rho[i][1])*(mvars.delta_rho[i][0] - mvars.delta_rho[i][1]) * \
            mcavars.volume_fraction[0]*mcavars.volume_fraction[0] * \
            mcavars.volume_fraction[1]*mcavars.volume_fraction[1]

        this_solution = numpy.linalg.solve(B, fit)
        solution.append(this_solution)
    log.debug('solution: ' + str(solution) + '\n')

    r12r22d2 = numpy.array(solution)
    mean_r12r22d2 = numpy.mean(r12r22d2, axis=0)
    r1r2d = numpy.sqrt(r12r22d2)
    mean_r1r2d = numpy.mean(r1r2d, axis=0)
#    print('mean r1r2d: ', mean_r1r2d)
#    print('r12r22d: ', r12r22d2)
#    print('mean r12r22d2: ', mean_r12r22d2)
#    print('r1r2d: ', r1r2d)
#    print('mean r1r2d: ', mean_r1r2d)

# these are not used further; kept for debugging
#    initial_rg1_squared = mean_r12r22d2[0]
#    initial_rg2_squared = mean_r12r22d2[1]
#    initial_cm_distance_squared = mean_r12r22d2[2]
#    initial_rg1_stuhrmann = mean_r1r2d[0]
#    initial_rg2_stuhrmann = mean_r1r2d[1]
#    initial_cm_distance_stuhrmann = mean_r1r2d[2]
    log.debug('initial mean R1**2, R2**2, D**2: '+str(mean_r12r22d2)+'\n')
    log.debug('initial mean R1, R2, D: '+str(mean_r1r2d)+'\n')

# MC error analysis
    pgui('Performing Monte Carlo error analysis')
    ntrials = 100
    r12r22d2_list, r1r2d_list, number_of_points = do_mc_error_analysis(
        mean_r12r22d2, mean_r1r2d, mvars.number_of_contrast_points, rg_infinite_contrast_squared, rg_infinite_contrast_squared_error, alpha, alpha_error, beta, beta_error, mcavars.volume_fraction, mvars.delta_rho, ntrials)

    number_of_points = len(r1r2d_list)
    log.debug('number of points for MC error analysis: ' +
              str(number_of_points) + '\n')

    final_mean_r12r22d2 = numpy.mean(r12r22d2_list, axis=0)
#    print('final_mean_r12r22d2: ', final_mean_r12r22d2)
    final_mean_r1r2d = numpy.mean(r1r2d_list, axis=0)
#    print('final_mean_r1r2d: ', final_mean_r1r2d)
    stddev_r1r2d = numpy.std(r1r2d_list, axis=0)
#    print('stddev: ', stddev_r1r2d)
    stderr_r1r2d = stddev_r1r2d/numpy.sqrt(number_of_points)
#    print('stderr: ', stderr_r1r2d)
    log.debug('final mean R1**2, R2**2, D**2: '+str(final_mean_r12r22d2)+'\n')
    log.debug('final mean R1, R2, D: '+str(final_mean_r1r2d)+'\n')
    log.debug('standard deviation of R1, R2, D: '+str(stddev_r1r2d)+'\n')
    log.debug('standard error of the mean of R1, R2, D: ' +
              str(stderr_r1r2d)+'\n')

# these are not used further; kept for debugging
#    rg1_squared = final_mean_r12r22d2[0]
#    rg2_squared = final_mean_r12r22d2[1]
#    cm_distance_squared = final_mean_r12r22d2[2]

    rg1_stuhrmann = final_mean_r1r2d[0]
    rg2_stuhrmann = final_mean_r1r2d[1]
    cm_distance_stuhrmann = final_mean_r1r2d[2]
    rg1_error_stuhrmann = stderr_r1r2d[0]
    rg2_error_stuhrmann = stderr_r1r2d[1]
    cm_distance_error_stuhrmann = stderr_r1r2d[2]
# If we want a 95% confidence interval, it is defined by 1.96*standard_error

# screen output

    pgui('Errors on R1, R2 and D are standard errors of the mean based on Monte Carlo calculations using alpha +/- error, beta +/- error and Rm**2 +/- error\n')
    pgui('R1: ' + str(round(rg1_stuhrmann, 4)) +
         ' +/- ' + str(round(rg1_error_stuhrmann, 4)) + '\n')
    pgui('R2: ' + str(round(rg2_stuhrmann, 4)) +
         ' +/- ' + str(round(rg2_error_stuhrmann, 4)) + '\n')
    pgui('D: ' + str(round(cm_distance_stuhrmann, 4)) +
         ' +/- ' + str(round(cm_distance_error_stuhrmann, 4)) + '\n')

# output file
    mcavars.outfile.write('Results from Stuhrmann Analysis: \n')
    mcavars.outfile.write('Errors are one standard deviation\n')
    mcavars.outfile.write('beta: ' + str(round(beta, 4)) +
                          ' +/- ' + str(round(beta_error, 4)) + '\n')
    mcavars.outfile.write('alpha: ' + str(round(alpha, 4)) +
                          ' +/- ' + str(round(alpha_error, 4)) + '\n')
    if (alpha > 0):
        mcavars.outfile.write(
            'alpha is positive, indicating that the component with a higher scattering length density lies toward the periphery of the complex.\n')
    elif (alpha < 0):
        mcavars.outfile.write(
            'alpha is negative, indicating that the component with a higher scattering length density lies toward the interior of the complex.\n')
    mcavars.outfile.write('Rm: ' + str(round(rg_infinite_contrast, 4)) +
                          ' +/- ' + str(round(rg_infinite_contrast_error, 4)) + '\n')
    if (chi_squared_stuhrmann >= 0):
        mcavars.outfile.write('reduced chi-squared Stuhrmann: ' +
                              str(round(chi_squared_stuhrmann, 4)) + '\n\n')
    else:
        mcavars.outfile.write(
            'reduced chi-squared Stuhrmann: N/A (number of contrast points = number of unknowns for the Stuhrmann equation)\n\n')
    mcavars.outfile.write(
        'Errors on R1, R2 and D are standard errors of the mean based on Monte Carlo calculations using alpha +/- error, beta +/- error and Rm**2 +/- error\n')
    mcavars.outfile.write('R1: ' + str(round(rg1_stuhrmann, 4)) +
                          ' +/- ' + str(round(rg1_error_stuhrmann, 4)) + '\n')
    mcavars.outfile.write('R2: ' + str(round(rg2_stuhrmann, 4)) +
                          ' +/- ' + str(round(rg2_error_stuhrmann, 4)) + '\n')
    mcavars.outfile.write('D: ' + str(round(cm_distance_stuhrmann, 4)) +
                          ' +/- ' + str(round(cm_distance_error_stuhrmann, 4)) + '\n')

    mcavars.outfile.write(
        '\n1/delta_rho\t Rg^2 exp\t Rg^2err\t Rg^2 calc\t diff\n')
    for i in range(mvars.number_of_contrast_points):
        mcavars.outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\n' % (
            delta_rho_inverse[i], rg_squared[i], rg_squared_error[i], rg_squared_calculated[i], diff[i]))
    mcavars.outfile.close()

    save_data_to_plot_as_json(other_self, delta_rho_inverse,
                              rg_squared, rg_squared_error, rg_squared_calculated, diff)

    time.sleep(1.0)

    return

# TODO: since the MC error analysis method isn't used by other methods, it can be rewritten to use mvars and mcavars


def do_mc_error_analysis(mean_r12r22d2, mean_r1r2d, number_of_contrast_points, rg_infinite_contrast_squared, rg_infinite_contrast_squared_error, alpha, alpha_error, beta, beta_error, volume_fraction, delta_rho, ntrials):
    r'''
    Method to perform Monte Carlo error analysis to take into account the errors on :math:`\alpha`, :math:`\beta` and :math:`R_{m}` that were determined from the polynomial fit to the Stuhrmann equation.  The three equations that relate :math:`\alpha`, :math:`\beta` and :math:`R_{m}` to :math:`R_{g}` for the individual components (:math:`R_{1}` and :math:`R_{2}`) and the distance between the centers of mass of the components (:math:`D`), are solved many times with different values of :math:`\alpha`, :math:`\beta` and :math:`R_{m}` constrained by +/- their errors.  A **uniform** random number generator is used so that any value within +/- error is equally likely to be generated.  The final :math:`R_{1}`, :math:`R_{2}` and :math:`D` values are the mean values from all of the Monte Carlo trials, with errors calculated as the standard errors of the mean.

    This method calls **numpy.random.default_rng().uniform** and **numpy.linalg.solve**

    Parameters
    ----------
    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    mean_r12r22d2: float array (dimension = 3)
        The initial mean value of :math:`R_{g}^2` for the individual components (:math:`R_{1}^2` and :math:`R_{2}^2`) and the distance between the centers of mass of the components (:math:`D^2`) determined from calculating the values at each contrast
    mean_r1r2d: float array (dimension = 3)
        The intial value of :math:`R_{g}` for the individual components (:math:`R_{1}` and :math:`R_{2}`) and the distance between the centers of mass of the components (:math:`D`) determined from :math:`R_{1}^2`, :math:`R_{2}^2` and :math:`D^2`
    alpha:  float
        The Stuhrmann equation coefficient, :math:`\alpha`
    alpha_error:  float
        The error on :math:`\alpha`
    beta:  float
        The Stuhrmann equation coefficient, :math:`\beta`
    beta_error:  float
        The error on :math:`\beta`
    rg_infinite_contrast:  float
        The Stuhrmann equation coefficient, :math:`R_{m}^2`
    rg_infinite_contrast_error:  float
        The error on :math:`R_{m}^2`
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
    volume_fraction: float array (dimension = number_of_components)
        The volume fraction of each component
    ntrials: int
        The number of Monte Carlo trials

    Returns
    -------
    r12r22d2_list: 2D float array (dimension = ntrials+1 x 3)
        The mean values of :math:`R_{g}^2` for the individual components (:math:`R_{1}^2` and :math:`R_{2}^2`) and the distance between the centers of mass of the components (:math:`D^2`) determined from calculating the values at each contrast; includes the initial values and those from each Monte Carlo trial
    r1r2d_list: 2D float array (dimension = ntrials+1 x 3)
        :math:`R_{g}` for the individual components (:math:`R_{1}` and :math:`R_{2}`) and the distance between the centers of mass of the components (:math:`D`); includes the initial values and those from each Monte Carlo trial        
    number_of_points: float
        the length of the r1r2d_list array; used in determining the final averaged :math:`R_{1}^2`, :math:`R_{2}^2` and :math:`D^2` values

    '''

    r12r22d2_list = numpy.zeros((ntrials+1, 3))
    r12r22d2_list[0] = mean_r12r22d2
    r1r2d_list = numpy.zeros((ntrials+1, 3))
    r1r2d_list[0] = mean_r1r2d
#    print(r12r22d2_list)
#    print(r1r2d_list)

# define the limits for the random number generator based on the errors in rm2, alpha and beta
    low_rm2 = rg_infinite_contrast_squared - rg_infinite_contrast_squared_error
    high_rm2 = rg_infinite_contrast_squared + rg_infinite_contrast_squared_error
    high_alpha = alpha + alpha_error
    low_alpha = alpha - alpha_error
    high_beta = -beta + beta_error
    low_beta = -beta - beta_error

# Big loop over ntrials:
    for j in range(ntrials):
        new_fit = numpy.zeros(3)
        B = numpy.zeros((3, 3))
        solution = []

# generate new values for the x coeffients based on the propagated errors; we want a new random number for each ntrial
        new_fit[0] = numpy.random.default_rng().uniform(low_rm2, high_rm2)
        new_fit[1] = numpy.random.default_rng().uniform(low_alpha, high_alpha)
        new_fit[2] = numpy.random.default_rng().uniform(low_beta, high_beta)

# get new mean values for R1, R2, D using delta_rho at all contrasts
# TODO: since the B matrix only changes as a function of contrast and not as a function of ntrials, maybe the entire matrix can be passed to this Monte Carlo method.  Or it can at least be taken out of the ntrials loop? Test to see if B is the same for each ntrial.
        for i in range(number_of_contrast_points):
            B[0][0] = volume_fraction[0]
            B[0][1] = volume_fraction[1]
            B[0][2] = volume_fraction[0]*volume_fraction[1]
            B[1][0] = (delta_rho[i][0] - delta_rho[i][1]) * \
                volume_fraction[0]*volume_fraction[1]
            B[1][1] = -(delta_rho[i][0] - delta_rho[i][1]) * \
                volume_fraction[0]*volume_fraction[1]
            B[1][2] = (delta_rho[i][0] - delta_rho[i][1])*volume_fraction[0]*volume_fraction[1] * \
                (volume_fraction[1]*volume_fraction[1] -
                 volume_fraction[0]*volume_fraction[0])
            B[2][0] = 0.0
            B[2][1] = 0.0
            B[2][2] = -(delta_rho[i][0] - delta_rho[i][1])*(delta_rho[i][0] - delta_rho[i][1]) * \
                volume_fraction[0]*volume_fraction[0] * \
                volume_fraction[1]*volume_fraction[1]

            this_solution = numpy.linalg.solve(B, new_fit)
            solution.append(this_solution)
#           print('solution: ', solution)

        r12r22d2 = numpy.array(solution)
        mean_r12r22d2 = numpy.mean(r12r22d2, axis=0)
        r1r2d = numpy.sqrt(r12r22d2)
        mean_r1r2d = numpy.mean(r1r2d, axis=0)
#       print('r12r22d: ', r12r22d2)
#       print('mean r12r22d2: ', mean_r12r22d2)
#       print('r1r2d: ', r1r2d)
#       print('mean r1r2d: ', mean_r1r2d)

# append mean values to master list for all ntrials
        r12r22d2_list[j+1] = mean_r12r22d2
        r1r2d_list[j+1] = mean_r1r2d

# end loop over ntrials

#    print(r12r22d2_list)
#    print(r1r2d_list)

    number_of_points = len(r1r2d_list)
#    print('number of points: ', number_of_points)

    return r12r22d2_list, r1r2d_list, number_of_points
