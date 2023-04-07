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
#       MULTI COMPONENT ANALYSIS: STOICHIOMETRY
#
#       10/22/2021       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''

    **Stoichiometry** contains the **Get Molecular Weights** method and related functions that calculate the molecular weights of the components in a complex containing multiple copies of two components.
    
    **Inputs:**
    
    I(0), concentration, fraction D\ :sub:`2`\ O for each contrast
    
    contrast values (delta_rho) for each component at each contrast
    
    **Outputs:**
    
    optimized molecular weights, weight fractions and volume fractions for each component
    
    Called by the **Multi-component Analysis** module.
    
'''


import io
import os
import sys
import string
import locale
import time
import numpy
import scipy.optimize


def get_molecular_weights(other_self):
    r'''
    **Get Molecular Weights** is the **Stoichiometry Analysis** method that calculates the molecular weights of the components in a complex containing multiple copies of two components.

    Executed if stoichiometry_flag == True.

    This method calls **Molecular Weight Function**.

    Notes:

        Output file with chosen output_file_name is stored in the output file directory.
        TODO: Need boilerplate language to add variable name/path for the directory. 

        Output file contains:
            - module variable input parameters
            - calculated molecular weights of the components in kDa
            - total molecular weight of the complex in kDa
            - weight fractions of each component

    Parameters
    ----------
    run_name:  string
        Name of directory containing the outputs
    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    number_of_components: int
        The number of components in the molecule with different scattering length densities
    fraction_d2o:   float array (dimension = number_of_contrast_points)
        The fraction D\ :sub:`2`\ O values that define the contrasts
    izero:  float array (dimension = number_of_contrast_points)
        I(0) value at each contrast in cm\ :sup:`-1`\
    concentration:  float array (dimension = number_of_contrast_points)
        concentration at each contrast in mg/mL
    partial_specific_volume:  float array (dimension = number_of_components)
        partial specific volume of each component in cm\ :sup:`3`\ /g
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
    output_file_name:  string
        name of output file

    Returns
    -------
    molecular_weight_1: float
        optiminzed molecular weight, :math:`(M_{w})_{1}`, for component 1 in kDa
    molecular_weight_2: float
        optimized molecular weight, :math:`(M_{w})_{2}`, for component 2 in kDa
    total_molecular_weight:  float
        total molecular weight of the complex in kDa
    weight_fraction_1:  float
        weight fraction of component 1
    weight_fraction_2:  float
        weight fraction of component 2
    volume_fraction_1:  float
        volume fraction of component 1
    volume_fraction_2:  float
        volume fraction of component 2

    Note
    ----

    - Optimized_molecular_weights are the optimal values for the molecular weights such that the sum of the squared error of molecular_weight_function(x, optimized_molecular_weights) - y is minimized. 

    - The y array = 0 since molecular_weight_function is defined to be equal to zero

    - The x array = the values of the coefficients that consist of the known parameters in molecular_weight_function once it is written in terms of the individual :math:`M_{w}` values.

    - scipy.optimize.curve_fit is being used to solve the set of simultaneous equations using non-linear least squares fitting.  It uses scipy 1.7 under Python 3. No syntax changes were necessary when converting from scipy version 1.2.1 in Python 2.7.


    '''
# mvars used:  fraction_d2o, number_of_contrast_points, partial_specific_volume, izero, concentration, delta_rho
# mcavars used:  outfile

    log = other_self.log
    log.debug('in get_molecular_weights')
    pgui = other_self.run_utils.print_gui

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in range(60)])

    pgui('\n%s \n' % (st))
    pgui('DATA FROM RUN: %s \n\n' % (ttxt))

    Na = 6.023  # 10**23

    partial_specific_volume_1 = mvars.partial_specific_volume[0]
    partial_specific_volume_2 = mvars.partial_specific_volume[1]
#   print ('partial_specific_volume_1,partial_specific_volume_2: ', partial_specific_volume_1, partial_specific_volume_2)

#   print ("number of contrast points: ", number_of_contrast_points)

    pgui('Stoichiometry method\n')
    pgui('setting up I(0) equation coefficients')

    # if mvars.read_from_file == True:
    #    #pgui('\n')
    #    pgui('contrast values read from input file: %s' % (mvars.input_file_name))

    x = numpy.zeros(shape=(mvars.number_of_contrast_points, 4))

    # Here we are defining the x coeficients at each contrast for the set of simultaneous I(0) equations.
    for i in range(mvars.number_of_contrast_points):
        delta_rho_1 = mvars.delta_rho[i][0]
        delta_rho_2 = mvars.delta_rho[i][1]
        izero_1 = mvars.izero[i]
        concentration_1 = mvars.concentration[i]
#        print('delta_rho_1, delta_rho_2, izero, concentration: ', delta_rho_1, delta_rho_2, izero_1, concentration_1)
        x0 = delta_rho_1**2*partial_specific_volume_1**2
        x1 = 2*delta_rho_1*delta_rho_2*partial_specific_volume_1*partial_specific_volume_2
        x2 = delta_rho_2**2*partial_specific_volume_2**2
        x3 = izero_1*Na/concentration_1
        x[i] = (x0, x1, x2, x3)

    #print ('x: ', x)
    #print ('length of x: ', len(x))

    y = numpy.zeros(len(x))
    # print ('y: ', y)

    # transpose x to have correct inputs for curve_fit.
    x = x.T

    #print ('transposed x: ', x)
    #print ('length of transposed x: ', len(x))

    pgui('calculating molecular weights')

    optimized_molecular_weights, covariance = scipy.optimize.curve_fit(
        molecular_weight_function, x, y)

# TODO: Test some errors on the coefficients (hardwired) here to see if we should be allowing users to enter errors. Look at Watson et al power law paper for ranges of values for molecular volume. Would we ask for errors on each parameter, i.e., izero, delta_rho, concentration, etc. and then propagate them here when we calculate the coefficients? Curve_fit wants the errors in y, not the coefficients, x.  So what is the best way to test how the errors in x affect y?

#    pgui('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
#    mcavars.outfile.write('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
#    pgui('calculated covariance: '+str(covariance)+'\n')
#    mcavars.outfile.write('calculated covariance: '+str(covariance)+'\n')

    pgui('results written to output file: %s' %
         (mcavars.multi_component_analysis_path+mvars.output_file_name))
    pgui('-------------------------------')
    mcavars.outfile.write('--------------------------------\n')
    pgui('Final Results\n')
    mcavars.outfile.write('Final Results\n')

    # convert molecular weights to kDa
    molecular_weight_1 = optimized_molecular_weights[0]*10**3
    molecular_weight_2 = optimized_molecular_weights[1]*10**3
    pgui('molecular_weight_1, molecular_weight_2 (kDa): ' +
         str(molecular_weight_1)+'\t'+str(molecular_weight_2)+'\n')
    mcavars.outfile.write('molecular_weight_1, molecular_weight_2 (kDa): ' +
                          str(molecular_weight_1)+'\t'+str(molecular_weight_2)+'\n')


# Since we aren't taking into account the errors on the coefficients, the standard deviations shouldn't be reported?
# The residuals for I(0) are calcluated below.
    # The covariance matrix will be infinite if there are only two contrasts, so we don't want to calculate standard deviations in this case; numpy.all tests whether all array elements along a given axis evaluate to True.
#    if numpy.all(covariance != numpy.inf):

    # rescale values in covariance matrix since molecular weights were converted to kDa
#        rescaled_covariance = covariance*10**6
#        pgui('covariance matrix: ' + str(rescaled_covariance)+'\n')
#        mcavars.outfile.write('covariance matrix: ' + str(rescaled_covariance)+'\n')
#        # compute one standard deviation of the molecular weights
#        standard_deviation_1 = numpy.sqrt(rescaled_covariance[0][0])
#        standard_deviation_2 = numpy.sqrt(rescaled_covariance[1][1])
#        pgui('standard_deviation_1, standard_deviation_2 (kDa): ' +
#                  str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
#        mcavars.outfile.write('standard_deviation_1, standard_deviation_2 (kDa): ' +
#                  str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
#        correlation_coefficient = rescaled_covariance[0][1]/(standard_deviation_1*standard_deviation_2)
#        pgui('correlation coefficient: '+str(correlation_coefficient)+'\n')
#        mcavars.outfile.write('correlation coefficient: '+str(correlation_coefficient)+'\n')

    total_molecular_weight = molecular_weight_1 + molecular_weight_2
    weight_fraction_1 = molecular_weight_1/total_molecular_weight
    weight_fraction_2 = molecular_weight_2/total_molecular_weight

# Calculate volume fractions in case they are needed in other methods.

    volume_fraction_1 = partial_specific_volume_1*molecular_weight_1 / \
        (partial_specific_volume_1*molecular_weight_1 +
         partial_specific_volume_2*molecular_weight_2)
    volume_fraction_2 = partial_specific_volume_2*molecular_weight_2 / \
        (partial_specific_volume_1*molecular_weight_1 +
         partial_specific_volume_2*molecular_weight_2)

    pgui('Total Mw (kDa): '+str(total_molecular_weight))
    pgui('weight_fraction_1, weight_fraction_2: ' +
         str(weight_fraction_1)+'\t'+str(weight_fraction_2))
    pgui('volume_fraction_1, volume_fraction_2: ' +
         str(volume_fraction_1)+'\t'+str(volume_fraction_2))
    mcavars.outfile.write('Total Mw (kDa): '+str(total_molecular_weight)+'\n')
    mcavars.outfile.write('weight_fraction_1, weight_fraction_2: ' +
                          str(weight_fraction_1)+'\t'+str(weight_fraction_2)+'\n')
    mcavars.outfile.write('volume_fraction_1, volume_fraction_2: ' +
                          str(volume_fraction_1)+'\t'+str(volume_fraction_2)+'\n\n')

# Calculate I(0) using the optimized molecular weight values. Use the values before rescaling so I(0) will be in the correct units. Then, calculate the residuals for plotting and writing to file.

    izero_calc = numpy.zeros(mvars.number_of_contrast_points)
    diff = numpy.zeros(mvars.number_of_contrast_points)

    for i in range(mvars.number_of_contrast_points):
        izero_calc[i] = (mvars.concentration[i]*(optimized_molecular_weights[0] + optimized_molecular_weights[1])/Na)*(weight_fraction_1 *
                                                                                                                       mvars.delta_rho[i][0]*partial_specific_volume_1 + weight_fraction_2*mvars.delta_rho[i][1]*partial_specific_volume_2)**2
        diff[i] = mvars.izero[i] - izero_calc[i]

        #print('i, I(0), I(0)calc, diff: ', mvars.izero[i], izero_calc[i], diff[i])

    mcavars.outfile.write(
        'fraction_d2o        I(0)            I(0)_calc       I(0)-I(0)_calc\n')
    for i in range(mvars.number_of_contrast_points):
        mcavars.outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\n' % (
            mvars.fraction_d2o[i], mvars.izero[i], izero_calc[i], diff[i]))

    mcavars.outfile.close()

    return


def molecular_weight_function(x, molecular_weight_1, molecular_weight_2):
    # TODO: generalize this function any number of components? Or offer hardwired options for 2 and 3 components?
    r'''The **Molecular Weight Function**. Called from **Get Molecular Weights**. 

    Note
    ----

    Molecular Weight Function is the equation: 

    :math:`n(\sum_{i} \Delta  \rho_{i}V_{i})^{2} - I(0) = 0`, where 

    :math:`I(0)` is the intensity at :math:`q=0`

    :math:`\Delta \rho_{i}` is the contrast of :math:`i^{\text{th}}` component

    :math:`V_{i}` is the volume of :math:`i^{\text{th}}` component

    :math:`n = c\frac{N_{A}}{\sum_{i} (M_{w})_{i}}` is the number density of particles where

    :math:`c` is the concentration

    :math:`N_{A}` is Avogadro's number and

    :math:`\sum_{i} (M_{w})_{i}` is the total molecular weight of the complex.

    If the volume is rewritten in terms of the partial specific volume, :math:`\overline{v}`,

    :math:`V = \overline{v} \frac {M_{w}}{N_{A}}`,

    then for **two components** :math:`M_{w} = (M_{w})_{1} + (M_{w})_{2}`, then the equation can be rewritten in terms of the individual :math:`M_{w}` values as

    :math:`x_{0} (M_{w})_{1}^2 + x_{1} (M_{w})_{1} (M_{w})_{2} + x_{2} (M_{w})_{2}^2 - x_{3} (M_{w})_{1} - x_{3} (M_{w})_{2} = 0`, where

    :math:`x_{0} = (\Delta \rho)_{1}^{2} \overline{v}_{1}^{2}`

    :math:`x_{1} = 2 (\Delta \rho)_{1} (\Delta \rho)_{2} \overline{v}_{1} \overline{v}_{2}`

    :math:`x_{2} = (\Delta \rho)_{2}^{2} \overline{v}_{2}^{2}`

    :math:`x_{3} = \frac {I(0) N_{A}}{c}`

    Parameters
    ----------
    x:  float array (dimension = 4 x number_of_contrast_points)
        :math:`x` coefficients at each contrast (currently hardwired for 2 components)

    Returns
    -------
    molecular_weight_1: float
        optiminzed molecular weight, :math:`(M_{w})_{1}`, for component 1
    molecular_weight_2: float
        optimized molecular weight, :math:`(M_{w})_{2}`, for component 2

    '''
    f = x[0] * molecular_weight_1**2 + x[1] * molecular_weight_1*molecular_weight_2 + \
        x[2] * molecular_weight_2**2 - x[3] * \
        molecular_weight_1 - x[3] * molecular_weight_2
    return f
