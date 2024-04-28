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
#       3/2024           --      added MC error analysis:   Susan Krueger
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

import time
import numpy
import scipy.optimize
import json

def save_data_to_plot_as_json(other_self, izero, izero_error, izero_calc, diff):

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    data_dict = {
        'fraction_d2o': [],
        'izero': [],
        'izero_error': [],
        'izero_calc': [],
        'I(0)-calculated_I(0)': []
    }

    for i in range(mvars.number_of_contrast_points):
    # Append new values to the lists in the dictionary
        data_dict['fraction_d2o'].append(mvars.fraction_d2o[i])
        data_dict['izero'].append(izero[i])
        data_dict['izero_error'].append(izero_error[i])
        data_dict['izero_calc'].append(izero_calc[i])
        data_dict['I(0)-calculated_I(0)'].append(diff[i])

    json_data = json.dumps(data_dict)

    mcavars.json_outfile.write(json_data)
    mcavars.json_outfile.close()

    return





def get_molecular_weights(other_self):
    r'''
    **Get Molecular Weights** is the **Stoichiometry Analysis** method that calculates the molecular weights of the components in a complex containing multiple copies of two components.

    Executed if stoichiometry_flag == True.

    This method calls **Molecular Weight Function**, **Get MC Errors** and **scipy.optimize.curve_fit.

    Notes:

        Output file with chosen output_file_name is stored in the output file directory.
        TODO: Need boilerplate language to add variable name/path for the directory. 

        Output file contains:
            - module variable input parameters
            - calculated molecular weights of the components in kDa
            - total molecular weight of the complex in kDa
            - weight fractions of each component
            - volume fractions of each component

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
    partial_specific_volume:  float array (dimension = number_of_components)
        partial specific volume of each component in cm\ :sup:`3`\ /g
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
    outfile: string
        output file name (with full path): path + output_file_name

    Returns
    -------
    molecular_weight_1: float
        average optimized molecular weight, :math:`(M_{w})_{1}`, for component 1 in kDa
    molecular_weight_1_error:  float
        standard error of the mean for molecular_weight_1
    molecular_weight_2: float
        average optimized molecular weight, :math:`(M_{w})_{2}`, for component 2 in kDa
    molecular_weight_2_error:  float
        standard error of the mean for molecular_weight_2
    total_molecular_weight:  float
        total molecular weight of the complex in kDa (from the averaged optimized molecular weights)
    total_molecular_weight_error:  float
        error on the total_molecular_weight propagated from the errors on molecular_weight_1 and molecular_weight_2
    weight_fraction_1:  float
        weight fraction of component 1
    weight_fraction_1_error:  float
        error on weight_fraction_1 propagated from the errors on molecular_weight_1 and total_molecular_weight
    weight_fraction_2:  float
        weight fraction of component 2
    weight_fraction_2_error:  float
        error on weight_fraction_2 propagated from the errors on molecular_weight_2 and total_molecular_weight
    volume_fraction_1:  float
        volume fraction of component 1
    volume_fraction_1_error:  float
        error on volume_fraction_1 propagated from the errors on the individual and total volumes
    volume_fraction_2:  float
        volume fraction of component 2
    volume_fraction_2_error:  float
        error on volume_fraction_2 propagated from the errors on the individual and total volumes

    Note
    ----

    - Optimized_molecular_weights are the optimal values for the molecular weights such that the sum of the squared error of molecular_weight_function(x, optimized_molecular_weights) - y is minimized. 

    - The y array = 0 since molecular_weight_function is defined to be equal to zero

    - The x array = the values of the coefficients that consist of the known parameters in molecular_weight_function once it is written in terms of the individual :math:`M_{w}` values.

    '''
# mvars used:  fraction_d2o, number_of_contrast_points, partial_specific_volume, izero, izero_error, concentration, concentration_error, delta_rho
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
    log.debug('partial_specific_volume_1,partial_specific_volume_2: ' +
              str(partial_specific_volume_1) + ',' + str(partial_specific_volume_2) + '\n')

    log.debug('number of contrast points: ' +
              str(mvars.number_of_contrast_points) + '\n')

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

    # print ('transposed x: ', x)
    # print ('length of transposed x: ', len(x))

    pgui('initial molecular weights calculation')
#    print('calculating molecular weights')
    optimized_molecular_weights, covariance = scipy.optimize.curve_fit(
        molecular_weight_function, x, y)

    log.debug('initial calculated molecular_weight_1,molecular_weight_2: ' +
         str(optimized_molecular_weights)+'\n')

# Now recalculate m1 and m2 using coefficients that randomly vary within their errors and save them in an array
# TODO:  Do we want number of MC trials to be a user input?

    pgui('Performing Monte Carlo error analysis')
    ntrials = 1000
    molecular_weights_list, number_of_points = get_mc_errors(
        optimized_molecular_weights, mvars.number_of_contrast_points, mvars.partial_specific_volume, mvars.izero, mvars.izero_error, mvars.concentration, mvars.concentration_error, mvars.delta_rho, ntrials)

    log.debug('molecular_weights_list: ' + str(molecular_weights_list) + '\n')
    number_of_points = len(molecular_weights_list)
    log.debug('number of points: ' + str(number_of_points) + '\n')

# Take average and find standard error of the mean
    pgui('averaging results')
    average_molecular_weight = numpy.mean(molecular_weights_list, axis=0)
    log.debug('average: ' + str(average_molecular_weight) + '\n')
    standard_deviation = numpy.std(molecular_weights_list, axis=0)
    log.debug('std deviation: ' + str(standard_deviation) + '\n')
    standard_error = standard_deviation/numpy.sqrt(number_of_points)
    log.debug('standard error of the mean: ' + str(standard_error) + '\n')
# If we want a 95% confidence interval, it is defined by 1.96*standard_error

    pgui('results written to output file: %s' %
         (mcavars.multi_component_analysis_path+mvars.output_file_name))
    pgui('-------------------------------')
    mcavars.outfile.write('--------------------------------\n')
    pgui('Final Results\n')
    mcavars.outfile.write('Final Results\n')

# molecular weights in kDa
    pgui('Errors are standard errors of the mean based on ' +
         str(number_of_points)+' molecular weight calculations\n')
    mcavars.outfile.write('Errors are standard errors of the mean based on ' +
                          str(number_of_points)+' molecular weight calculations\n')

    molecular_weight_1 = average_molecular_weight[0]*10**3
    molecular_weight_2 = average_molecular_weight[1]*10**3
    standard_error_1 = standard_error[0]*10**3
    standard_error_2 = standard_error[1]*10**3
    pgui('average molecular weight 1 (kDa): ' + ('%.3f' %
         molecular_weight_1) + ' +/- ' + ('%.3f' % standard_error_1) + '\n')
    pgui('average molecular weight 2 (kDa): ' + ('%.3f' %
         molecular_weight_2) + ' +/- ' + ('%.3f' % standard_error_2) + '\n')
    mcavars.outfile.write('average molecular weight 1 (kDa): ' + ('%.3f' %
                                                                  molecular_weight_1) + ' +/- ' + ('%.3f' % standard_error_1) + '\n')
    mcavars.outfile.write('average molecular weight 2 (kDa): ' + ('%.3f' %
                                                                  molecular_weight_2) + ' +/- ' + ('%.3f' % standard_error_2) + '\n')

# weight fractions
    total_molecular_weight = molecular_weight_1 + molecular_weight_2
    total_molecular_weight_error = numpy.sqrt(
        standard_error_1**2 + standard_error_2**2)
    pgui('Total Mw (kDa): ' + ('%.3f' % total_molecular_weight) +
         ' +/- ' + ('%.3f' % total_molecular_weight_error)+'\n')
    mcavars.outfile.write('Total Mw (kDa): ' + ('%.3f' % total_molecular_weight) +
                          ' +/- ' + ('%.3f' % total_molecular_weight_error)+'\n')

    weight_fraction_1 = molecular_weight_1/total_molecular_weight
    weight_fraction_2 = molecular_weight_2/total_molecular_weight
    weight_fraction_1_error = numpy.sqrt((standard_error_1/total_molecular_weight)**2 + (
        molecular_weight_1*total_molecular_weight_error/total_molecular_weight**2)**2)
    weight_fraction_2_error = numpy.sqrt((standard_error_2/total_molecular_weight)**2 + (
        molecular_weight_2*total_molecular_weight_error/total_molecular_weight**2)**2)

    pgui('weight fraction 1 (kDa): ' + ('%.3f' % weight_fraction_1) +
         ' +/- ' + ('%.3f' % weight_fraction_1_error) + '\n')
    pgui('weight fraction 2 (kDa): ' + ('%.3f' % weight_fraction_2) +
         ' +/- ' + ('%.3f' % weight_fraction_2_error) + '\n')
    mcavars.outfile.write('weight fraction 1 (kDa): ' + ('%.3f' % weight_fraction_1) +
                          ' +/- ' + ('%.3f' % weight_fraction_1_error) + '\n')
    mcavars.outfile.write('weight fraction 2 (kDa): ' + ('%.3f' % weight_fraction_2) +
                          ' +/- ' + ('%.3f' % weight_fraction_2_error) + '\n')

# volume fractions
# First, calculate the volumes using the mean Mw values. This will be used to propagate the errors in the volume fraction. Na was omitted since it will divide out when calculating the volume fraction.
    volume_1 = partial_specific_volume_1*molecular_weight_1
    volume_1_error = partial_specific_volume_1*standard_error_1
    volume_2 = partial_specific_volume_2*molecular_weight_2
    volume_2_error = partial_specific_volume_2*standard_error_2

    log.debug('volume_1, volume_2: ' + ('%.6f' %
              volume_1)+',\t'+('%.6f' % volume_2) + '\n')
    log.debug('volume_1_error, volume_2_error: ' + ('%.6f' %
              volume_1_error)+',\t'+('%.6f' % volume_2_error) + '\n,')

    total_volume = volume_1 + volume_2
    total_volume_error = numpy.sqrt(volume_1_error**2 + volume_2_error**2)
    log.debug('total volume: ' + str(total_volume) + '\n')
    log.debug('total volume error: ' + str(total_volume_error) + '\n')
    volume_fraction_1 = volume_1/(volume_1 + volume_2)
    volume_fraction_2 = volume_2/(volume_1 + volume_2)
    volume_fraction_1_error = numpy.sqrt(
        (volume_1_error/total_volume)**2 + (volume_1*total_volume_error/total_volume**2)**2)
    volume_fraction_2_error = numpy.sqrt(
        (volume_2_error/total_volume)**2 + (volume_2*total_volume_error/total_volume**2)**2)
    pgui('volume fraction 1 (kDa): '+('%.3f' % volume_fraction_1) +
         ' +/- '+('%.3f' % volume_fraction_1_error) + '\n')
    pgui('volume fraction 2 (kDa): ' + ('%.3f' % volume_fraction_2) +
         ' +/- '+('%.3f' % volume_fraction_2_error) + '\n')
    mcavars.outfile.write('volume fraction 1 (kDa): ' + ('%.3f' %
                          volume_fraction_1) + ' +/- ' + ('%.3f' % volume_fraction_1_error) + '\n')
    mcavars.outfile.write('volume fraction 2 (kDa): ' + ('%.3f' %
                          volume_fraction_2) + ' +/- ' + ('%.3f' % volume_fraction_2_error) + '\n')

# Calculate I(0) using the mean molecular weight values. Use the values before rescaling so I(0) will be in the correct units for comparison to the experimental values. Then, calculate the residuals for plotting and writing to file.
    print('calculating residuals')

    izero_calc = numpy.zeros(mvars.number_of_contrast_points)
    diff = numpy.zeros(mvars.number_of_contrast_points)

    for i in range(mvars.number_of_contrast_points):
        izero_calc[i] = (mvars.concentration[i]*(average_molecular_weight[0] + average_molecular_weight[1])/Na) *\
            (weight_fraction_1 * mvars.delta_rho[i][0]*partial_specific_volume_1 +
             weight_fraction_2*mvars.delta_rho[i][1] * partial_specific_volume_2)**2
        diff[i] = mvars.izero[i] - izero_calc[i]

    pgui('\n\nfraction_d2o        I(0)            I(0)_calc       I(0)-I(0)_calc\n')
    for i in range(mvars.number_of_contrast_points):
        pgui('%9.4f\t%9.4f\t%9.4f\t%9.4f' %
             (mvars.fraction_d2o[i], mvars.izero[i], izero_calc[i], diff[i]))
    mcavars.outfile.write(
        '\n\nfraction_d2o        I(0)            I(0)_calc       I(0)-I(0)_calc\n')
    for i in range(mvars.number_of_contrast_points):
        mcavars.outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\n' % (
            mvars.fraction_d2o[i], mvars.izero[i], izero_calc[i], diff[i]))

    mcavars.outfile.close()

    time.sleep(0.1)

    save_data_to_plot_as_json(other_self, mvars.izero, mvars.izero_error, izero_calc, diff)


    return


def molecular_weight_function(x, molecular_weight_1, molecular_weight_2):
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

# TODO: since the MC error analysis method isn't used by other methods, it can be rewritten to use mvars and mcavars


def get_mc_errors(optimized_molecular_weights, number_of_contrast_points, partial_specific_volume, izero, izero_error, concentration, concentration_error, delta_rho, ntrials):
    '''
    Method to perform Monte Carlo error analysis to take into account the errors on I(0) and concentration in the final values of the molecular weights.  The **Molecular Weight Function** is solved many times with different values of the coefficients constrained by +/- coefficient error.  The final molecular weights are the mean molecular weight values from all of the Monte Carlo trials, with errors calculated as the standard errors of the mean.

    This method calls numpy.random.default_rng().uniform and scipy.optimize.curve_fit

    Parameters
    ----------
    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    izero:  float array (dimension = number_of_contrast_points)
        I(0) value at each contrast in cm\ :sup:`-1`\ 
    izero_error:  float array (dimension = number_of_contrast_points)
        I(0) error value at each contrast
    concentration:  float array (dimension = number_of_contrast_points)
        concentration at each contrast in mg/mL
    concentration_error:  float array (dimension = number_of_contrast_points)
        concentration error at each contrast
    partial_specific_volume:  float array (dimension = number_of_components)
        partial specific volume of each component in cm\ :sup:`3`\ /g
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
    optimized_molecular_weights:  float array (dimensions = number_of_components)
        The initial calculated molecular weights for each component before the Monte Carlo error analysis

    Returns
    -------
    molecular_weights_list: 2D float array (dimension = ntrials+1 x number_of_components)
        optimized molecular weights, :math:`(M_{w})_{1}` and :math:`(M_{w})_{2}`, for components 1 and 2 in kDa; includes the initial optimized molecular weights and those from each Monte Carlo trial
    number_of_points: float
        the length of the molecular_weights_list array; used in determining the average optimzed molecular weights

    '''
# Big loop over ntrials:
#    print('Monte Carlo error analysis')
#    print('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')

    Na = 6.023  # 10**23
    partial_specific_volume_1 = partial_specific_volume[0]
    partial_specific_volume_2 = partial_specific_volume[1]
    molecular_weights_list = numpy.zeros((ntrials+1, 2))
    molecular_weights_list[0] = [
        optimized_molecular_weights[0], optimized_molecular_weights[1]]
    for j in range(ntrials):
        new_x = numpy.zeros(shape=(number_of_contrast_points, 4))
# Here we are defining the x coeficients at each contrast for the set of simultaneous I(0) equations and propagating the errors in I(0) and concentration.
        for i in range(number_of_contrast_points):
            delta_rho_1 = delta_rho[i][0]
            delta_rho_2 = delta_rho[i][1]
            izero_1 = izero[i]
            izero_1_error = izero_error[i]
            concentration_1 = concentration[i]
            concentration_1_error = concentration_error[i]
#            print('i, delta_rho_1, delta_rho_2, izero, izero_error, concentration, concentration_error: ', i, delta_rho_1, delta_rho_2, izero_1, izero_1_error, concentration_1, concentration_1_error)
            x0 = delta_rho_1**2*partial_specific_volume_1**2
            x1 = 2*delta_rho_1*delta_rho_2*partial_specific_volume_1*partial_specific_volume_2
            x2 = delta_rho_2**2*partial_specific_volume_2**2
            x3 = izero_1*Na/concentration_1
            x3_error = numpy.sqrt((Na*izero_1_error/concentration_1) **
                                  2 + (-Na*izero_1*concentration_1_error/concentration_1)**2)

# define the limits for the random number generator based on the propagated errors
            low_x3 = x3 - x3_error
            high_x3 = x3 + x3_error

# generate new values for the x coeffients based on the propagated errors
            new_x0 = x0
            new_x1 = x1
            new_x2 = x2
            new_x3 = numpy.random.default_rng().uniform(low_x3, high_x3)

            new_x[i] = (new_x0, new_x1, new_x2, new_x3)

#        print ('new_x: ', new_x)
#        print ('length of new_x: ', len(new_x))

        new_y = numpy.zeros(len(new_x))
#        print ('new_y: ', new_y)
#        print ('length of new_y: ', len(new_y))

    # transpose new_x to have correct inputs for curve_fit, which accepts a 1 column matrix with number of rows = number of contrast points.  We have a 1 row matrix of length = number of contrast points.
        new_x = new_x.T

#        print('transposed new x: ', new_x)
#        print ('length of transposed new x: ', len(new_x))

        new_optimized_molecular_weights, covariance = scipy.optimize.curve_fit(
            molecular_weight_function, new_x, new_y)

#        print('new calculated molecular_weight_1,molecular_weight_2: '+str(new_optimized_molecular_weights)+'\n')

        molecular_weights_list[j+1] = [new_optimized_molecular_weights[0],
                                       new_optimized_molecular_weights[1]]

# end of big loop

#    print('molecular_weights_list: ', molecular_weights_list)
    number_of_points = len(molecular_weights_list)
#    print('number of points: ', number_of_points)

    return molecular_weights_list, number_of_points
