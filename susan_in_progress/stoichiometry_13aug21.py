# -*- coding: utf-8 -*-

"""
Stoichiometry: Method to calculate the Mw of the components in a two-component system for which a contrast variation series of measurements has been performed. It should be used for systems that have multiple copies of the two components. While this module can be used for experiment planning, it is more likely that it would be used for data analysis.

TODO: LIST methods in file and what they do. Top to bottom. Resembles the flow chart. 
Inputs and outputs (what does it return, what does it do?)  Look at zazmol as an example. 

delta_rho (contrast) is in units of 10**10 cm**-2
concentration is in units of 10^-3 g/cm**3 (mg/ml)
izero is in cm^-1
partial_specfic_volume is in units of cm**3/g
molecular_weight_1 and molecular_weight_2 are given in Da * 10^-6 and are converted to kDa
Avogadro's number, Na, is in units of 10**23

@author: susan.krueger
"""

# TODO: hard-wired for two components at this time; adding components would require changing the function definition on the fly, which could be done by going back to the basic I(0) equation, which is a "square of the sum"
# autopep8 run 3aug21 


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import io
import os
import string
import locale
import numpy
import scipy.optimize
import read_contrast_output_files

def get_molecular_weights(read_from_file, number_of_contrast_points, delta_rho, partial_specific_volume, izero, concentration, file_name_1, file_name):

    # TODO: this function could be outside of this method. If keeping it inside, put it right above the curve fit. NOTE:  I tried the latter and it didn't work.  Error said that molecular_weight_1 was referenced before assignment.
    # This function is the I(0) equation for two components rearranged such that the x array represents the known coefficients of the equation and the y array is == 0. This function is satisfied at each contrast.
    #TODO: put the actual equation in the doc string
    def izero_function(x, molecular_weight_1, molecular_weight_2):
        return x[0] * molecular_weight_1**2 + x[1] * molecular_weight_1*molecular_weight_2 + x[2] * molecular_weight_2**2 - x[3] * molecular_weight_1 - x[3] * molecular_weight_2


    Na = 6.023

    partial_specific_volume_1 = locale.atof(partial_specific_volume[0])
    partial_specific_volume_2 = locale.atof(partial_specific_volume[1])
    print ('partial_specific_volume_1,partial_specific_volume_2: ', partial_specific_volume_1, partial_specific_volume_2)

    # print ("number of contrast points: ", number_of_contrast_points)

    outfile = io.open(file_name, 'w')
    if read_from_file == 1:
        print('input file: ', file_name_1)
        # this is just to write the contrast calculator filename into the output file so we know this option was used
        outfile.write('input file: '+file_name_1+'\n')
    else:
        print('input file: None')
        outfile.write('input file: None\n')

    x = numpy.zeros(shape=(number_of_contrast_points, 4))

    # Here we are defining the x coeficients at each contrast for the set of simultaneous I(0) equations
    # NOTE that range is being used here since xrange in Python 2 == range in Python 3. So, we aren't doing exactly the same thing here as will be done in Python 3, as range will behave like Python 2's xrange in Python 3.  But, for this case, it shouldn't matter (if I understand correctly).
    for i in range(number_of_contrast_points):
        delta_rho_1 = locale.atof(delta_rho[i][0])
        delta_rho_2 = locale.atof(delta_rho[i][1])
        izero_1 = locale.atof(izero[i])
        concentration_1 = locale.atof(concentration[i])
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

    # optimized_molecular_weights are the optimal values for the molecular weights so that the sum of the squared error of izero_function(x, *optimized_molecular_weights) - y is minimized
    # covariance is the estimated covariance of molecular weights. The diagonals are used to obtained the standard deviations (below).
    optimized_molecular_weights, covariance = scipy.optimize.curve_fit(izero_function, x, y)

    print ('calculated molecular_weight_1,molecular_weight_2: ', optimized_molecular_weights)
    outfile.write('calculated molecular_weight_1,molecular_weight_2: '+str(optimized_molecular_weights)+'\n')
    print ('calculated covariance: ', covariance)
    outfile.write('calculated covariance: '+str(covariance)+'\n')

    print ('-------------------------------')
    outfile.write('--------------------------------\n')
    print ('Final Results:')
    outfile.write('Final Results\n')

    # convert molecular weights to kDa
    molecular_weight_1 = optimized_molecular_weights[0]*10**3
    molecular_weight_2 = optimized_molecular_weights[1]*10**3
    print ('molecular_weight_1, molecular_weight_2 (kDa): ', molecular_weight_1, molecular_weight_2)
    outfile.write('molecular_weight_1, molecular_weight_2 (kDa): '+str(molecular_weight_1)+'\t'+str(molecular_weight_2)+'\n')

    # The covariance matrix will be infinite if there are only two contrasts, so we don't want to calculate it in this case; numpy.all tests whether all array elements along a given axis evaluate to True.
    if numpy.all(covariance != numpy.inf):

        # rescale values in covariance matrix since molecular weights were converted to kDa
        rescaled_covariance = covariance*10**6
        print ('covariance matrix: ', rescaled_covariance)
        outfile.write('covariance matrix: ' + str(rescaled_covariance)+'\n')
        # compute one standard deviation of the molecular weights
        standard_deviation_1 = numpy.sqrt(rescaled_covariance[0][0])
        standard_deviation_2 = numpy.sqrt(rescaled_covariance[1][1])
        print ('standard deviations, standard_deviation_1, standard_deviation_2: ', standard_deviation_1, standard_deviation_2)
        outfile.write('standard deviations, standard_deviation_1, standard_deviation_2: ' +
                      str(standard_deviation_1)+'\t'+str(standard_deviation_2)+'\n')
        correlation = rescaled_covariance[0][1]/(standard_deviation_1*standard_deviation_2)
        print ('correlation: ', correlation)
        outfile.write('correlation: '+str(correlation)+'\n')

    total_molecular_weight = molecular_weight_1 + molecular_weight_2
    weight_fraction_1 = molecular_weight_1/total_molecular_weight
    weight_fraction_2 = molecular_weight_2/total_molecular_weight

    # Do we want to output the covariance matrix? How much meaning do the standard deviations have in this case since we are performing an "unweighted" fit? We aren't specifying errors on x, as they are hard to determine. While errors on I(0) and concentration can be estimated, it is harder to estimate the errors on delta rho and partial specific volume. 
    # Do we want to calculate volume fractions here as well? Will they be needed?
    # TODO: We will want to write all of the input values to the output file so that we know what values were used in the calculation.

    print ('Mw (kDa), weight_fraction_1, weight_fraction_2: ', total_molecular_weight, weight_fraction_1, weight_fraction_2)
    outfile.write('Mw (kDa), weight_fraction_1, weight_fraction_2: '+str(total_molecular_weight)+'\t'+str(weight_fraction_1)+'\t'+str(weight_fraction_2))
    outfile.close()


if __name__ == "__main__":

    # There are some inputs that will always be needed:
    #   1) number of contrasts (integer >/= 2)
    #   2) partial specific volumes (real, array)(We can make this part of the contrast calculator output; for now I am assuming the values will be input manually; we can provide default values for protein, RNA, DNA, but this will require similar input as for the contrast calculator where we ask whether the components are protein, DNA, RNA or other molecule)
    #   3) I(0) values (real, array)
    #   4) total concentration of complex (real array)
    # Users will have a choice for other inputs.
    #   1) Get contrasts (and perhaps partial specific volumes) from contrast calculator output file
    #       a) %D2O (real, array) (for consistency with contrast calculator and SasCalc; converted to fraction D2O)
    #       b) contrast calculator output filename  and read contrast values
    #   2) Input contrasts
    #       a) contrast (real, array)
    #       b) %D2O (real, array) (not really needed but maybe would help the user if they reattached to the job later)

    # right now, we will just read the values from an input file that lists the fraction D2O and the contrasts in the same way as the contrast calculator output file OR input the everything by hand.

    # SASSIE-web inputs
    # 0 is no and 1 is yes to reading the contrasts from a contrast calculator output file
    read_from_file = 1
    # must be >= number of components (2 in this case)
    number_of_contrast_points = 3
    number_of_components = 2
    fraction_d2o = [u'0.99', u'0.12', u'0.41']  # 1 values for each contrast
    partial_specific_volume = [u'0.745', u'0.903']  # 1 value for each component
    izero = [u'11.8', u'0.6', u'0.17']  # 1 value for each contrast
    concentration = [u'3.7', u'3.6', u'3.1']  # 1 value for each contrast
    path = ('./')
    # output file (will always be specified) We may only want to ask for a file prefix like in contrast calculator?
    file_name = os.path.join(path, '99_12_41.out')

    print ('Output file: ', file_name)

    # Specify whether contrast values will be read from a file.
    if read_from_file == 1:
        # input file (specify only if reading contrasts from file; default file can be '')
        file_name_1 = os.path.join(path, 'input_contrast.txt')
        delta_rho = read_contrast_output_files.read_contrast_file(
            file_name_1, fraction_d2o, number_of_contrast_points, number_of_components)
    else:
        file_name_1 = ''  # won't be used but needs to be defined
        # 2 values for each contrast point (since there are  2 components)
        delta_rho = [[u'-3.2', u'-5.7'], [u'1.6', u'0.26'], [u'0.031', u'-1.74']]

    print ('input file: ', file_name_1)

    # print ('delta_rho: ', delta_rho)

    get_molecular_weights(read_from_file, number_of_contrast_points,
            delta_rho, partial_specific_volume, izero, concentration, file_name_1, file_name)
