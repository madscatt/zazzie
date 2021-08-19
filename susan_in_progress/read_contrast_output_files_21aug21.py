# -*- coding: utf-8 -*-

"""
Read Contrast Output Files:  Reads the output files from the contrast calculator to obtain contrast, Mw, I(0) or other information for further calculation.  It can be used by both experiment planning and data analysis tools.

@author: susan.krueger
"""

# Right now, I am only reading contrast values from the *_contrast.txt output files in order to obtain contrast values at a particular %D2O.  Values MUST BE READ the same order as the other variables are entered. I am not sure how we are going to do this in SASSIE-web.  We need some way to know which partial specific volume, concentration, I(0), etc. goes with which contrast value.  In contrast calculator, we may need to ask users to input a name for the component along with the sequence or chemical formula so that this metadata can be used when info needs to be read from the files.

# autopep8 run 6aug21

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import io
import os
import sys
import string
import locale
import numpy
import scipy.stats



def read_contrast_file(file_name_1, fraction_d2o, number_of_contrast_points, number_of_components):

    # This will be in the input filter for the module.
    if(number_of_contrast_points == len(fraction_d2o)):
        print ('number of contrasts: ', number_of_contrast_points)
    else:
        print ('number of contrasts is not equal to length of fraction_d2o')
        sys.exit()

    if(number_of_contrast_points < number_of_components):
        print ('number of contrast points must be >= number of components')
        sys.exit()

    with io.open(file_name_1, 'r') as infile:
        lines = infile.readlines()
    # infile.close()


    # Read the fraction_d2o values of interest and associate them with numpy array x (fraction D2O). so don't use x
    x = numpy.zeros(shape=(number_of_contrast_points))
    for i in range(number_of_contrast_points):
        x[i] = locale.atof(fraction_d2o[i])
    # print ('frac D2O: ', x)

    number_of_lines = len(lines)
    # print ("number of lines ", number_of_lines)

# Find linear equation to calculate desired delta_rho values. This equation can eventually be part of the contrast calculator output file so it can be read directly instead of having to re-calculate it.

# First, find out where the desired delta_rho values are in the file.
    for i in range(number_of_lines):
        # print(i, lines[i][0:4])
        if(lines[i][0:4] == u'# fr'):
            starting_line = i
            # print('start here: ', i)

    # print(lines[starting_line])
    l = string.split(lines[starting_line], ',')
    # print(l)

# Associate the contrast arrays with a title [protein, RNA, DNA, Molecule 1, etc.], as we may have to ask users which contrast goes first to be in agreement with the way they are inputing I(0), concentration, etc.
    title = []
    title.append(u' frac D2O')
    for i in range(number_of_components+1):
        title.append(l[i+1])
    # print(title)
    # print(title[0])
    # print(title[1])
    # print(title[2])


# Read the lines and associate them with numpy array data [[fraction d2o from file] [contrast1] [contrast2]]

    number_of_data_points = number_of_lines - starting_line - 1
    # print('number of data points for linear fit: ', number_of_data_points)

    data = numpy.zeros(shape=(number_of_components + 1, number_of_data_points))
    # print('data: ', data)

    for i in range(starting_line+1, number_of_lines):
        # print('i: ', i)
        l = string.split(lines[i], '\t')
        # print(l)
        # print(l[0])
        # print(l[1])
        # print(l[2])
        for j in range(number_of_components+1):
            # print('j,line(j): ',j, l[j])
            data[j][i-starting_line-1] = locale.atof(l[j])
            # print(data[j][i-starting_line-1])

    # print('data: ', data)
    # print('data[0]: ', data[0])
    # print('data[1]: ', data[1])
    # print('data[0][0]: ', data[0][0])


# Find the slope and intercept for each component and put them into arrays. (I'm doing it this way in case we end up just reading the values from an equation that is included in the contrast calculator output file.)
# Match points are printed as a check.
    slope_values = numpy.zeros(shape=number_of_components)
    intercept_values = numpy.zeros(shape=number_of_components)
    # print('slope_values: ', slope_values)
    # print('intercept_values: ', intercept_values)
    # for mass in mass_array:  auto_ranges, only use range when you are skipping or need to start at a certain number, etc.
    for i in range(number_of_components):
        slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(
            data[0], data[i+1])
        x_intercept = -intercept / slope
        # print ('slope, intercept, x_intercept: ', slope, intercept, x_intercept)
        match_point = x_intercept * 100
        print('match point '+title[i+1]+': ', match_point)
        slope_values[i] = slope
        intercept_values[i] = intercept
    # print('slope_values: ', slope_values)
    # print('intercept_values: ', intercept_values)

# Find the delta_rho values and put them into a delta_rho array
# Converted the delta_rho values into list of lists with unicode values to pass to other methods, i.e., delta_rho = [[u'value', u'value'],[u'value',u'value'],[u'value',u'value]]
# The unicode conversion not be necessary for Python 3? Will it be unicode by default if just str is used?
    delta_rho = [[] for n in range(number_of_contrast_points)]
    # print('delta_rho: ', delta_rho)
    for i in range(number_of_contrast_points):
        for j in range(number_of_components):
            this_delta_rho = '%.3f' % (slope_values[j]*x[i] + intercept_values[j])
            unicode_delta_rho = unicode(this_delta_rho)
            # print(unicode_delta_rho,type(unicode_delta_rho))
            delta_rho[i].append(unicode_delta_rho)
    # print('delta_rho: ', delta_rho)

    return delta_rho


if __name__ == "__main__":

    path = ('./')
    file_name_1 = os.path.join(path, 'input_contrast.txt')
    print (file_name_1)
    number_of_contrast_points = 3
    number_of_components = 2
    fraction_d2o = [u'0.99', u'0.12', u'0.41']

    delta_rho = read_contrast_file(
        file_name_1, fraction_d2o, number_of_contrast_points, number_of_components)

    print('delta_rho after return to main: ', delta_rho)
