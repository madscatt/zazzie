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
#       CALCULATE MATCH POINT is the method that calculates a contrast
#       match point from SANS contrast variation data
#
#       10/5/2021       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''


    INPUTS:
      
        
    
    OUTPUTS:
       

    Requires line_fit
        
'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import io
import numpy
import sys
import linear_fit as linear_fit


def calculate_match_point(izero, izero_error, concentration, concentration_error, fraction_d2o, number_of_data_points, initial_match_point_guess, mode):

#   calculate the match point from SANS contrast variation data

#   In most cases, this routine will only be called when analyzing data (since planning is done with the contrast calculator), so it is assumed that the user has input I(0), c and an error on both I(0) and c, and a weighted fit will be used.  Values of "0" for the errors aren't allowed, so this will have to be handled in the module filter.

# calculate sqrt[I(0)/c] and propagate the error on I(0)/c to make the relationship between I(0) and fraction D2O linear.  The user needs to to input a frac_d2o value for which the sign changes. In other words, they will need to make an initial guess for the match point for their data and may have to try a few times before they get a good fit.  An initial guess should be able to be made just from the data (perhaps with help from contrast calculator output).

    normalized_izero = numpy.zeros(number_of_data_points, numpy.float)
    normalized_izero_error = numpy.zeros(number_of_data_points, numpy.float)
    square_root_izero = numpy.zeros(number_of_data_points, numpy.float)
    square_root_izero_error = numpy.zeros(number_of_data_points, numpy.float)
 
#    print('fraction d2o, izero, error: ', fraction_d2o, izero, izero_error)
    
    for i in range(number_of_data_points):
        normalized_izero[i] = izero[i]/concentration[i]
        normalized_izero_error[i] = numpy.sqrt(izero_error[i]**2/concentration[i]**2 + izero[i]**2*concentration_error[i]**2/concentration[i]**4)
#To obtain a straight line, set sqrt(I(0)/c) as negative if fraction D2O is <= to the initial matchpoint guess
        if (fraction_d2o[i] <= initial_match_point_guess):
#            print('fraction d2o is less than guess: ', fraction_d2o[i])
            square_root_izero[i] = numpy.sqrt(normalized_izero[i])
        else:
            square_root_izero[i] = -numpy.sqrt(normalized_izero[i])
        square_root_izero_error[i] = normalized_izero_error[i]/(2.0*numpy.sqrt(normalized_izero[i]))

#    print('normalized I(0), error: ', normalized_izero, normalized_izero_error)
#    print('sqrt(I(0)), error: ', square_root_izero, square_root_izero_error)


# calculate the linear fit of sqrt[I(0)] vs fraction D2O

    slope, slope_error, intercept, intercept_error, r_value, reduced_chi_squared, ycalc, diff = linear_fit.line_fit(fraction_d2o,square_root_izero, square_root_izero_error,number_of_data_points, mode)

    print('slope, error, intercept, error, r_value, reduced_chi_squared: ', slope, slope_error, intercept, intercept_error, r_value, reduced_chi_squared)

# calculate the match point from the slope and intercept

    match_point = -intercept/slope
    match_point_error = numpy.sqrt((intercept_error/slope)**2 + ((slope_error*intercept)/(slope**2))**2)

    print('match point, error: ', match_point, match_point_error)
    

    return match_point, match_point_error, square_root_izero, square_root_izero_error, ycalc, diff, r_value, reduced_chi_squared


if __name__ == "__main__":

    import matplotlib.pyplot as plt 
    path = ('./')
    input_file_name = os.path.join(path,'test_matchpoint1.txt')
    print (input_file_name)
    output_file_name = os.path.join(path,'matchpoint_testdata.out')
    initial_match_point_guess = 0.4
    mode = 1


    with io.open(input_file_name, 'r') as infile:
        line = infile.readlines()
#    print(a)
    number_of_data_points = len(line)
    print ("number of data points: ", number_of_data_points)

#split the lines into words
    words = []
    for i in range(number_of_data_points):
        value = line[i].split()
#        print(value)
        words.append(value)
#    print(words)

    fraction_d2o = numpy.zeros(number_of_data_points, numpy.float)
    izero = numpy.zeros(number_of_data_points, numpy.float)
    izero_error = numpy.zeros(number_of_data_points, numpy.float)
    concentration = numpy.zeros(number_of_data_points, numpy.float)
    concentration_error = numpy.zeros(number_of_data_points, numpy.float)

    
    for i in range(number_of_data_points):
        fraction_d2o[i] = words[i][0]
        izero[i] = words[i][1]
        izero_error[i] = words[i][2]
        if(i == 3):
            concentration[i] = 26.9
            concentration_error[i] = 1.3
        else:
            concentration[i] = 11.9
            concentration_error[i] = 0.6
        print('i, fraction_d2o[i], izero[i], izero_error[i], concentration, concentration_error: ', i, fraction_d2o[i], izero[i], izero_error[i], concentration[i], concentration_error[i])
        
    match_point, match_point_error, square_root_izero, square_root_izero_error, ycalc, diff, r_value, reduced_chi_squared = calculate_match_point(izero, izero_error, concentration, concentration_error, fraction_d2o, number_of_data_points, initial_match_point_guess, mode)

    outfile = io.open(output_file_name, 'w')
    outfile.write('input file: ' + input_file_name +'\n')
    outfile.write('number of points fit: ' + str(number_of_data_points) + '\n')
    outfile.write('match point: ' + str(round(match_point,4)) + ' +/- ' + str(round(match_point_error,4)) + '\n')
    outfile.write('r value: ' + str(round(r_value,4)) + '\n')    
    outfile.write('reduced chi-squared: ' + str(round(reduced_chi_squared,4)) + '\n')
    outfile.write('fraction_d2o  sqrt[I(0)/c]  sqrt[I(0)/c]_error  sqrt[I(0)/c]_calc  sqrt[I(0)/c]-sqrt[I(0)/c]_calc\n')
    for i in range(number_of_data_points):
        outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\t\t%9.4f\n' % (fraction_d2o[i], square_root_izero[i], square_root_izero[i], ycalc[i], diff[i]))
    outfile.close()

    fig, (ax0,ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.errorbar(fraction_d2o, square_root_izero, yerr=square_root_izero_error, fmt='bo', markersize = 2)
    ax0.plot(fraction_d2o,ycalc,'r--')
    ax0.set_title('data')
    ax0.set_ylabel('sqrt[I(0)/c]')
#
    ax1.plot(fraction_d2o, diff, 'ko', markersize = 5)
    ax1.set_title('residuals')
    ax1.set_xlabel('fraction D2O')
    ax1.set_ylim([-0.5, 0.5])
    plt.show()


