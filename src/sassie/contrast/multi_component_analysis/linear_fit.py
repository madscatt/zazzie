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
#       LINE FIT
#
#       10/4/2021       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Linear Fit** contains the method that performs a weighted linear fit to a straight line using the method of determinants. Adapted from the original FORTRAN program from Data Analysis and Error Reduction for the Physical Sciences, by Philip R. Bevington and D. Keith Robinson, McGraw-Hill, New York, NY 2003, 1992, 1969 (Chapter 6).

    **Inputs:**
    
        mode, number of data points, x array, y array, error in y array

    **Outputs:**
    
        slope, error in slope, intercept, error in intercept, calculated y values (ycalc), difference between y and ycalc, R value, reduced chi-squared
        
        
    Although a fit with weighting = 1.0/yerr\ :sup:`2`  is assumed here, mode is an
    input variable. Thus, this method is applicable to fits with different
    weightings if desired.
    
    Can be called by any method performing a fit to a straight line.
    
'''

import os
import io
import numpy
import sys


def line_fit(x, y, yerr, number_of_data_points, mode):
    '''
    Performs a weighted fit to a straight line :math:`y = a + bx`.

    Note
    ----

    The value of mode determines the method of weighting for the fit:

        1:  instrumental weighting = 1.0/yerr\ :sup:`2`

        0:  no weighting = 1.0

        -1:  statistical weighting = 1.0/abs(y)

    Parameters
    ----------

    mode: int (values = -1, 0, 1)
        determines the type of weighting for the fit
    number_of_data_points:  int
        number of data points    
    X: float array (dimension = number_of_data_points)
        x data
    Y: float array (dimension = number_of_data_points)
        y data
    yerr:float array (dimension = number_of_data_points)
        error in y data

    Returns
    -------

    slope: float
        slope of the straight line
    slope_error: float
        error in the slope
    intercept: float
        intercept of the straight line
    intercept_error: float
        error in the intercept
    r_value: float
        the R value
    reduced chi_squared: float
        chi-squared divided by: the number of data points - 2 degrees of freedom
    ycalc: float array (dimension = number_of_data_points)
        the y values calculated from the fitted straight line
    diff: float array (dimension = number_of_data_points)
        y - ycalc at each x value

    '''

#   fit to y = a + bx; b = slope, a = intercept

#    print('in line fit\n')
#    print(x, y, yerr)

    weight = numpy.zeros(number_of_data_points, numpy.float)
    ycalc = numpy.zeros(number_of_data_points, numpy.float)
    diff = numpy.zeros(number_of_data_points, numpy.float)

    for i in range(number_of_data_points):
        if mode == -1:
            if y[i] != 0:
                weight[i] = 1.0/numpy.abs(y[i])
            elif y[i] == 0:
                weight[i] = 1.0
        elif mode == 0:
            weight[i] = 1.0
        elif mode == 1:
            weight[i] = 1.0/yerr[i]**2
# accumulate weighted sums
    sumwt = numpy.sum(weight)
    sumx = numpy.sum(weight*x)
    sumy = numpy.sum(weight*y)
    sumx2 = numpy.sum(weight*x*x)
    sumy2 = numpy.sum(weight*y*y)
    sumxy = numpy.sum(weight*x*y)

#    print('sumwt, sumx, sumy, sumx2, sumy2, sumxy: ', sumwt, sumx, sumy, sumx2, sumy2, sumxy)

# calculate determinant
    determinant = sumwt*sumx2 - sumx*sumx
    if determinant > 0:
        intercept = (sumx2*sumy - sumx*sumxy)/determinant
        slope = (sumxy*sumwt - sumx*sumy)/determinant
    else:
        #       how to handle this error in the GUI
        sys.exit('error: determinant < or = 0 in line_fit')

# calculate slope, intercept, errors, R value, reduced chi-squared
    if mode == -1:
        variance = 1.0
    elif mode == 0:
        variance = (sumy2 + intercept*intercept*sumwt + slope*slope*sumx2 - 2.00 *
                    (intercept*sumy + slope*sumxy - intercept*slope*sumx))/(number_of_data_points - 2)
    elif mode == 1:
        variance = 1.0
    intercept_error = numpy.sqrt(variance*sumx2/determinant)
    slope_error = numpy.sqrt(variance*sumwt/determinant)
    r_value = (sumwt*sumxy - sumx*sumy) / \
        numpy.sqrt(determinant*(sumwt*sumy2 - sumy*sumy))

    for i in range(number_of_data_points):
        ycalc[i] = intercept + slope*x[i]
        diff[i] = y[i] - ycalc[i]
#        print('i, ycalc[i], diff[i]: ', i, ycalc[i], diff[i])

    chi_squared = numpy.sum(weight*diff*diff)
#    print('chi_squared: ', chi_squared)

#   reduced x2 assuming 2 degrees of freedom in the linear fit
    reduced_chi_squared = chi_squared/(number_of_data_points - 2)
#    print('reduced_chi_squared: ', reduced_chi_squared)

    return slope, slope_error, intercept, intercept_error, r_value, reduced_chi_squared, ycalc, diff


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import scipy.optimize
    path = ('./')
    input_file_name = os.path.join(path, 'test_linefit3.txt')
    print(input_file_name)
    output_file_name = os.path.join(path, 'linefit3.out')
    mode = 1

    with io.open(input_file_name, 'r') as infile:
        line = infile.readlines()
#    print(a)
    number_of_data_points = len(line)
    print("number of data points: ", number_of_data_points)

# split the lines into words
    words = []
    for i in range(number_of_data_points):
        value = line[i].split()
        print(value)
        words.append(value)
#    print(words)

    x = numpy.zeros(number_of_data_points, numpy.float)
    y = numpy.zeros(number_of_data_points, numpy.float)
    yerr = numpy.zeros(number_of_data_points, numpy.float)

    for i in range(number_of_data_points):
        x[i] = words[i][0]
        y[i] = words[i][1]
        yerr[i] = words[i][2]
#        print('i, x[i], y[i], yerr[i]: ', i, x[i], y[i], yerr[i])
#    print('x,y,yerr: ', x, y, yerr)

    slope, slope_error, intercept, intercept_error, r_value, reduced_chi_squared, ycalc, diff = line_fit(
        x, y, yerr, number_of_data_points, mode)
    print('slope, error, intercept, error, r_value, reduced_chi_squared: ', slope,
          slope_error, intercept, intercept_error, r_value, reduced_chi_squared)

    outfile = io.open(output_file_name, 'w')
    outfile.write('input file: ' + input_file_name + '\n')
    outfile.write('number of points fit: ' + str(number_of_data_points) + '\n')
    outfile.write('slope: ' + str(round(slope, 4)) +
                  ' +/- ' + str(round(slope_error, 4)) + '\n')
    outfile.write('intercept: ' + str(round(intercept, 4)) +
                  '+/-' + str(round(intercept_error, 4)) + '\n')
    outfile.write('reduced chi-squared: ' +
                  str(round(reduced_chi_squared, 4)) + '\n')
    outfile.write('   x\t\t   y\t\t   yerr\t\t   ycalc\t   y-ycalc\n')
    for i in range(number_of_data_points):
        outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\n' %
                      (x[i], y[i], yerr[i], ycalc[i], diff[i]))
    outfile.close()

#Test results from stats linear regression and chisquare here#####
    def func(x, m, c):
        f = m * x + c
        return f
# If we use absolute_sigma = True, the standard deviations are the same as the errors on the slope and intercept obtained above.
    popt, covariance = scipy.optimize.curve_fit(
        func, x, y, sigma=yerr, absolute_sigma=True)
#    popt, covariance = scipy.optimize.curve_fit(func, x, y)
    print('curve_fit results: ', popt, covariance)
    if numpy.all(covariance != numpy.inf):
        #   compute one standard deviation of the molecular weights
        standard_deviation_slope = numpy.sqrt(covariance[0][0])
        standard_deviation_intercept = numpy.sqrt(covariance[1][1])
        print('standard_deviation_slope, standard_deviation_intercept: ',
              standard_deviation_slope, standard_deviation_intercept)
        correlation_coefficient = covariance[0][1] / \
            (standard_deviation_slope*standard_deviation_intercept)
        print('correlation coefficient: ', correlation_coefficient)

#    plt.figure()
    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.errorbar(x, y, yerr=yerr, fmt='bo', markersize=2)
    ax0.plot(x, ycalc, 'r--')
    ax0.set_title('data')
#
    ax1.plot(x, diff, 'ko', markersize=5)
    ax1.set_title('residuals')
    ax1.set_ylim([-0.5, 0.5])
#    plt.subplot(2,1,1)
#    plt.errorbar(x,y,yerr, fmt='bo')
#    plt.plot(x,ycalc,'r--')

#    plt.subplot(2,1,2)
    #plt.plot(x, diff, 'bo', x, diff, 'k')
#    plt.plot(x,diff, 'bo')
    plt.show()
