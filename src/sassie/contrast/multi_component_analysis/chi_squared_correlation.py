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
#    CHI_SQUARED_CORRELATION calculates the reduced chi-squared for a
#    weighted fit as well as the standard deviation array and the
#    correlation matrix from the covariance matrix
#
#       2/2024       --      initial coding and testing:   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Chi Squared Correlation** contains the **Get Reduced Chi Squared** method that calculates the reduced chi-squared for a weighted fit as well as the **Correlation from Covariance** method that calculates the standard deviation array and the correlation matrix from the covariance matrix.
    
'''
import numpy


def get_reduced_chi_squared(y, yerr, yfit, number_of_data_points, number_of_unknowns):
    r'''
    **Get Reduced Chi Squared** calculates the reduced chi-squared value of a weighted fit to an equation.

    Note
    ----

    :math:`\chi^2 = \sum_{i=1}^{n} \frac{(y_{i} - yfit_{i})^2}{yerr_i^2}` where

    :math:`y_{i}` is the :math:`i^{th}` experimental y value

    :math:`yfit_{i}` is the :math:`i^{th}` calculated y value  

    :math:`yerr_{i}` is the error on the :math:`i^{th}` experimental y value

    :math:`n` is the number of data points

    The reduced :math:`\chi^2` is defined as 

    :math:`\chi^2_{\text{red}} = \frac{\chi^2}{\nu}` where

    :math:`\chi^2` is the chi-squared statistic

    :math:`\nu` is the number of degrees of freedom 



    Parameters
    ----------
    number_of_data_points:  int
        The number of equations that are being solved simultaneously
    number_of_unknowns:  int
        The number of parameters that are being determined from the weighted fit 
    y:  float array (dimension = number_of_data_points)
        The :math:`y` values
    yerr:  float array (dimension = number_of_data_points)
        The errors on the :math:`y` values
    yfit: float array (dimension = number_of_data_points)
        The :math:`y` values calculated using the parameters from the weighted fit


    Returns
    -------
    reduced_chi_squared:  float
        The reduced chi-squared for the fit to the equation, i.e, chi-squared divided by: the number of data points - the number of unknowns; if the number of data points = the number of unknowns, the reduced chi-squared is set to -1 (which can be tested for outside of this method)

    '''

    diff = y - yfit
#    print('diff: ', diff)
    chi_squared = numpy.sum((diff/yerr)**2)
#    print('chi_squared: ', chi_squared)
    degrees_of_freedom = number_of_data_points - number_of_unknowns
#    print('degrees_of_freedom: ', degrees_of_freedom)
    if (degrees_of_freedom > 0):
        reduced_chi_squared = chi_squared/degrees_of_freedom
    else:
        reduced_chi_squared = -1
#    print('reduced chi-squared: ', reduced_chi_squared)
    return reduced_chi_squared


def correlation_from_covariance(covariance):
    r'''
    **Correlation from Covariance** calculates the standard deviation array and the correlation matrix from the covariance matrix.  The correlation coefficient is calculated by dividing the covariance between two variables by the product of their individual standard deviations. The correlation coefficient is a standardized measure that ranges from -1 to 1, where -1 indicates a perfect negative linear relationship, 1 indicates a perfect positive linear relationship, and 0 indicates no linear relationship.  The standard deviation of a variable is calculated by taking the square root of its corresponding diagonal element of the covariance matrix.

    Note
    ----

    The standard deviation of each variable is obtained from the covariance matrix as

    :math:`\sigma_{i} = \sqrt{COV_{ii}}` where

    :math:`COV_{ii}` is the diagonal element of the covariance matrix that corresonds to the :math:`i^{th}` variable

    The correlation matrix can be written in terms of the covariance matrix and standard deviations as

    :math:`COR_{ij} = \frac{COV_{ij}}{\sigma_{i} \sigma_{j}}`


    Parameters
    ----------
    covariance:  2D float array (dimension = number of fit parameters * number of fit parameters)
        The covariance matrix from a fit to an equation

    Returns
    -------
    standard_deviation:  float array (dimension = number of fit parameters)
        The standard deviations of the fitted parameters
    correlation:  2D float array (dimension = number of fit parameters * number of fit parameters)
        The correlation matrix calculated from the covariance matrix, i.e., the normalized covariance matrix

    '''

    standard_deviation = numpy.sqrt(numpy.diag(covariance))
#    print('standard_deviation: ', standard_deviation)
# take the outer product of standard_deviation * standard_deviation
    outer_standard_deviation = numpy.outer(
        standard_deviation, standard_deviation)
    correlation = covariance / outer_standard_deviation
# set the elements in the correlation matrix to 0 if the corresponding elements in the covariance matrix are 0.  This avoids very small numbers in the correlation matrix
    correlation[covariance == 0] = 0
#    print('correlation: ', correlation)
    return standard_deviation, correlation
