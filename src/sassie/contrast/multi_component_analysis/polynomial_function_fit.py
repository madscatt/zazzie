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
#       POLYNOMIAL FIT
#
#       2/2006       --      initial coding         :   Andrew Whitten
#       6/2021       --      converted to Python    :   Kathryn Sarachan
#       3/2022       --      revised for SASSIE 2.0 :   Susan Krueger
#       3/2023       --      revised for SASSIE 3.0 :   Joseph E. Curtis
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **

'''

    **Polynomial Function Fit** contains the method that performs a weighted fit to a polynomial function. This is the version that was first written by Andrew Whitten (2/2006) for a 2nd order polynomial as part of the MulCh program. Rewritten in Python by Kathryn Sarachan (6/2021). 
    
    **Reference:** Whitten, A.E., Cai, S, Trewhella, J. (2008). "MULCh: modules for the analysis of small-angle neutron contrast variation data from biomolecular assemblies", *J. Appl. Cryst.* **41**, 222 - 226.

    **Inputs:**
    
        order of polynomial, number of data points, x array, y array, error in y array

    **Outputs:**
    
        coefficients from the weighted polynomial fit, correlation matrix, reduced chi-squared

    Can be called by any method performing a fit to a polynomial.

'''

import numpy


def polynomial_fit(order, X, Y, yerr, n):
    '''Performs a weighted fit (1.0/yerr\ :sup:`2` ) to a polynomial function of order n.

    Parameters
    ----------

    order: int
        order of the polynomial
    n:  int
        number of data points    
    X:  float array (dimension = n)
        x data
    Y:  float array (dimension = n)
        y data
    yerr: float array (dimension = n)
        error in y data

    Returns
    -------

    M: 2D float array (dimensions = order+1 x 1)
        coefficients from the polynomial fit with weighting 1.0/yerr\ :sup:`2` 
    Bi: 2D float array (dimensions = order+1 x order+1)
        correlation matrix
    reduced chi_squared: float
        chi-squared divided by the number of accepted data points for which yerr > 0 - (order + 1) degrees of freedom

    '''

    vector_shape = (order+1, 1)
    A = numpy.zeros(vector_shape)
    M = numpy.zeros(vector_shape)
    shape = (order+1, order+1)
    B = numpy.zeros(shape)
    Bi = numpy.zeros(shape)

    for i in range(0, n):
        if yerr[0] > 0.0:
            w = 1.0/(yerr[i]*yerr[i])
            for j in range(0, order+1):
                A[j] += w*Y[i]*X[i]**(order-j)
                for k in range(0, order+1):
                    B[j][k] += w*X[i]**(order-j)*X[i]**(order-k)

    Bi = numpy.linalg.inv(numpy.matrix(B))  # correlation matrix
    M = Bi*A    # fit coefficients

    chi_squared = 0.0
    nDisregarded = 0

    if n > order+1:
        for i in range(0, n):
            if yerr[i] > 0.0:
                delta = 0.0
                for j in range(0, order+1):
                    delta += M[j]*X[i]**(order-j)
                chi_squared += (delta - Y[i])*(delta - Y[i])/yerr[i]/yerr[i]
            else:
                nDisregarded = nDisregarded+1

    reduced_chi_squared = chi_squared/float((n-nDisregarded)-(order+1))

    return reduced_chi_squared, M, Bi
