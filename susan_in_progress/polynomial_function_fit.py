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
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **

'''
    POLYNOMIAL FIT performs a weighted polynomial fit to a polynomial function. This
    is the version that was first written by Andrew Whitten (2/2006) for a 2nd order
    polynomial as part of the MulCh program. Rewritten in Python by Kathryn Sarachan (6/2021). 
    REFERENCE: Whitten, A.E., Cai, S, Trewhella, J. (2008). "MULCh: modules for the
    analysis of small-angle neutron contrast variation data from biomolecular assemblies",
    J. Appl. Cryst. 41, 222 - 226.

    INPUTS:
        order, order of polynomial
        x array
        y array
        yerr, error in y array
        number of data points
        
    
    OUTPUTS:
        M, a vector with the coefficients from the polynomial fit with weighting 1.0/yerr**2
        Bi, the correlation matrix
        reduced chi-squared, number of accepted data points for which yerr > 0 - (order + 1) degrees of freedom
        

    Can be called by any method performing a fit to a polynomial.
    
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy

def polynomial_fit(order, X, Y, yerr, n):
    vector_shape = (order+1,1)
    A = numpy.zeros(vector_shape)
    M = numpy.zeros(vector_shape)
    shape = (order+1,order+1)
    B = numpy.zeros(shape)
    Bi = numpy.zeros(shape)
    
    for i in range (0,n):
        if yerr[0] > 0.0:
            w = 1.0/(yerr[i]*yerr[i])
            for j in range (0,order+1):
                A[j] += w*Y[i]*X[i]**(order-j)
                for k in range(0,order+1):
                    B[j][k] += w*X[i]**(order-j)*X[i]**(order-k)
    
    Bi = numpy.linalg.inv(numpy.matrix(B))  #correlation matrix
    M = Bi*A    #fit coefficients
    
    chi_squared = 0.0
    nDisregarded = 0
    
    if n > order+1:
        for i in range(0,n):
            if yerr[i] > 0.0:
                delta = 0.0
                for j in range(0,order+1):
                    delta += M[j]*X[i]**(order-j)
                chi_squared += (delta - Y[i])*(delta - Y[i])/yerr[i]/yerr[i]
            else:
                nDisregarded = nDisregarded+1
            
    reduced_chi_squared = chi_squared/float((n-nDisregarded)-(order+1))
    
    return reduced_chi_squared, M, Bi
        
