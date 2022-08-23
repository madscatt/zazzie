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
#       STUHRMANN is the method that fits Rg**2 vs 1/delta_rho data to a
#       2nd order polynomial and then calculates the radii of gyration for
#       each of two components and the distance between their center of mass
#
#       6/3/2021        --      initial coding         :   Kathryn Sarachan
#       2/14/2022       --      revised for SASSIE 2.0 :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''


    INPUTS:
      
        
    
    OUTPUTS:
       

    Requires numpy.linalg.inv
        
'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import io
import numpy

# data import assumes a perfectly formatted MulCh-Contrast output file
# line 1: description; line 2: number of contrast points; lines 3+: Rg, RgErr, Drho1, Drho2, fV1
# there is NO exception handling for this right now
file_path = "rg.txt"

inputContrastData = numpy.loadtxt(file_path, skiprows = 2)
# Rg, RgErr, Drho1, Drho2, fV1

numPoints = inputContrastData.shape[0]

# vectors to hold fV2 and Drho
fv2_arr = numpy.zeros(numPoints)
drho_arr = numpy.zeros(numPoints)

count = 0
for item in inputContrastData:
    # eliminate possible divide-by-zero errors in RgErr
    if item[1] == 0.0: item[1] = 0.01
    # add fV2 values to vector
    fv2_arr[count] = 1.0 - item[4]
    # add Drho total values to vector
    drho_arr[count] = item[4]*item[2] + fv2_arr[count]*item[3]
    count += 1
# add fV2 and Drho vectors to data matrix
inputContrastData = numpy.column_stack((inputContrastData, fv2_arr))
inputContrastData = numpy.column_stack((inputContrastData, drho_arr))
print('inputContrastData: ', inputContrastData)

#Stuhrmann analysis

# this function is the second-order polynomial fitting function used in MulCH
# it takes X, the x-coordinates (1/Drho), Y, the y-coordinates (Rg^2), and
# sigY, the error in the y-coordinates (sig(Rg^2))
# it returns fit, a vector with the cooeficients from the polynomial fit; corr, the correlation matrix, and chi2 
def polynomialFit(order, X, Y, sigY, n):
    if order > 2: print("Polynomial fit can only fit a 2nd order polynomial or less.")
    
    vectShape = (order+1,1)
    A = numpy.zeros(vectShape)
    M = numpy.zeros(vectShape)
    shape = (order+1,order+1)
    B = numpy.zeros(shape)
    Bi = numpy.zeros(shape)
    
    for i in range (0,n):
        if sigY[0] > 0.0:
            w = 1.0/(sigY[i]*sigY[i])
            for j in range (0,order+1):
                A[j] += w*Y[i]*X[i]**(order-j)
                for k in range(0,order+1):
                    B[j][k] += w*X[i]**(order-j)*X[i]**(order-k)
    
    Bi = numpy.linalg.inv(numpy.matrix(B))
    M = Bi*A
    
    chi2 = 0.0
    nDisregarded = 0
    
    if n > order+1:
        for i in range(0,n):
            if sigY[i] > 0.0:
                delta = 0.0
                for j in range(0,order+1):
                    delta += M[j]*X[i]**(order-j)
                chi2 += (delta - Y[i])*(delta - Y[i])/sigY[i]/sigY[i]
            else: nDisregarded = nDisregarded+1
            
    chi2 = chi2/float((n-nDisregarded)-(order+1))
    
    return [chi2, M, Bi]
        

vectShape = (numPoints,1)
X = numpy.zeros(vectShape)
Y = numpy.zeros(vectShape)
sigY = numpy.zeros(vectShape)

count = 0

# X is the x-coordinates (1/Drho), Y, the y-coordinates (Rg^2), and
# sigY, the error in the y-coordinates (sig(Rg^2))

for item in inputContrastData:
    X[count] = 1.0/item[6]
    Y[count] = item[0]*item[0]
    sigY[count] = 2.0*item[0]*item[1]
    count += 1

print('X, Y, sigY: ', X, Y, sigY)
    
# the next statement fits a second-order polynomial to the data

chi2,fit,corr = polynomialFit(2,X,Y,sigY,numPoints)

print('chi2, fit, corr: ', chi2, fit, corr)

print('inputContrastData: ', inputContrastData)
#B is the RHS of the equations 5a-c given by Olah 1994; Bi is the inverse of B; C contains the solution (i.e. R1, R2 and D)

B = numpy.zeros((3,3))
print('inputContrastData[0]: ',inputContrastData[0])
B[2][0] = inputContrastData[0][4]
B[2][1] = inputContrastData[0][5]
B[2][2] = inputContrastData[0][4]*inputContrastData[0][5]
B[1][0] = (inputContrastData[0][2] - inputContrastData[0][3])*inputContrastData[0][4]*inputContrastData[0][5]
B[1][1] = -(inputContrastData[0][2] - inputContrastData[0][3])*inputContrastData[0][4]*inputContrastData[0][5]
B[1][2] = (inputContrastData[0][2] - inputContrastData[0][3])*inputContrastData[0][4]*inputContrastData[0][5]*(inputContrastData[0][5]*inputContrastData[0][5] - inputContrastData[0][4]*inputContrastData[0][4])
B[0][0] = 0.0
B[0][1] = 0.0
B[0][2] = -(inputContrastData[0][2] - inputContrastData[0][3])*(inputContrastData[0][2] - inputContrastData[0][3])*inputContrastData[0][4]*inputContrastData[0][4]*inputContrastData[0][5]*inputContrastData[0][5]

print('B matrix: ', B)

Bi = numpy.linalg.inv(numpy.matrix(B))
C = Bi*fit


print("\nStuhrmann plot data")
print("Delta-rho^-1, Rg^2 (exp), sigma(Rg^2), Rg^2 (calc)")

for i in range(0, numPoints):
    print(str(X[i].item()) + ", " + str(Y[i].item()) + ", " + str(sigY[i].item()) + ", " + str(fit[2].item() + (fit[1].item()*X[i].item()) + (fit[0].item()*X[i].item()*X[i].item()))) 

# DR1R2D are the error estimates (accounting for parameter correlations of the derived parameters)

DR1R2D = numpy.zeros((3,1))

for i in range(0,3):
    for j in range (0,3):
        for k in range (0,3):
            DR1R2D[i] += chi2.item()*Bi.item((i,j))*corr.item((k,j))*Bi.item((i,k))

print("\nStuhrmann analysis")
print("Parameter: Value, ESD")

print("chi2: " + str(chi2.item()))

print("\nSturhmann plot coefficients")
print("Rm^2: " + str(fit[2].item()) + ", " + str(numpy.sqrt(chi2.item()*corr.item((2,2)))))
print("alpha(dagger): " + str(fit[1].item()) + ", " + str(numpy.sqrt(chi2.item()*corr.item((1,1)))))
print("beta: " + str(-fit[0].item()) + ", " + str(numpy.sqrt(chi2.item()*corr.item((0,0)))) + "\n")

print("Extracted parameters")
print("R1^2: " + str(C[0].item()) + ", " + str(numpy.sqrt(DR1R2D[0].item())))
print("R2^2: " + str(C[1].item()) + ", " + str(numpy.sqrt(DR1R2D[1].item())))
print("D^2: " + str(C[2].item()) + ", " + str(numpy.sqrt(DR1R2D[2].item())))
print("R1: " + str(numpy.sqrt(C[0].item())) + ", " + str(0.5/numpy.sqrt(C[0].item())*numpy.sqrt(DR1R2D[0].item())))
print("R2: " + str(numpy.sqrt(C[1].item())) + ", " + str(0.5/numpy.sqrt(C[1].item())*numpy.sqrt(DR1R2D[1].item())))
print("D: " + str(numpy.sqrt(C[2].item())) + ", " + str(0.5/numpy.sqrt(C[2].item())*numpy.sqrt(DR1R2D[2].item())))
print("Rm: " + str(numpy.sqrt(fit[2].item())) + ", " + str(0.5/numpy.sqrt(fit[2].item())*numpy.sqrt(chi2.item()*corr.item((2,2)))))
