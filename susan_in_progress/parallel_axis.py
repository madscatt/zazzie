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
#       PARALLEL AXIS is the method that calculates the radii of gyration for
#       each of two components and the distance between their center of mass
#
#       7/18/2021       --      initial coding         :   Kathryn Sarachan
#       2/10/2022       --      revised for SASSIE 2.0 :   Susan Krueger
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
import sys
import math

# data import assumes a perfectly formatted MulCh-Contrast output file
# line 1: description; line 2: number of contrast points; lines 3+: Rg, RgErr, Drho1, Drho2, fV1
# there is NO exception handling for this right now
file_path = "rg.txt"

inputContrastData = numpy.loadtxt(file_path, skiprows = 2)
# 0,  1,     2,     3,     4
# Rg, RgErr, Drho1, Drho2, fV1

#inputContrastData = numpy.append(inputContrastData, [1.0, 2.0, 3.0, 3.0, 3.0])
numPoints = inputContrastData.shape[0]

# vectors to hold fV2 and Drho
fv2_arr = numpy.zeros(numPoints)
drho_arr = numpy.zeros(numPoints)

count = 0
for item in inputContrastData:
    # eliminate possible divide-by-zero errors in RgErr
    if item[1] == 0.0: item[1] = 0.01
    # add fV2 to vector
    fv2_arr[count] = 1.0 - item[4]
    # add Drho total to vector
    drho_arr[count] = item[4]*item[2] + fv2_arr[count]*item[3]
# drho_arr = fV1*Drho1 + fV2*Drho2
    count += 1
# add fV2 and Drho vectors to data matrix
inputContrastData = numpy.column_stack((inputContrastData, fv2_arr))
inputContrastData = numpy.column_stack((inputContrastData, drho_arr))

print('inputContrastData: ', inputContrastData)
print('length of inputContrastData: ', len(inputContrastData))

#Parallel axis analysis

# X is the Hessian; Xi is the inverse Hessian; Y is the LSQ vector; I is the result
# of the LSQ; chi2 is chi^2; paaVector contains products of contrasts and volumes, and is
# the coefficients in the parallel axis theorem; w is the weight

vectShape = (3,1)
I = numpy.zeros(vectShape)
Y = numpy.zeros(vectShape) #LSQ vector
paaVector = numpy.zeros(vectShape) #contains products of contrasts and volumes - coefficients for parallel axis theorem
shape = (3,3)    
X = numpy.zeros(shape) #Hessian
Xi = numpy.zeros(shape) #inverse Hessian

chi2 = 0.0

for k in range(0,numPoints):

    paaVector[0] = inputContrastData[k][2]*inputContrastData[k][4]/inputContrastData[k][6]
    paaVector[1] = 1.0 - paaVector[0] # calculate the coefficients
    paaVector[2] = paaVector[0]*paaVector[1]
    
    w = 1.0/(2.0*inputContrastData[k][0]*inputContrastData[k][1]) #convert 1/sig(Rg) to 1/sig(R^2)
    
    for j in range(0,3):
        Y[j] = Y[j] + w*w*paaVector[j]*inputContrastData[k][0]*inputContrastData[k][0]
        
        for l in range(0,3):
            X[j][l] = X[j][l] + w*w*paaVector[j]*paaVector[l]
            
print('paaVector: ', paaVector)


Xi = numpy.linalg.inv(numpy.matrix(X)) # solve LSQ problem
#note that numpy manual advises using regular arrays, numpy.array, instead of numpy.matrix.

print('Xi: ', Xi)
I = Xi*Y

print('I: ', I)

for k in range(0,len(inputContrastData)):

    paaVector[0] = inputContrastData[k][2]*inputContrastData[k][4]/inputContrastData[k][6]
    paaVector[1] = 1.0 - paaVector[0]
    paaVector[2] = paaVector[0]*paaVector[1]
    
    w = 1.0/(2.0*inputContrastData[k][0]*inputContrastData[k][1]) 

    chi2 = chi2 + w*w*(inputContrastData[k][0]*inputContrastData[k][0] - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])*(inputContrastData[k][0]*inputContrastData[k][0] - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])/(inputContrastData.shape[0]-3)
    
print('paaVector: ', paaVector)
print('w: ', w)
print('chi2: ', chi2)
    
print("\nParallel Axis Theorem Analysis")
print("chi^2: " + str(chi2.item()))
print("R1^2: " + str(I[0].item()) + ", " + str(math.sqrt(chi2*Xi.item((0,0)))))
print("R2^2: " + str(I[1].item()) + ", " + str(math.sqrt(chi2*Xi.item((1,1)))))
print("D^2: " + str(I[2].item()) + ", " + str(math.sqrt(chi2*Xi.item((2,2)))))
print("R1: " + str(math.sqrt(I[0].item())) + ", " + str(0.5/math.sqrt(I[0].item())*math.sqrt(chi2*Xi.item((0,0)))))
print("R2: " + str(math.sqrt(I[1].item())) + ", " + str(0.5/math.sqrt(I[1].item())*math.sqrt(chi2*Xi.item((1,1)))))
print("D: " + str(math.sqrt(I[2].item())) + ", " + str(0.5/math.sqrt(I[2].item())*math.sqrt(chi2*Xi.item((2,2)))))

