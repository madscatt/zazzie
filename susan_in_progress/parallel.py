# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 14:17:43 2018

@author: kathryn.sarachan
"""
import math
import numpy as np
from numpy.linalg import inv

# data import assumes a perfectly formatted MulCh-Contrast output file
# line 1: description; line 2: number of contrast points; lines 3+: Rg, RgErr, Drho1, Drho2, fV1
# there is NO exception handling for this right now
file_path = "rg.txt"

inputContrastData = np.loadtxt(file_path, skiprows = 2)
# Rg, RgErr, Drho1, Drho2, fV1

#inputContrastData = np.append(inputContrastData, [1.0, 2.0, 3.0, 3.0, 3.0])
numPoints = inputContrastData.shape[0]

# vectors to hold fV2 and Drho
fv2_arr = np.zeros(numPoints)
drho_arr = np.zeros(numPoints)

count = 0
for item in inputContrastData:
    # eliminate possible divide-by-zero errors in RgErr
    if item[1] == 0.0: item[1] = 0.01
    # add fV2 to vector
    fv2_arr[count] = 1.0 - item[4]
    # add Drho total to vector
    drho_arr[count] = item[4]*item[2] + fv2_arr[count]*item[3]
    count += 1
# add fV2 and Drho vectors to data matrix
inputContrastData = np.column_stack((inputContrastData, fv2_arr))
inputContrastData = np.column_stack((inputContrastData, drho_arr))

#Parallel axis analysis

# X is the Hessian; Xi is the inverse Hessian; Y is the LSQ vector; I is the result
# of the LSQ; chi2 is chi^2; paaVector contains products of contrasts and volumes, and is
# the coefficients in the parallel axis theorem; w is the weight

vectShape = (3,1)
I = np.zeros(vectShape)
Y = np.zeros(vectShape) #LSQ vector
paaVector = np.zeros(vectShape) #contains products of contrasts and volumes - coefficients for parallel axis theorem
shape = (3,3)    
X = np.zeros(shape) #Hessian
Xi = np.zeros(shape) #inverse Hessian

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
            
Xi = inv(np.matrix(X)) # solve LSQ problem
#note that numpy manual advises using regular arrays, numpy.array, instead of numpy.matrix.

I = Xi*Y


for k in range(0,len(inputContrastData)):

    paaVector[0] = inputContrastData[k][2]*inputContrastData[k][4]/inputContrastData[k][6]
    paaVector[1] = 1.0 - paaVector[0]
    paaVector[2] = paaVector[0]*paaVector[1]
    
    w = 1.0/(2.0*inputContrastData[k][0]*inputContrastData[k][1]) 

    chi2 = chi2 + w*w*(inputContrastData[k][0]*inputContrastData[k][0] - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])*(inputContrastData[k][0]*inputContrastData[k][0] - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])/(inputContrastData.shape[0]-3)
    
print("\nParallel Axis Theorem Analysis")
print("chi^2: " + str(chi2.item()))
print("R1^2: " + str(I[0].item()) + ", " + str(math.sqrt(chi2*Xi.item((0,0)))))
print("R2^2: " + str(I[1].item()) + ", " + str(math.sqrt(chi2*Xi.item((1,1)))))
print("D^2: " + str(I[2].item()) + ", " + str(math.sqrt(chi2*Xi.item((2,2)))))
print("R1: " + str(math.sqrt(I[0].item())) + ", " + str(0.5/math.sqrt(I[0].item())*math.sqrt(chi2*Xi.item((0,0)))))
print("R2: " + str(math.sqrt(I[1].item())) + ", " + str(0.5/math.sqrt(I[1].item())*math.sqrt(chi2*Xi.item((1,1)))))
print("D: " + str(math.sqrt(I[2].item())) + ", " + str(0.5/math.sqrt(I[2].item())*math.sqrt(chi2*Xi.item((2,2)))))

