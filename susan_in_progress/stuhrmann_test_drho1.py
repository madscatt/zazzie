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

def stuhrmann(radius_of_gyration, radius_of_gyration_error, delta_rho, volume_fraction, fraction_d2o, number_of_contrast_points, number_of_components):



# data import assumes a perfectly formatted MulCh-Contrast output file
# line 1: description; line 2: number of contrast points; lines 3+: Rg, RgErr, Drho1, Drho2, fV1
# there is NO exception handling for this right now
#file_path = "rg.txt"

#inputContrastData = numpy.loadtxt(file_path, skiprows = 2)
# 0=Rg, 1=RgErr, 2=Drho1, 3=Drho2, 4=fV1

##number_of_contrast_points = inputContrastData.shape[0]

# vectors to hold fV2 and Drho
#fv2_arr = numpy.zeros(number_of_contrast_points)
#drho_arr = numpy.zeros(number_of_contrast_points)

#count = 0
#for item in inputContrastData:
    # eliminate possible divide-by-zero errors in RgErr
#    if item[1] == 0.0: item[1] = 0.01
    # add fV2 values to vector
##    fv2_arr[count] = 1.0 - item[4]
    # add Drho total values to vector
#    drho_arr[count] = item[4]*item[2] + fv2_arr[count]*item[3]
#    count += 1
# add fV2 and Drho vectors to data matrix
#inputContrastData = numpy.column_stack((inputContrastData, fv2_arr))
#inputContrastData = numpy.column_stack((inputContrastData, drho_arr))

#Stuhrmann analysis

# this function is the second-order polynomial fitting function used in MulCH
# it takes X, the x-coordinates (1/Drho), Y, the y-coordinates (Rg^2), and
# sigY, the error in the y-coordinates (sig(Rg^2))
# it returns fit, a vector with the cooeficients from the polynomial fit; corr, the correlation matrix, and chi2 
    def polynomialFit(order, X, Y, sigY, n):
        if order > 2: print("Polynomial fit can only fit a 2nd order polynomial or less.")
    
        vector_shape = (order+1,1)
        A = numpy.zeros(vector_shape)
        M = numpy.zeros(vector_shape)
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
    
        return chi2, M, Bi
        

    vector_shape = (number_of_contrast_points,1)
    delta_rho_inverse = numpy.zeros(vector_shape)
    rg_squared = numpy.zeros(vector_shape)
    rg_squared_error = numpy.zeros(vector_shape)


    # Rg**2 = Rm**2 + alpha/delta_rho - beta/delta_rho
    # X is the x-coordinates (1/delta_rho), Y is the y-coordinates (Rg**2), and
    # sigY is the error in the y-coordinates (sig(Rg**2), propagated from sig(Rg))
    
    for i in range(number_of_contrast_points):
        delta_rho_inverse[i] = 1.0/(delta_rho[i][0]*volume_fraction[0] + delta_rho[i][1]*volume_fraction[1])
        rg_squared[i] = radius_of_gyration[i]*radius_of_gyration[i]
        rg_squared_error[i] = 2.0*radius_of_gyration[i]*radius_of_gyration_error[i]

    
    # the next statement fits a second-order polynomial to the data to get alpha, beta and Rm

    chi2,fit,corr = polynomialFit(2,delta_rho_inverse,rg_squared,rg_squared_error,number_of_contrast_points)
    print('chi2, fit, corr: ', chi2, fit, corr)

    #B is the RHS of the equations 5a-c given by Olah 1994; Bi is the inverse of B; C contains the solution (i.e. R1, R2 and D). It looks like only the delta_rho values from the first contrast point are used. 

    B = numpy.zeros((3,3))

    B[2][0] = volume_fraction[0]
    B[2][1] = volume_fraction[1]
    B[2][2] = volume_fraction[0]*volume_fraction[1]
    B[1][0] = (delta_rho[1][0] - delta_rho[1][1])*volume_fraction[0]*volume_fraction[1]
    B[1][1] = -(delta_rho[1][0] - delta_rho[1][1])*volume_fraction[0]*volume_fraction[1]
    B[1][2] = (delta_rho[1][0] - delta_rho[1][1])*volume_fraction[0]*volume_fraction[1]*(volume_fraction[1]*volume_fraction[1] - volume_fraction[0]*volume_fraction[0])
    B[0][0] = 0.0
    B[0][1] = 0.0
    B[0][2] = -(delta_rho[1][0] - delta_rho[1][1])*(delta_rho[1][0] - delta_rho[1][1])*volume_fraction[0]*volume_fraction[0]*volume_fraction[1]*volume_fraction[1]

    print('B matrix: ', B)
    Bi = numpy.linalg.inv(numpy.matrix(B))
    print('inverse of B matrix: ', Bi)
    print('fit: ', fit)
    C = Bi*fit
    print('C matrix: ', C)
    


    print("\nStuhrmann plot data")
    print("Delta-rho^-1, Rg^2 (exp), sigma(Rg^2), Rg^2 (calc)")

    for i in range(0, number_of_contrast_points):
        print(str(delta_rho_inverse[i].item()) + ", " + str(rg_squared[i].item()) + ", " + str(rg_squared_error[i].item()) + ", " + str(fit[2].item() +     (fit[1].item()*delta_rho_inverse[i].item()) + (fit[0].item()*delta_rho_inverse[i].item()*delta_rho_inverse[i].item()))) 

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

    return chi2

if __name__ == "__main__":

    import matplotlib.pyplot as plt 
    path = ('./')
    input_file_name = os.path.join(path,'test_stuhrmann_data.txt')
    print (input_file_name)
    output_file_name = os.path.join(path,'stuhrmann_testdata.out')
    read_from_file = False
    number_of_components = 2
    molecular_weight = [50.7,11.7] #kDa
    partial_specific_volume = [0.73,0.73]
    delta_rho = [[2.31,6.06],[1.75,5.49],[1.18,4.93],[0.06,3.80],[-2.21,1.54],[-2.78,0.98],[-3.34,0.41]]
#contrast calculator numbers for delta_rho:
# 0.00	  2.313	  6.055
# 0.10	  1.749	  5.492
# 0.20	  1.184	  4.929
# 0.40	  0.055	  3.801
# 0.80	 -2.209	  1.542
# 0.90	 -2.776	  0.976
# 1.00	 -3.343	  0.410
    print('delta_rho: ', delta_rho)


#volume fraction:  volume = molecular_weight*partial_specific_volume/Na; volume fraction = v1/(v1+v2); Na cancels out and units of volume can be kept in kDa, as we just want the volume fraction
    total_volume = 0.0
    volume_fraction = numpy.zeros(number_of_components, numpy.float)
    for i in range(number_of_components):
        total_volume = total_volume + molecular_weight[i]*partial_specific_volume[i]

    for i in range(number_of_components):
        volume_fraction[i] = molecular_weight[i]*partial_specific_volume[i]/total_volume
    
    print('volume fraction: ', volume_fraction)

    with io.open(input_file_name, 'r') as infile:
        line = infile.readlines()
#    print(line)
    number_of_contrast_points = len(line)
    print ("number of data points: ", number_of_contrast_points)
    for k in range(0,number_of_contrast_points):
        print('k, delta_rho[k][0], delta_rho[k][1]: ', k, delta_rho[k][0], delta_rho[k][1])

# fraction_d2o, Rg, RgErr    
# split the lines into words
    words = []
    for i in range(number_of_contrast_points):
        value = line[i].split()
#        print(value)
        words.append(value)
#    print(words)

    fraction_d2o = numpy.zeros(number_of_contrast_points, numpy.float)
    radius_of_gyration = numpy.zeros(number_of_contrast_points, numpy.float)
    radius_of_gyration_error = numpy.zeros(number_of_contrast_points, numpy.float)
    concentration = numpy.zeros(number_of_contrast_points, numpy.float)
    concentration_error = numpy.zeros(number_of_contrast_points, numpy.float)

    
    for i in range(number_of_contrast_points):
        fraction_d2o[i] = words[i][0]
        radius_of_gyration[i] = words[i][1]
        radius_of_gyration_error[i] = words[i][2]
        if(i == 3):
            concentration[i] = 26.9
            concentration_error[i] = 1.3
        else:
            concentration[i] = 11.9
            concentration_error[i] = 0.6

    print('fraction_d2o: ', fraction_d2o)
    print('radius_of_gyration: ', radius_of_gyration)
    print('radius_of_gyration_error: ', radius_of_gyration_error)
    print('concentration: ', concentration)
    print('concentration_error: ', concentration_error)

    chi2 = stuhrmann(radius_of_gyration, radius_of_gyration_error, delta_rho, volume_fraction, fraction_d2o, number_of_contrast_points, number_of_components)
    
    
