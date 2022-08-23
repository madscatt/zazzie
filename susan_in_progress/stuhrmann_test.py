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

def stuhrmann(radius_of_gyration, radius_of_gyration_error, delta_rho, volume_fraction, fraction_d2o, number_of_contrast_points):

#Stuhrmann analysis
#TODO: This should be a separate function outside of stuhrmann (perhaps even a separate file like linear fit)
#TODO: Look into replacing this with scipy routine?
# this function is the second-order polynomial fitting function used in MulCH
# it takes X, the x-coordinates (1/Drho), Y, the y-coordinates (Rg^2), and
# sigY, the error in the y-coordinates (sig(Rg^2))
# it returns fit, a vector with the cooeficients from the polynomial fit; corr, the correlation matrix, and chi2 
    def polynomialFit(order, X, Y, sigY, n):
        #the routine seems general; is this just because the stuhrmann analysis equation is of 2nd order?
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

    #B is the RHS of the equations 5a-c given by Olah 1994; Bi is the inverse of B; C contains the solution (i.e. R1, R2 and D). 
    #NOTE: It looks like only the delta_rho values from the first contrast point are used. The results should be the same no matter which contrast values are used.
    #TODO: test results using delta_rho for all of the contrasts to see how close they are. Should all of these results be reported? Should the values be averaged and the standard error of the mean reported? 
    
    B = numpy.zeros((3,3))

    B[2][0] = volume_fraction[0]
    B[2][1] = volume_fraction[1]
    B[2][2] = volume_fraction[0]*volume_fraction[1]
    B[1][0] = (delta_rho[0][0] - delta_rho[0][1])*volume_fraction[0]*volume_fraction[1]
    B[1][1] = -(delta_rho[0][0] - delta_rho[0][1])*volume_fraction[0]*volume_fraction[1]
    B[1][2] = (delta_rho[0][0] - delta_rho[0][1])*volume_fraction[0]*volume_fraction[1]*(volume_fraction[1]*volume_fraction[1] - volume_fraction[0]*volume_fraction[0])
    B[0][0] = 0.0
    B[0][1] = 0.0
    B[0][2] = -(delta_rho[0][0] - delta_rho[0][1])*(delta_rho[0][0] - delta_rho[0][1])*volume_fraction[0]*volume_fraction[0]*volume_fraction[1]*volume_fraction[1]

    print('B matrix: ', B)
    Bi = numpy.linalg.inv(numpy.matrix(B))
    print('inverse of B matrix: ', Bi)
    print('fit: ', fit)
    C = Bi*fit
    print('C matrix: ', C)
    

    #If we decide to report the values for all contrasts, then Rg^2 (calc) would be reported for all contrasts
    print("\nStuhrmann plot data")
    print("Delta-rho^-1, Rg^2 (exp), sigma(Rg^2), Rg^2 (calc)")

    for i in range(0, number_of_contrast_points):
        print('i: ', i)
        print(str(delta_rho_inverse[i].item()) + ", " + str(rg_squared[i].item()) + ", " + str(rg_squared_error[i].item()) + ", " + str(fit[2].item() +     (fit[1].item()*delta_rho_inverse[i].item()) + (fit[0].item()*delta_rho_inverse[i].item()*delta_rho_inverse[i].item()))) 


    # DR1R2D are the error estimates (accounting for parameter correlations of the derived parameters)

    DR1R2D = numpy.zeros((3,1))

    for i in range(0,3):
        for j in range (0,3):
            for k in range (0,3):
                DR1R2D[i] += chi2.item()*Bi.item((i,j))*corr.item((k,j))*Bi.item((i,k))


    chi_squared_stuhrmann = chi2.item()
    beta = -fit[0].item()
    beta_error = numpy.sqrt(chi_squared_stuhrmann*corr.item((0,0)))
    alpha = fit[1].item()
    alpha_error = numpy.sqrt(chi_squared_stuhrmann*corr.item((1,1)))
    rg_infinite_contrast_squared = fit[2].item()
    rg_infinite_contrast_squared_error = numpy.sqrt(chi_squared_stuhrmann*corr.item((2,2)))
    rg_infinite_contrast = numpy.sqrt(rg_infinite_contrast_squared)
    rg_infinite_contrast_error = 0.5/numpy.sqrt(fit[2].item())*numpy.sqrt(chi_squared_stuhrmann*corr.item((2,2)))
    print(beta, beta_error)
    print(alpha, alpha_error)
    print(rg_infinite_contrast_squared, rg_infinite_contrast_squared_error)
    print(rg_infinite_contrast, rg_infinite_contrast_error)
    print(chi_squared_stuhrmann)
    rg1_squared = C[0].item()
    rg1_squared_error = numpy.sqrt(DR1R2D[0].item())
    rg2_squared = C[1].item()
    rg2_squared_error = numpy.sqrt(DR1R2D[1].item())
    cm_distance_squared = C[2].item()
    cm_distance_squared_error = numpy.sqrt(DR1R2D[2].item())
    print(rg1_squared, rg1_squared_error)
    print(rg2_squared, rg2_squared_error)
    print(cm_distance_squared, cm_distance_squared_error)
    rg1_stuhrmann = numpy.sqrt(rg1_squared)
    rg1_error_stuhrmann = 0.5/numpy.sqrt(C[0].item())*numpy.sqrt(DR1R2D[0].item())
    rg2_stuhrmann = numpy.sqrt(rg2_squared)
    rg2_error_stuhrmann = 0.5/numpy.sqrt(C[1].item())*numpy.sqrt(DR1R2D[1].item())
    cm_distance_stuhrmann = numpy.sqrt(cm_distance_squared)
    cm_distance_error_stuhrmann = 0.5/numpy.sqrt(C[2].item())*numpy.sqrt(DR1R2D[2].item())
    print(rg1_stuhrmann, rg1_error_stuhrmann)
    print(rg2_stuhrmann, rg2_error_stuhrmann)
    print(cm_distance_stuhrmann, cm_distance_error_stuhrmann)

    return alpha, alpha_error, beta, beta_error, rg_infinite_contrast, rg_infinite_contrast_error, chi_squared_stuhrmann, rg1_stuhrmann,rg1_error_stuhrmann, rg2_stuhrmann, rg2_error_stuhrmann, cm_distance_stuhrmann, cm_distance_error_stuhrmann 

if __name__ == "__main__":

    import matplotlib.pyplot as plt 
    path = ('./')
    input_file_name = os.path.join(path,'test_stuhrmann_data.txt')
    print (input_file_name)
    output_file_name = os.path.join(path,'stuhrmann_testdata.out')
    read_from_file = False
    number_of_components = 2    #hardwired at this time
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
#    for k in range(0,number_of_contrast_points):
#        print('k, delta_rho[k][0], delta_rho[k][1]: ', k, delta_rho[k][0], delta_rho[k][1])

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


#Note that RgErr = 0 check will be done at the input filter checking level in SASSIE.    
    for i in range(number_of_contrast_points):
        fraction_d2o[i] = words[i][0]
        radius_of_gyration[i] = words[i][1]
        radius_of_gyration_error[i] = words[i][2]

    print('fraction_d2o: ', fraction_d2o)
    print('radius_of_gyration: ', radius_of_gyration)
    print('radius_of_gyration_error: ', radius_of_gyration_error)


    alpha, alpha_error, beta, beta_error, rg_infinite_contrast, rg_infinite_contrast_error, chi_squared_stuhrmann, rg1_stuhrmann,rg1_error_stuhrmann, rg2_stuhrmann, rg2_error_stuhrmann, cm_distance_stuhrmann, cm_distance_error_stuhrmann = stuhrmann(radius_of_gyration, radius_of_gyration_error, delta_rho, volume_fraction, fraction_d2o, number_of_contrast_points)

    print('after call to stuhrmann')    
    print(beta, beta_error)
    print(alpha, alpha_error)
    print(rg_infinite_contrast, rg_infinite_contrast_error)
    print(chi_squared_stuhrmann)    
    print(rg1_stuhrmann, rg1_error_stuhrmann)
    print(rg2_stuhrmann, rg2_error_stuhrmann)
    print(cm_distance_stuhrmann, cm_distance_error_stuhrmann)