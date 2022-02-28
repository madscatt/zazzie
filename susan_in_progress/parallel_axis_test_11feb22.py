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
import locale
import numpy



def parallel_axis(radius_of_gyration, radius_of_gyration_error, delta_rho, volume_fraction, fraction_d2o, number_of_contrast_points, number_of_components):
#Parallel axis analysis
#inputContrastData
#Rg [0], RgErr [1], Drho1 [2], Drho2 [3], fV1 [4], fV2 [5], (fV1*Drho1 + fV2*Drho2) [6]

# X is the Hessian; Xi is the inverse Hessian; Y is the LSQ vector; I is the result
# of the LSQ; chi2 is chi^2; paaVector contains products of contrasts and volumes, and is
# the coefficients in the parallel axis theorem; w is the weight

    print('fraction_d2o: ', fraction_d2o)
    print('radius_of_gyration: ', radius_of_gyration)
    print('radius_of_gyration_error: ', radius_of_gyration_error)
    print('volume_fraction: ', volume_fraction)
    print('delta_rho: ', delta_rho)

    vectShape = (3,1) #3 unknowns? rg1, rg2 and D?
    I = numpy.zeros(vectShape)
    Y = numpy.zeros(vectShape) #LSQ vector
    paaVector = numpy.zeros(vectShape) #contains products of contrasts and volumes - coefficients for parallel axis theorem
    shape = (3,3)    
    X = numpy.zeros(shape) #Hessian
    Xi = numpy.zeros(shape) #inverse Hessian

    chi2 = 0.0

    for k in range(0,number_of_contrast_points):
        print('k, delta_rho[k][0], delta_rho[k][1]: ', k, delta_rho[k][0], delta_rho[k][1])

        paaVector[0] = delta_rho[k][0]*volume_fraction[0]/(delta_rho[k][0]*volume_fraction[0] + delta_rho[k][1]*volume_fraction[1])
#                           Drho1*fV1/(fV1*Drho1 + fV2*Drho2)
        paaVector[1] = 1.0 - paaVector[0] # calculate the coefficients  (Where does fV2 come in?)
        paaVector[2] = paaVector[0]*paaVector[1]
        
        w = 1.0/(2.0*radius_of_gyration[k]*radius_of_gyration_error[k]) #convert 1/sig(Rg) to 1/sig(R^2)
#               1/2*Rg*Rgerr       
        for j in range(0,3):
            Y[j] = Y[j] + w*w*paaVector[j]*radius_of_gyration[k]*radius_of_gyration[k]
#                                               Rg*Rg            
            for l in range(0,3):
                X[j][l] = X[j][l] + w*w*paaVector[j]*paaVector[l]
            
    Xi = numpy.linalg.inv(numpy.matrix(X)) # solve LSQ problem
#note that numpy manual advises using regular arrays, numpy.array, instead of numpy.matrix.

    print('Xi: ', Xi)

    I = Xi*Y
    print('I: ', I)


    for k in range(0,number_of_contrast_points):
        print('k, delta_rho[k][0], delta_rho[k][1]: ', k, delta_rho[k][0], delta_rho[k][1])

        paaVector[0] = delta_rho[k][0]*volume_fraction[0]/(delta_rho[k][0]*volume_fraction[0] + delta_rho[k][1]*volume_fraction[1])
        paaVector[1] = 1.0 - paaVector[0] 
        paaVector[2] = paaVector[0]*paaVector[1]
        
        w = 1.0/(2.0*radius_of_gyration[k]*radius_of_gyration_error[k]) #convert 1/sig(Rg) to 1/sig(R^2)

        chi2 = chi2 + w*w*(radius_of_gyration[k]*radius_of_gyration[k] - paaVector[0]*I[0] - paaVector[1]*I[1] -     paaVector[2]*I[2])*(radius_of_gyration[k]*radius_of_gyration[k] - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])/(number_of_contrast_points-3)
        
    print('chi2: ', chi2.item())
    print('I[0]: ',I[0].item())
    print('Xi(00): ', Xi.item(0,0))
    chi_squared = chi2.item()
    rg1_squared = I[0].item()
    rg1_squared_error = numpy.sqrt(chi_squared*Xi.item(0,0))
    rg2_squared = I[1].item()
    rg2_squared_error = numpy.sqrt(chi_squared*Xi.item(1,1))
    cm_distance_squared = I[2].item()
    cm_distance_squared_error = numpy.sqrt(chi_squared*Xi.item(2,2))
    print(rg1_squared, rg1_squared_error)
    print(rg2_squared, rg2_squared_error)
    print(cm_distance_squared, cm_distance_squared_error)
    
    
    print("\nParallel Axis Theorem Analysis")
    print("chi^2: " + str(chi2.item()))
    print("R1^2: " + str(I[0].item()) + ", " + str(numpy.sqrt(chi2.item()*Xi.item((0,0)))))
    print("R2^2: " + str(I[1].item()) + ", " + str(numpy.sqrt(chi2.item()*Xi.item((1,1)))))
    print("D^2: " + str(I[2].item()) + ", " + str(numpy.sqrt(chi2.item()*Xi.item((2,2)))))
    print("R1: " + str(numpy.sqrt(I[0].item())) + ", " + str(0.5/numpy.sqrt(I[0].item())*numpy.sqrt(chi2.item()*Xi.item((0,0)))))
    print("R2: " + str(numpy.sqrt(I[1].item())) + ", " + str(0.5/numpy.sqrt(I[1].item())*numpy.sqrt(chi2.item()*Xi.item((1,1)))))
    print("D: " + str(numpy.sqrt(I[2].item())) + ", " + str(0.5/numpy.sqrt(I[2].item())*numpy.sqrt(chi2.item()*Xi.item((2,2)))))

if __name__ == "__main__":

    import matplotlib.pyplot as plt 
    path = ('./')
    input_file_name = os.path.join(path,'test_parallel_axis_data.txt')
    print (input_file_name)
    output_file_name = os.path.join(path,'parallel_axis_testdata.out')
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

    
    chi2 = 0.0
    
    chi2 = parallel_axis(radius_of_gyration, radius_of_gyration_error, delta_rho, volume_fraction, fraction_d2o, number_of_contrast_points, number_of_components)
    
