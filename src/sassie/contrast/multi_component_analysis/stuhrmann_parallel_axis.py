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
#    STUHRMANN PARALLEL AXIS is the method that performs both a Sturhmann
#    and parallel axis theorem analysis to obtain the radii of gyration of
#    the individual components in a two component complex and the distance
#    between their center of mass.  The results from the two methods should
#    be similar and can be compared.#
#
#       STUHRMANN is the method that fits Rg**2 vs 1/mvars.delta_rho data to a
#       2nd order polynomial and then calculates the radii of gyration for
#       each of two components and the distance between their center of mass
#
#       PARALLEL AXIS is the method that calculates the radii of gyration for
#       each of two components and the distance between their center of mass
#       using the parallel axis theorem.
#
#       7/18/2021       --  Parallel Axis initial coding:   Kathryn Sarachan
#       2/10/2022       --  Revised for SASSIE 2.0      :   Susan Krueger
#       6/3/2021        --  Stuhrmann initial coding    :   Kathryn Sarachan
#       2/14/2022       --  Revised for SASSIE 2.0      :   Susan Krueger
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
import time
import numpy
import read_contrast_output_files
import polynomial_function_fit

def stuhrmann(other_self):

    '''
    put detailed description of module, inputs and outputs here
    '''

#mvars used:  fraction_d2o, number_of_contrast_points, number_of_components, partial_specific_volume, izero, concentration, molecular_weight, partial_specific_volume, delta_rho 
#mcavars used:  outfile, volume_fraction

#    print('in stuhrmann')

    log = other_self.log
    log.debug('in stuhrmann')
    pgui = other_self.run_utils.print_gui
    
    
    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables
    
    log.debug(vars(mvars))
    log.debug(vars(mcavars))


#Stuhrmann analysis

    pgui('Stuhrmann Analysis Method\n')
    pgui('Calculating alpha, beta and Rm')

    vector_shape = (mvars.number_of_contrast_points,1)
    delta_rho_inverse = numpy.zeros(vector_shape)
    rg_squared = numpy.zeros(vector_shape)
    rg_squared_error = numpy.zeros(vector_shape)
    rg_squared_calculated = numpy.zeros(vector_shape)
    diff = numpy.zeros(vector_shape)
    
#    print('volume fraction in stuhrmann: ', mcavars.volume_fraction)

    # Rg**2 = Rm**2 + alpha/mvars.delta_rho - beta/mvars.delta_rho
    # X is the x-coordinates (1/mvars.delta_rho), Y is the y-coordinates (Rg**2), and
    # yerr is the error in the y-coordinates (sig(Rg**2), propagated from sig(Rg))
    
    for i in range(mvars.number_of_contrast_points):
        delta_rho_inverse[i] = 1.0/(mvars.delta_rho[i][0]*mcavars.volume_fraction[0] + mvars.delta_rho[i][1]*mcavars.volume_fraction[1])
        rg_squared[i] = mvars.radius_of_gyration[i]*mvars.radius_of_gyration[i]
        rg_squared_error[i] = 2.0*mvars.radius_of_gyration[i]*mvars.radius_of_gyration_error[i]

    
    # the next statement fits a second-order polynomial to the data to get alpha, beta and Rm

    reduced_chi_squared,fit,correlation = polynomial_function_fit.polynomial_fit(2,delta_rho_inverse,rg_squared,rg_squared_error,mvars.number_of_contrast_points)
#    print('chi_squared, fit, correlation: ', chi_squared, fit, correlation)

    for i in range(mvars.number_of_contrast_points):
        rg_squared_calculated[i] = fit[2].item() + (fit[1].item()*delta_rho_inverse[i].item()) + (fit[0].item()*delta_rho_inverse[i].item()*delta_rho_inverse[i].item())
        diff[i] = rg_squared[i] - rg_squared_calculated[i]
#    print('diff: ', diff)        
        
#TODO:  the info below needs to go into the docstring
    #B is the RHS of the equations 5a-c given by Olah 1994; Bi is the inverse of B; C contains the solution (i.e. R1, R2 and D). 
    #NOTE: It looks like only the mvars.delta_rho values from the first contrast point are used, i.e. mvars.delta_rho[0]. The results should be the same no matter which contrast values are used.
    #TODO: test results using mvars.delta_rho for all of the contrasts to see how close they are. Should all of these results be reported? Should the values be averaged and the standard error of the mean reported? 
    
    pgui('Calculating individual Rg values and distance between centers of mass \n')

    B = numpy.zeros((3,3))

    B[2][0] = mcavars.volume_fraction[0]
    B[2][1] = mcavars.volume_fraction[1]
    B[2][2] = mcavars.volume_fraction[0]*mcavars.volume_fraction[1]
    B[1][0] = (mvars.delta_rho[0][0] - mvars.delta_rho[0][1])*mcavars.volume_fraction[0]*mcavars.volume_fraction[1]
    B[1][1] = -(mvars.delta_rho[0][0] - mvars.delta_rho[0][1])*mcavars.volume_fraction[0]*mcavars.volume_fraction[1]
    B[1][2] = (mvars.delta_rho[0][0] - mvars.delta_rho[0][1])*mcavars.volume_fraction[0]*mcavars.volume_fraction[1]*(mcavars.volume_fraction[1]*mcavars.volume_fraction[1] - mcavars.volume_fraction[0]*mcavars.volume_fraction[0])
    B[0][0] = 0.0
    B[0][1] = 0.0
    B[0][2] = -(mvars.delta_rho[0][0] - mvars.delta_rho[0][1])*(mvars.delta_rho[0][0] - mvars.delta_rho[0][1])*mcavars.volume_fraction[0]*mcavars.volume_fraction[0]*mcavars.volume_fraction[1]*mcavars.volume_fraction[1]

#    print('B matrix: ', B)
    Bi = numpy.linalg.inv(numpy.matrix(B))
#    print('inverse of B matrix: ', Bi)
#    print('fit: ', fit)
    C = Bi*fit
#    print('C matrix: ', C)

    # DR1R2D are the error estimates (accounting for parameter correlations of the derived parameters)

    DR1R2D = numpy.zeros((3,1))

    for i in range(0,3):
        for j in range (0,3):
            for k in range (0,3):
                DR1R2D[i] += reduced_chi_squared.item()*Bi.item((i,j))*correlation.item((k,j))*Bi.item((i,k))


    chi_squared_stuhrmann = reduced_chi_squared.item()
    beta = -fit[0].item()
    beta_error = numpy.sqrt(chi_squared_stuhrmann*correlation.item((0,0)))
    alpha = fit[1].item()
    alpha_error = numpy.sqrt(chi_squared_stuhrmann*correlation.item((1,1)))
    rg_infinite_contrast_squared = fit[2].item()
    rg_infinite_contrast_squared_error = numpy.sqrt(chi_squared_stuhrmann*correlation.item((2,2)))
    rg_infinite_contrast = numpy.sqrt(rg_infinite_contrast_squared)
    rg_infinite_contrast_error = 0.5/numpy.sqrt(fit[2].item())*numpy.sqrt(chi_squared_stuhrmann*correlation.item((2,2)))
    rg1_squared = C[0].item()
    rg1_squared_error = numpy.sqrt(DR1R2D[0].item())
    rg2_squared = C[1].item()
    rg2_squared_error = numpy.sqrt(DR1R2D[1].item())
    cm_distance_squared = C[2].item()
    cm_distance_squared_error = numpy.sqrt(DR1R2D[2].item())
    rg1_stuhrmann = numpy.sqrt(rg1_squared)
    rg1_error_stuhrmann = 0.5/numpy.sqrt(C[0].item())*numpy.sqrt(DR1R2D[0].item())
    rg2_stuhrmann = numpy.sqrt(rg2_squared)
    rg2_error_stuhrmann = 0.5/numpy.sqrt(C[1].item())*numpy.sqrt(DR1R2D[1].item())
    cm_distance_stuhrmann = numpy.sqrt(cm_distance_squared)
    cm_distance_error_stuhrmann = 0.5/numpy.sqrt(C[2].item())*numpy.sqrt(DR1R2D[2].item())

#    for i in range(0, mvars.number_of_contrast_points):
#        print('i, rho^-1, Rg^2, sigma Rg^2, Rg^2 calc, diff ', i,delta_rho_inverse[i], rg_squared[i], rg_squared_error[i], rg_squared_calculated[i], diff[i])
#    print(beta, beta_error)
#    print(alpha, alpha_error)
#    print(rg_infinite_contrast, rg_infinite_contrast_error)
#    print(chi_squared_stuhrmann)    
#    print(rg1_stuhrmann, rg1_error_stuhrmann)
#    print(rg2_stuhrmann, rg2_error_stuhrmann)
#    print(cm_distance_stuhrmann, cm_distance_error_stuhrmann)

    pgui('Results from Stuhrmann Analysis:')
    pgui('beta: ' + str(round(beta,4)) + ' +/- ' + str(round(beta_error,4)))
    pgui('alpha: ' + str(round(alpha,4)) + ' +/- ' + str(round(alpha_error,4)))
    if(alpha > 0):
        pgui('alpha is positive, indicating that the component with a higher scattering length density lies toward the periphery of the complex.')
    elif(alpha < 0):
        pgui('alpha is negative, indicating that the component with a higher scattering length density lies toward the interior of the complex.')
    pgui('Rm: ' + str(round(rg_infinite_contrast,4)) + ' +/- ' + str(round(rg_infinite_contrast_error,4)) + '\n')
    pgui('R1: ' + str(round(rg1_stuhrmann,4)) + ' +/- ' + str(round(rg1_error_stuhrmann,4)))
    pgui('R2: ' + str(round(rg2_stuhrmann,4)) + ' +/- ' + str(round(rg2_error_stuhrmann,4)))
    pgui('D: ' + str(round(cm_distance_stuhrmann,4)) + ' +/- ' + str(round(cm_distance_error_stuhrmann,4)))
    pgui('reduced chi-squared: ' + str(round(chi_squared_stuhrmann,4)) + '\n')


    mcavars.outfile.write('Results from Stuhrmann Analysis: \n')
    mcavars.outfile.write('beta: ' + str(round(beta,4)) + ' +/- ' + str(round(beta_error,4)) + '\n')
    mcavars.outfile.write('alpha: ' + str(round(alpha,4)) + ' +/- ' + str(round(alpha_error,4)) + '\n')
    if(alpha > 0):
        mcavars.outfile.write('alpha is positive, indicating that the component with a higher scattering length density lies toward the periphery of the complex.\n')
    elif(alpha < 0):
        mcavars.outfile.write('alpha is negative, indicating that the component with a higher scattering length density lies toward the interior of the complex.\n')
    mcavars.outfile.write('Rm: ' + str(round(rg_infinite_contrast,4)) + ' +/- ' + str(round(rg_infinite_contrast_error,4)) + '\n')
    mcavars.outfile.write('R1: ' + str(round(rg1_stuhrmann,4)) + ' +/- ' + str(round(rg1_error_stuhrmann,4)) + '\n')
    mcavars.outfile.write('R2: ' + str(round(rg2_stuhrmann,4)) + ' +/- ' + str(round(rg2_error_stuhrmann,4)) + '\n')
    mcavars.outfile.write('D: ' + str(round(cm_distance_stuhrmann,4)) + ' +/- ' + str(round(cm_distance_error_stuhrmann,4)) + '\n')
    mcavars.outfile.write('reduced chi-squared: ' + str(round(chi_squared_stuhrmann,4)) + '\n\n')
    mcavars.outfile.write('1/delta_rho\t Rg^2 exp\t Rg^2err\t Rg^2 calc\t diff\n')
    for i in range(mvars.number_of_contrast_points):
        mcavars.outfile.write('%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\n' % (delta_rho_inverse[i], rg_squared[i], rg_squared_error[i], rg_squared_calculated[i], diff[i]))

    mcavars.outfile.close()

    return


def parallel_axis(other_self):

    '''
    put detailed description of module, inputs and outputs here
    '''

#mvars used:  fraction_d2o, number_of_contrast_points, number_of_components, partial_specific_volume, izero, concentration, molecular_weight, partial_specific_volume, delta_rho 
#mcavars used:  outfile, volume_fraction

#    print('in parallel axis')

    log = other_self.log
    log.debug('in parallel axis')
    pgui = other_self.run_utils.print_gui
        
    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables
    
    log.debug(vars(mvars))
    log.debug(vars(mcavars))
    
#    print('volume fraction in parallel axis: ', mcavars.volume_fraction)

    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in xrange(60)])

    pgui('\n%s \n' % (st))
    pgui('DATA FROM RUN: %s \n\n' % (ttxt))

    pgui('results written to output file: %s' % (mcavars.multi_component_analysis_path+mvars.output_file_name))
    pgui('-------------------------------')
    pgui('\nNumber of points fit: ' + str(mvars.number_of_contrast_points) + '\n')

#Parallel axis analysis

# X is the Hessian; Xi is the inverse Hessian; Y is the LSQ vector; I is the result
# of the LSQ; chi_squared is chi^2; paaVector contains products of contrasts and volumes, and is
# the coefficients in the parallel axis theorem; w is the weight

    pgui('Parallel Axis Theorem Method\n')
    pgui('Calculating individual Rg values and distance between centers of mass \n')

    vector_shape = (3,1) #3 unknowns rg1, rg2 and D?
    I = numpy.zeros(vector_shape)
    Y = numpy.zeros(vector_shape) #LSQ vector
    paaVector = numpy.zeros(vector_shape) #contains products of contrasts and volumes - coefficients for parallel axis theorem
    shape = (3,3)    
    X = numpy.zeros(shape) #Hessian
    Xi = numpy.zeros(shape) #inverse Hessian

    reduced_chi_squared = 0.0
#    print('starting chi_squared: ', chi_squared)

    #define variables for parallel axis equation
    #Rg**2 = f1*Rg1**2 + f2*Rg2**2 + f1*f2*D**2, where f1 = drho1*vf1/(drho1*vf1 + drho2*vf2), f2 = 1 - f1

    for k in range(0,mvars.number_of_contrast_points):
#        print('k, mvars.delta_rho[k][0], mvars.delta_rho[k][1]: ', k, mvars.delta_rho[k][0], mvars.delta_rho[k][1])

        paaVector[0] = mvars.delta_rho[k][0]*mcavars.volume_fraction[0]/(mvars.delta_rho[k][0]*mcavars.volume_fraction[0] + mvars.delta_rho[k][1]*mcavars.volume_fraction[1]) #f1 = drho1*vf1/(vf1*drho1 + vf2*drho2)
        paaVector[1] = 1.0 - paaVector[0] #f2 = 1 - f1
        paaVector[2] = paaVector[0]*paaVector[1] #f1*f2        
        w = 1.0/(2.0*mvars.radius_of_gyration[k]*mvars.radius_of_gyration_error[k]) #1/sigma(Rg**2) = 1/(2*Rg*sigma(Rg))
#               1/2*Rg*Rgerr       
        for j in range(0,3):
            Y[j] = Y[j] + w*w*paaVector[j]*mvars.radius_of_gyration[k]*mvars.radius_of_gyration[k]
#                                               Rg*Rg            
            for l in range(0,3):
                X[j][l] = X[j][l] + w*w*paaVector[j]*paaVector[l]

#    print('X matrix: ', X)            
    Xi = numpy.linalg.inv(numpy.matrix(X)) # solve LSQ problem
#note that numpy manual advises using regular arrays, numpy.array, instead of numpy.matrix.

#    print('Xi: ', Xi)

    I = Xi*Y
#    print('I: ', I)


    for k in range(0,mvars.number_of_contrast_points):
#        print('k, mvars.delta_rho[k][0], mvars.delta_rho[k][1]: ', k, mvars.delta_rho[k][0], mvars.delta_rho[k][1])

        paaVector[0] = mvars.delta_rho[k][0]*mcavars.volume_fraction[0]/(mvars.delta_rho[k][0]*mcavars.volume_fraction[0] + mvars.delta_rho[k][1]*mcavars.volume_fraction[1])
        paaVector[1] = 1.0 - paaVector[0] 
        paaVector[2] = paaVector[0]*paaVector[1]
        
        w = 1.0/(2.0*mvars.radius_of_gyration[k]*mvars.radius_of_gyration_error[k]) #convert 1/sig(Rg) to 1/sig(R^2)

        reduced_chi_squared = reduced_chi_squared + w*w*(mvars.radius_of_gyration[k]*mvars.radius_of_gyration[k] - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])*(mvars.radius_of_gyration[k]*mvars.radius_of_gyration[k] - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])/(mvars.number_of_contrast_points-3)
        
    chi_squared_parallel = reduced_chi_squared.item()
    rg1_squared = I[0].item()
    rg1_squared_error = numpy.sqrt(chi_squared_parallel*Xi.item(0,0))
    rg2_squared = I[1].item()
    rg2_squared_error = numpy.sqrt(chi_squared_parallel*Xi.item(1,1))
    cm_distance_squared = I[2].item()
    cm_distance_squared_error = numpy.sqrt(chi_squared_parallel*Xi.item(2,2))
    rg1_parallel = numpy.sqrt(rg1_squared)
    rg1_error_parallel = 0.5/numpy.sqrt(I[0].item())*numpy.sqrt(chi_squared_parallel*Xi.item(0,0))
    rg2_parallel = numpy.sqrt(rg2_squared)
    rg2_error_parallel = 0.5/numpy.sqrt(I[1].item())*numpy.sqrt(chi_squared_parallel*Xi.item(1,1))
    cm_distance_parallel = numpy.sqrt(cm_distance_squared)
    cm_distance_error_parallel = 0.5/numpy.sqrt(I[2].item())*numpy.sqrt(chi_squared_parallel*Xi.item(2,2))

#    print(chi_squared_parallel)
#    print(rg1_parallel, rg1_error_parallel)
#    print(rg2_parallel, rg2_error_parallel)
#    print(cm_distance_parallel, cm_distance_error_parallel)  

    pgui('Results from Parallel Axis Theorem Analysis: ')
    pgui('R1: ' + str(round(rg1_parallel,4)) + ' +/- ' + str(round(rg1_error_parallel,4)))
    pgui('R2: ' + str(round(rg2_parallel,4)) + ' +/- ' + str(round(rg2_error_parallel,4)))
    pgui('D: ' + str(round(cm_distance_parallel,4)) + ' +/- ' + str(round(cm_distance_error_parallel,4)))
    pgui('reduced chi-squared: ' + str(round(chi_squared_parallel,4)) + '\n')
    
    mcavars.outfile.write('\nnumber of points fit: ' + str(mvars.number_of_contrast_points) + '\n\n')
    mcavars.outfile.write('Results from Parallel Axis Theorem Analysis: \n')
    mcavars.outfile.write('R1: ' + str(round(rg1_parallel,4)) + ' +/- ' + str(round(rg1_error_parallel,4)) + '\n')
    mcavars.outfile.write('R2: ' + str(round(rg2_parallel,4)) + ' +/- ' + str(round(rg2_error_parallel,4)) + '\n')
    mcavars.outfile.write('D: ' + str(round(cm_distance_parallel,4)) + ' +/- ' + str(round(cm_distance_error_parallel,4)) + '\n')
    mcavars.outfile.write('reduced chi-squared: ' + str(round(chi_squared_parallel,4)) + '\n\n')
    
    return 
