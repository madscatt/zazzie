# -*- coding: utf-8 -*-
'''
    Guinier analysis performs a Guinier fit to the data in the range qRgmin to qRgmax where qRgmin is defined by the first data point and qRgmax is based on the qRg limit, both as specified by the user.
    
'''
import io
import os
import math
import sys
import polynomial_function_fit_test as polynomial_function_fit
import numpy


def guinier_fit(number_of_contrast_points, number_of_data_points, scattering_data, q_rg_limit_guinier, starting_data_point_guinier, initial_points_to_use_guinier, refine_scale_factor_flag, scale_factor, delta_rho_v):

    fit = numpy.zeros((2, 1))
    correlation = numpy.zeros((2, 2))
    points_used_guinier = numpy.zeros(number_of_contrast_points)
    chi_squared_guinier = numpy.zeros(number_of_contrast_points)
    rg_guinier = numpy.zeros(number_of_contrast_points)
    rg_error_guinier = numpy.zeros(number_of_contrast_points)
    izero_guinier = numpy.zeros(number_of_contrast_points)
    izero_error_guinier = numpy.zeros(number_of_contrast_points)
    q_rg_min_guinier = numpy.zeros(number_of_contrast_points)
    q_rg_max_guinier = numpy.zeros(number_of_contrast_points)
    q_min_guinier = numpy.zeros(number_of_contrast_points)
    q_max_guinier = numpy.zeros(number_of_contrast_points)

#    print('qRg limit: ', q_rg_limit_guinier)
#    print('scale factor passed to Guinier: ', scale_factor)

    for j in range(number_of_contrast_points):

        # if the starting point is not the first point, adjust the q**2, ln(I) and ln(I) error arrays accordingly
        start = starting_data_point_guinier[j] - 1
        number_of_points = number_of_data_points - start
        q_squared = numpy.zeros(number_of_points)
        ln_i = numpy.zeros(number_of_points)
        ln_i_error = numpy.zeros(number_of_points)

        for i in range(number_of_points):
            index = start+i
#            print('i, index, q: ', i, index, scattering_data[j][index][0])
            q_squared[i] = scattering_data[j][index][0] * \
                scattering_data[j][index][0]
#            print('i, q: ', i, math.sqrt(q_squared[i]))
            # in C, ln(something negative) throws a floating point error that just goes into the array. Here I needed a place holder to serve the same purpose (for now): 999.999
# NOTE: numpy.log returns a RuntimeWarning and not a ValueError as math.log does and "except RuntimeWarning" doesn't set the new value as desired.
            try:
                ln_i[i] = math.log(scattering_data[j][index][1])
            except ValueError as e:
                ln_i[i] = 999.999
            ln_i_error[i] = scattering_data[j][index][2] / \
                scattering_data[j][index][1]

        status = 0
# points to use is somewhat arbitrary.  It is input by the user since it depends on delta_q as to how many points will be in the Guinier region.
        points_to_use = initial_points_to_use_guinier[j]
#        print('initial points to use: ', points_to_use)

        while (status == 0):
            reduced_chi_squared, fit, correlation = polynomial_function_fit.polynomial_fit(
                1, q_squared, ln_i, ln_i_error, points_to_use)
            reduced_chi_squared = reduced_chi_squared[0, 0]
            correlation00 = correlation[0, 0]
            correlation11 = correlation[1, 1]
            slope = fit[0, 0]
            intercept = fit[1, 0]
#            print('j, slope, intercept: ', j, slope, intercept)

# another spot where C exception handling is dopey - negative sqrts just populate the answer with a code, rather than throwing an exception. Here, we're going to call any negative sqrt 0.00, which should always be safely  < q_rg_limit_guinier
# NOTE: numpy.sqrt returns a RuntimeWarning and not a ValueError as math.sqrt does and "except RuntimeWarning" doesn't set the new value as desired.
#            print('qmax: ', scattering_data[j][start + points_to_use - 1][0])
            try:
                limit_test = math.sqrt(-3.0*slope) * \
                    scattering_data[j][start + points_to_use - 1][0]
            except ValueError as e:
                limit_test = 0.00
#            print('j, limit_test: ', j, limit_test)

# WHY?  slope should be negative, otherwise not in a valid Guinier region?
#            if (limit_test < q_rg_limit_guinier) or slope > 0.0:
            if (limit_test < q_rg_limit_guinier):
                points_used_guinier[j] = points_to_use
                chi_squared_guinier[j] = reduced_chi_squared
#                print('j, points used, chi squared: ', j, points_used_guinier[j], chi_squared_guinier[j])

                # Rg, where slope = -Rg**2/3; Rg = sqrt(-3.0*slope)
                try:
                    rg_guinier[j] = math.sqrt(-3.0*slope)
                except ValueError as e:
                    rg_guinier[j] = 0.00
                # RgErr
                try:
                    rg_error_guinier[j] = 3.0*math.sqrt(
                        reduced_chi_squared*correlation00)/2.0/math.sqrt(-3.0*slope)
                except ValueError as e:
                    rg_error_guinier[j] = 0.01
#                print('j, Rg, error: ', j, rg_guinier[j], rg_error_guinier[j])

                # I0, where intercept = ln[I(0)]
                izero_guinier[j] = math.exp(intercept)
                # I0Err
                izero_error_guinier[j] = math.exp(
                    intercept)*math.sqrt(reduced_chi_squared*correlation11)

                q_min_guinier[j] = scattering_data[j][start][0]
                q_max_guinier[j] = scattering_data[j][start +
                                                      points_to_use - 1][0]
                q_rg_min_guinier[j] = rg_guinier[j] * \
                    scattering_data[j][start][0]
                q_rg_max_guinier[j] = rg_guinier[j] * \
                    scattering_data[j][start + points_to_use - 1][0]
#                print('j, qmin, qmax, rg, qrgmin, qrgmax: ', j, q_min_guinier[j], q_max_guinier[j], rg_guinier[j], q_rg_min_guinier[j], q_rg_max_guinier[j])

                points_to_use += 1

            else:
                status = 1

        if(refine_scale_factor_flag == True):

            # scale factor = c[0]/C[j] = I0[0]*delta_rho_v[j]**2/(I0[j]/delta_rho_v[0]**2)
            scale_factor[j] = izero_guinier[0]/delta_rho_v[0] / \
                delta_rho_v[0]/(izero_guinier[j]/delta_rho_v[j]/delta_rho_v[j])
#            print(j, scale_factor[j])
#    print('refined scale factor: ', scale_factor)
# the refined scale factor gets passed back to the main program, no need to explicitely return it

    return points_used_guinier, chi_squared_guinier, rg_guinier, rg_error_guinier, izero_guinier, izero_error_guinier, q_min_guinier, q_max_guinier, q_rg_min_guinier, q_rg_max_guinier


if __name__ == "__main__":

#    import matplotlib.pyplot as plt
    import read_contrast_output_files_test as read_contrast_output_files
    import read_data_file_test as read_data_file

# input variables
    path = ('./')
    run_name = 'KinA:SDA_2:2'
    output_file_name = os.path.join(path, 'guinier_analysis_output_file.out')
    data_file_name = [os.path.join((path), '0.dat'), os.path.join((path), '10.dat'), os.path.join((path), '20.dat'), os.path.join(
        (path), '40.dat'), os.path.join((path), '80.dat'), os.path.join((path), '90.dat'), os.path.join((path), '100.dat')]
#    data_file_name = [os.path.join((path),'0.dat'),os.path.join((path),'100.dat')]
#    data_file_name = [os.path.join((path),'0.dat'),os.path.join((path),'40.dat')]
#    print('data_file_name: ', data_file_name)
    number_of_contrast_points = 7
#    number_of_contrast_points = 2
    print('number_of_contrast_points: ', number_of_contrast_points)
    fraction_d2o = [0.0, 0.1, 0.2, 0.4, 0.8, 0.9, 1.0]
#    fraction_d2o = [0.0, 1.0]
#    fraction_d2o = [0.0, 0.4]
    print('fraction D2O: ', fraction_d2o)
    number_of_components = 2  # hardwired at this time
    component_name = ['KinA', 'SDA']
    molecular_weight = [50.7, 11.7]  # kDa
    partial_specific_volume = [0.73, 0.73]
    q_rg_limit_guinier = 1.3
#    print('q_rg_limit_guinier: ', q_rg_limit_guinier)
    starting_data_point_guinier = [1, 1, 1, 1, 1, 1, 1]
    initial_points_to_use_guinier = [6, 6, 6, 6, 6, 6, 6]
    concentration = [11.9, 11.9, 11.9, 26.9, 11.9, 11.9, 11.9]
    # not currently used (but could be used to put an error on the scale factor)
    concentration_error = [0.6, 0.6, 0.6, 1.3, 0.6, 0.6, 0.6]
#    starting_data_point_guinier = [5, 10]
#    initial_points_to_use_guinier = [6, 6]
#    concentration = [11.9, 11.9]
#    concentration_error = [0.6, 0.6]
#    concentration = [11.9, 26.9]
#    concentration_error = [0.6, 1.3]
    read_from_contrast_calulator_output_file = False
    refine_scale_factor_flag = True

# initialization

# volume fraction:  volume = molecular_weight*partial_specific_volume/Na; volume fraction = v1/(v1+v2); Na cancels out and units of volume can be kept in kDa, as we just want the volume fraction
    total_volume = 0.0
    volume_fraction = numpy.zeros(number_of_components, float)
    for i in range(number_of_components):
        #        total_volume = total_volume + molecular_weight[i]*partial_specific_volume[i]
        total_volume += molecular_weight[i]*partial_specific_volume[i]
    for i in range(number_of_components):
        volume_fraction[i] = molecular_weight[i] * \
            partial_specific_volume[i]/total_volume
    print('volume fraction: ', volume_fraction)

    if read_from_contrast_calulator_output_file == True:
        contrast_file_name = os.path.join(path, 'mulch_test_data_contrast.txt')
        delta_rho = []
    else:
        contrast_file_name = ''
    print('contrast file name: ', contrast_file_name)

    if read_from_contrast_calulator_output_file == True:
        delta_rho = read_contrast_output_files.read_contrast_file(
            contrast_file_name, fraction_d2o, number_of_contrast_points, number_of_components)
    else:
        # the values below are exactly from the original decomposition input file
        #        delta_rho = [[2.270,6.227],[1.706,5.686],[1.141,5.146],[0.013,4.066],[-2.245,1.905],[-2.809,1.364],[-3.373,0.824]]
        # the values below are the same as those obtained from the contrast calculator output file
        delta_rho = [[2.31, 6.06], [1.75, 5.49], [1.18, 4.93], [
            0.055, 3.80], [-2.21, 1.54], [-2.78, 0.98], [-3.34, 0.41]]
#        delta_rho = [[2.270,6.227],[-3.373,0.824]]
#        delta_rho = [[2.270,6.227],[0.013,4.066]]
    print('delta_rho: ', delta_rho)

# Calculate an initial scale factor for each contrast based on the ratio of the concentrations.
# TODO: Test to make sure this works. The results are different when the scale factors aren't refined, especially if some of the concentrations are sufficiently different. Try making all concentrations the same first and see if the scale factor is refined like it was in the initial code.
# The scale factor is normalized to that of the first data set, i.e., scale_factor[i] = c[0]/c[i].
# TODO: Now the normalization is to the first data set by default. Allow the user to specify which data set?

# calculate the initial scale factor from the ratio of the concentrations.
    # want to keep track of the initial scale factor to write to output file
    initial_scale_factor = []
    scale_factor = []
    for i in range(number_of_contrast_points):
        ratio = concentration[0]/concentration[i]
        initial_scale_factor.append(ratio)
        scale_factor.append(ratio)
    print('initial scale factor from concentrations: ', scale_factor)

# need delta_rho_v = drho1*vf1 + drho2*vf2 for the scale factor refinement; I(0) = n*(delta_rho*volume)**2 is found from the Guinier analysis and then the scale factor is found using the ratio I(0)[0]/I(0)[i]; volume fraction is used here instead of volume since vf1 = v1/(v1+v2) but the denominator cancels out when calculating the ratio.

    delta_rho_v = numpy.zeros(number_of_contrast_points, float)
#    print('delta_rho_v: ', delta_rho_v)
#   delta_rho_v = vf1*delta_rho1 + vf2*delta_rho2
    for i in range(number_of_contrast_points):
        delta_rho_v[i] = volume_fraction[0]*delta_rho[i][0] + \
            volume_fraction[1]*delta_rho[i][1]
#        print(i, delta_rho_v[i])
#    print('delta_rho_v: ', delta_rho_v)

# read the scattering data files
    scattering_data = []
    number_of_data_lines = []

    for item in data_file_name:
        #        print('item: ', item)
        #        print('data file: ', item)
        q = []
        i = []
        ierr = []

        numpts, q, i, ierr = read_data_file.read_file(item)
#        print('numpts: ', numpts)
#        print('q: ', q)
#        print('i: ', i)
#        print('ierr: ', ierr)
        data = numpy.zeros((numpts, 3))
        for j in range(numpts):
            data[j][0] = q[j]
            data[j][1] = i[j]
            data[j][2] = ierr[j]
#        print(data)
        number_of_data_lines.append(numpts)
        scattering_data.append(data)

    print('number of data lines: ', number_of_data_lines)
#    print('scattering_data1: ', scattering_data)

# NOTE: conversion to numpy array doesn't seem necessary.  Results are the same if the next line is commented out.
    scattering_data = numpy.asarray(scattering_data)

# scattering_data is, with the sample data, a 7x182x3 array. For each pt (182) in each dataset (7), there are three values recorded: q, I(q), and Ierr
#    print('scattering_data2: ', scattering_data)

# the number of data points are the same for all files 
    number_of_data_points = number_of_data_lines[0]
    print('number_of_data_points: ', number_of_data_points)


# Guinier analysis
    print('scale factor before Guinier: ', scale_factor)

    points_used_guinier, chi_squared_guinier, rg_guinier, rg_error_guinier, izero_guinier, izero_error_guinier, q_min_guinier, q_max_guinier, q_rg_min_guinier, q_rg_max_guinier = guinier_fit(
        number_of_contrast_points, number_of_data_points, scattering_data, q_rg_limit_guinier, starting_data_point_guinier, initial_points_to_use_guinier, refine_scale_factor_flag, scale_factor, delta_rho_v)

    print('scale factor after Guinier: ', scale_factor)
    print('points used: ', points_used_guinier)
    print('chi squared: ', chi_squared_guinier)
    print('Rg, Rg_error', rg_guinier, rg_error_guinier)
    print('I(0), I(0) error: ', izero_guinier, izero_error_guinier)
    print('qmin, qmax: ', q_min_guinier, q_max_guinier)
    print('qRgmin, qRgmax: ', q_rg_min_guinier, q_rg_max_guinier)

# screen output
    print("Run name: " + run_name)
    print("\nGuinier analysis and scale factor calculation")
    print("fD2O\tRg\terrRg \tI(0)\terrI(0)\tqmin\tqmax\tqRgmin\tqRgmax\tX^2    #pts  Init Scale  Final Scale  delta-rho^2(10^-20 cm^-4)")
    for j in range(number_of_contrast_points):
        drho_sq = delta_rho_v[j]*delta_rho_v[j]
        print(str("%.2f" % fraction_d2o[j])+"\t"+str("%.2f" % rg_guinier[j])+"\t"+str("%.2f" % rg_error_guinier[j])+"\t"+str("%.3f" % izero_guinier[j])+"\t"+str("%.3f" % izero_error_guinier[j])+"\t"+str("%.4f" % q_min_guinier[j])+"\t"+str("%.4f" % q_max_guinier[j])+"\t"+str(
            "%.2f" % q_rg_min_guinier[j])+"\t"+str("%.2f" % q_rg_max_guinier[j])+"\t"+str("%.2f" % chi_squared_guinier[j])+"\t"+str("%2i" % points_used_guinier[j])+"\t"+str("%.2f" % initial_scale_factor[j])+"\t   "+str("%.2f" % scale_factor[j])+"\t\t"+str("%.3f" % drho_sq))


# output file
    with io.open(output_file_name, 'w') as outfile:
        outfile.write("Run name: "+run_name+"\n")
        outfile.write("\nGuinier analysis and scale factor calculation\n")
        outfile.write(
            "fD2O\tRg\terrRg \tI(0)\terrI(0)\tqmin\tqmax\tqRgmin\tqRgmax\tX^2    #pts  Init Scale  Final Scale  delta-rho^2(10^-20 cm^-4)\n")
        for j in range(number_of_contrast_points):
            drho_sq = delta_rho_v[j]*delta_rho_v[j]
            outfile.write(str("%.2f" % fraction_d2o[j])+"\t"+str("%.2f" % rg_guinier[j])+"\t"+str("%.2f" % rg_error_guinier[j])+"\t"+str("%.3f" % izero_guinier[j])+"\t"+str("%.3f" % izero_error_guinier[j])+"\t"+str("%.4f" % q_min_guinier[j])+"\t"+str("%.4f" % q_max_guinier[j])+"\t"+str(
                "%.2f" % q_rg_min_guinier[j])+"\t"+str("%.2f" % q_rg_max_guinier[j])+"\t"+str("%.2f" % chi_squared_guinier[j])+"\t"+str("%2i" % points_used_guinier[j])+"\t"+str("%.2f" % initial_scale_factor[j])+"\t   "+str("%.2f" % scale_factor[j])+"\t\t"+str("%.3f" % drho_sq) + '\n')
    outfile.close()
