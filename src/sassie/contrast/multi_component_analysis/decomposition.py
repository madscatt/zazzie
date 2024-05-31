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
#    DECOMPOSITION is the method that calculates the composite scattering
#    functions, I11, I12, and I22 for a two-component complex from a set of
#    contrast variation data. I11 and I22 are the scattering intensities of
#    components 1 and 2, respectively, and I12 is the cross-term.
#    A Guinier analysis is performed at each contrast to rescale the data if
#    desired by the user. The composite scattering functions are then found
#    using multiple linear regression at each q value:
#    I=delta_rho_1**2*I11+delta_rho_1*delta_rho_2*I12+delta_rho_2**2*I22
#
#
#       5/25/2022       --  Initial coding              :   Kathryn Sarachan
#       4/23/2023       --  Python 2 code finished      :  Kathryn Sarachan
#       4/25/2023       --  Converted to Python 3       :   Susan Krueger
#       5/23/2023       --  SASSIE 3.0 test program     :   Susan Krueger
#       6/27/2023       --  SASSIE 3.0 production       :   Susan Krueger
#       2/2024          --  Errors added to model data  :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''

    **Decomposition** contains the **Get Composite Scattering Intensities** method and related functions that calculate the scattering intensities of the individual components in a two component complex and the cross-term intensity, which gives an idea of the spatial relationship between the two components, from SANS contrast variation data. The original code was written by Andrew Whitten (2/2006) as part of the MULCh program. Rewritten in Python by Kathryn Sarachan (4/2023). 
    
    **Reference:** Whitten, A.E., Cai, S, Trewhella, J. (2008). "MULCh: modules for the analysis of small-angle neutron contrast variation data from biomolecular assemblies", *J. Appl. Cryst.* **41**, 222 - 226. 
    
    **Inputs:**
    
    SANS data, fraction D\ :sub:`2`\ O for each contrast
    
    contrast values (delta_rho) for each component at each contrast
    
    **Outputs:**
    
    Intensities of the individual components, cross-term intensity
    
    Called by the **Multi-component Analysis** module.

    Requires **Get Composite Scattering Intensities**, **Read Data File**, **Gunier Analysis**, **Do Multiple Linear Regression**, **Calculate Reduced Chi-squared** and **numpy.linalg.inv**
        
'''

import io
import time
import numpy
import sassie.contrast.multi_component_analysis.read_sas_data_file_add_error as read_data_file
import sassie.contrast.multi_component_analysis.guinier_analysis as guinier_analysis
#import read_sas_data_file_add_error as read_data_file
#import guinier_analysis as guinier_analysis




def save_data_to_plot_as_json(other_self, square_root_izero, square_root_izero_error, square_root_izero_calculated, diff, match_point, match_point_error):

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    data_dict = {
        'fraction_d2o': [],
        'sqrt[I(0)/c]': [],
        'sqrt[I(0)/c]_error': [],
        'sqrt[I(0)/c]_calc': [],
        'sqrt[I(0)/c]-sqrt[I(0)/c]_calc': [],
        'match_point': [],
        'match_point_error': [],
        'match_point_y_value': []
    }

    for i in range(mvars.number_of_contrast_points):
        # Append new values to the lists in the dictionary
        data_dict['fraction_d2o'].append(mvars.fraction_d2o[i])
        data_dict['sqrt[I(0)/c]'].append(square_root_izero[i])
        data_dict['sqrt[I(0)/c]_error'].append(square_root_izero_error[i])
        data_dict['sqrt[I(0)/c]_calc'].append(square_root_izero_calculated[i])
        data_dict['sqrt[I(0)/c]-sqrt[I(0)/c]_calc'].append(diff[i])

    data_dict['match_point'].append(match_point)
    data_dict['match_point_error'].append(match_point_error)
    data_dict['match_point_y_value'].append(0)

    json_data = json.dumps(data_dict)

    mcavars.json_outfile.write(json_data)
    mcavars.json_outfile.close()

    return


def get_composite_scattering_intensities(other_self):
    r'''

    **Get Composite Scattering Intensities** is the **Decomposition Analysis** method that calculates the composite scattering functions, :math:`I_{11}`, :math:`I_{12}`, and :math:`I_{22}` for a two-component complex from a set of SANS contrast variation data. :math:`I_{11}`, and :math:`I_{22}` are the scattering intensities of components 1 and 2, respectively, and :math:`I_{12}` is the cross-term. First, the data are rescaled with respect to each other to account for differences in concentration as a function of contrast (fraction D\ :sub:`2`\ O  in the solvent).  Then, a Guinier analysis is performed to obtain :math:`R_{g}` and :math:`I_{0}` as a function of contrast.  If desired by the user, the :math:`I_{0}` values are used to further refine the scale factors at each contrast to account for possible inaccuracies in the concentrations. See the **Guinier Analysis** helper program for details.  

    Executed if decomposition_flag == True.

    Notes:

        Output file with chosen output_file_name is stored in the output file directory.
        TODO: Need boilerplate language to add variable name/path for the directory. 

        Output file contains:
            - module variable input parameters
            - results from Guinier analysis at each fraction D\ :sub:`2`\ O value
                - Guinier R\ :sub:`g`\ , Guinier R\ :sub:`g`\  error, Guinier I(0), Guinier I(0) error
                - q\ :sub:`min`\ , q\ :sub:`max`\ , qR\ :sub:`gmin`\ , qR\ :sub:`gmax`\ 
                - number of points used, reduced chi-squared, :math:`\Delta \rho^2` 
                - initial and final scale factors
            - results from decomposition analysis
                - mean square differences at each q value for all fraction D\ :sub:`2`\ O values
                - global reduced chi-squared value at each q value

        Additional output files stored in the output file directory:
            - composition scattering intensities (3 files)
            - rescaled data files (number of files = number of contrast points)
            - calculated data files (number of files = number of contrast points)

    This method calls **Read Data File**, **Gunier Analysis**, **Do Multiple Linear Regression** and **Calculate Reduced Chi-squared**.


    Parameters
    ----------

    run_name:  string
        Name of directory containing the outputs

    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values

    fraction_d2o:   float array (dimension = number_of_contrast_points)
        The fraction D\ :sub:`2`\ O values that define the contrasts

    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )

    data_file_name:  string array (dimension = number_of_contrast_points)
        The names of the contrast variation data files

    q_rg_limit_guinier: float array (dimension = number_of_contrast_points)
        qR\ :sub:`g`\  limit for the Guinier analysis at each fraction D\ :sub:`2`\ O

    starting_data_point_guinier: int array (dimension = number_of_contrast_points)
        The index of the starting data point for the Guinier fit at each fraction D\ :sub:`2`\ O  (index of the first data point = 1)

    initial_points_to_use_guinier: int array (dimension = number_of_contrast_points)
        The number of data points to use initially for the Guinier fit at each fraction D\ :sub:`2`\ O  (the final number of points used depends on the qR\ :sub:`g`\  limit)

    initial_guess_guinier: float array (dimension = 2 for a line fit)
        The initial guess for the Guinier fit parameters (default = [1., 1.])

    refine_scale_factor_flag: boolean
        Indicates whether the scale factor at each fraction D\ :sub:`2`\ O  will be adjusted based on the I(0) values calculated from the Guinier fit and delta_rho_v

    multi_component_analysis_path: string
        sub-path where output file will be written: run_name + \'multi_component_analysis\' + method-dependent sub-path

    outfile: string
        output file name (with full path): path + output_file_name

    initial_scale_factor:  float array (dimension = number_of_contrast_points)
        initial scale factor for the data at each fraction D\ :sub:`2`\ O 

    scale_factor:  float array (dimension = number_of_contrast_points)
        scale factor for the data at each fraction D\ :sub:`2`\ O  that is the same as the initial scale factor before the Guinier analysis is performed

    delta_rho_v: float array (dimension = number_of_contrast_points)
        :math:`\Delta \rho V` as defined above at each fraction D\ :sub:`2`\ O  as defined in the Guinier analysis helper program

    composite_intensity_file_name:  string array (dimension = 3)
        names of the composite scattering intensity output files

    rescaled_data_file_name:  string array (dimension = number_of_contrast_points)
        names of the rescaled data output files

    calculated_data_file_name:  string array (dimension = number_of_contrast_points)
        names of the calculated data output files

    signal_to_noise_amplitude: float array (dimension = number_of_contrast_points)
        amplitude of the Gaussian equation that describes the signal-to-noise (S/N) vs q behavior of SANS data; used when adding noise to model SANS data

    signal_to_noise_mean: float
        mean of the Gaussian equation that describes the signal-to-noise (S/N) vs q behavior of SANS data; used when adding noise to model SANS data    

    signal_to_noise_standard_deviation: float
        standard deviation of the Gaussian equation that describes the signal-to-noise (S/N) vs q behavior of SANS data; used when adding noise to model SANS data

    signal_to_noise_background: float
        background term of the Gaussian equation that describes the signal-to-noise (S/N) vs q behavior of SANS data; used when adding noise to model SANS data  


    Returns
    -------

    rg_guinier: float array( dimension = number_of_contrast_points)
        The radius of gyration from the Guinier fit at each fraction D\ :sub:`2`\ O
        
    rg_guinier_error: float array (dimension = number_of_contrast_points)
        The error in the Gunier radius of gyration at each fraction D\ :sub:`2`\ O

    izero_guinier: float array 
        The I(0) value from the Guinier fit at each fraction D\ :sub:`2`\ O 

    izero_error_guinier: float array
        The error in the Guinier I(0) value at each fraction D\ :sub:`2`\ O 

    points_used_guinier: int array (dimension = number_of_contrast_points)
        The final number of points used in the Guinier fit at each fraction D\ :sub:`2`\ O 

    chi_squared_guinier: float array (dimension = number_of_contrast_points)
        The chi-square value for the Guinier fit at each fraction D\ :sub:`2`\ O 

    q_min_guinier: float array (dimension = number_of_contrast_points)
        The minimum q value for the Guinier fit at each fraction D\ :sub:`2`\ O 

    q_max_guinier: float array (dimension = number_of_contrast_points)
        The maximum q value for the Guinier fit at each fraction D\ :sub:`2`\ O 

    q_rg_min_guinier: float array (dimension = number_of_contrast_points)
        The minimum qR\ :sub: `g`\  value for the Guinier fit at each fraction D\ :sub:`2`\ O 
        
    q_rg_max_guinier: float array (dimension = number_of_contrast_points)
        The maximum qR\ :sub: `g`\  value for the Guinier fit at each fraction D\ :sub:`2`\ O 

    scale_factor: float array (dimension = number_of_contrast_points)
        The final scale factor at each fraction D\ :sub:`2`\ O  (rescaling is only preformed if refine_scale_factor_flag is True) 

    rescaled_scattering_data:  3D array (dimensions = number_of_contrast_points x number_of_data_points x 3)
        q, I(q) and I(q) error for the rescaled data at at each fraction D\ :sub:`2`\ O

    mean_square_difference:  2D float array (dimensions = number_of_data_points x number_of_contrast_points)
        a complete list (all q values, all datasets) of mean square differences

    reduced_chi_squared_list: float array (dimension = number_of_data_points)
        global reduced chi-squared at all q values

    composite_scattering_intensity:  3D float array (dimensions = 3 x number_of_data_points x 2)
        I(q) and I(q) error for each of the 3 composite scattering intensities, :math:`I_{11}`, :math:`I_{12}` and :math:`I_{22}`

    calculated_scattering_data:  2D float array (dimensions = number_of_data_points x number_of_contrast_points)
        calculated I(q) values from the above :math:`I_{exp}` equation at each fraction D\ :sub:`2`\ O

    '''

# mvars used:  run_name, fraction_d2o, number_of_contrast_points, delta_rho, data_file_name, q_rg_limit_guinier, starting_data_point_guinier, initial_points_to_use_guinier, initial_guess_guinier, refine_scale_factor_flag, signal_to_noise_amplitude, signal_to_noise_mean, signal_to_noise_standard_deviation, signal_to_noise_background
# mcavars used:  outfile, initial_scale_factor, scale_factor, delta_rho_v, composite_intensity_file_name, rescaled_data_file_name, calculated_data_file_name, multi_component_analysis_path

    log = other_self.log
    log.debug('\nin get_composite_scattering_intensities\n')
    pgui = other_self.run_utils.print_gui

    mvars = other_self.module_variables
    mcavars = other_self.multi_component_analysis_variables

    log.debug('variables before decomposition analysis\n')
    log.debug(vars(mvars))
    log.debug(vars(mcavars))

    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in range(60)])

    pgui('\n%s \n' % (st))
    pgui('DATA FROM RUN: %s \n\n' % (ttxt))

    pgui('Decomposition method\n')
    pgui('reading the contrast variation data files\n')

# read the scattering data files
    mcavars.scattering_data = []
    number_of_data_lines = []
    error_flag_message = []
# NOTE:  error_flag_message indicates whether errors had to be added (or modified in the case of an error of 0.0) to the data, as in the case of model data.  If the message is blank, then no changes were made to the errors.  The errors must be nonzero since weighted fits are performed.
# TODO:  Should errors on the parameters be reported in the case of model data?  Need to test on model data to determine how changing the errors affect the results.


# loop to read multiple files for contrast variation data
# NOTE: the files are checked in the module filter to make sure that they are data files containing at least q and I(q); if no errors on I(q) are found, errors are added as described in read_data_file

    for i in range(mvars.number_of_contrast_points):
        q = []
        iq = []
        ierr = []

        numpts, message, q, iq, ierr = read_data_file.read_file(
            mvars.data_file_name[i], mvars.signal_to_noise_amplitude[i], mcavars.signal_to_noise_mean, mcavars.signal_to_noise_standard_deviation, mcavars.signal_to_noise_background)

#        print('numpts: ', numpts)
#        print('q: ', q)
#        print('iq: ', iq)
#        print('ierr: ', ierr)
        data = numpy.zeros((numpts, 3))
        for j in range(numpts):
            data[j][0] = q[j]
            data[j][1] = iq[j]
            data[j][2] = ierr[j]
#        print(data)
        number_of_data_lines.append(numpts)
        mcavars.scattering_data.append(data)
        error_flag_message.append(message)

    log.debug('number of data lines: ' + str(number_of_data_lines) + '\n')
    log.debug('error_flag_message: ' + str(error_flag_message) + '\n')
#    print('scattering_data1: ', mcavars.scattering_data)

# NOTE: conversion to numpy array (from Katie's original version) doesn't seem necessary.  Results are the same if the next line is commented out.
    mcavars.scattering_data = numpy.asarray(mcavars.scattering_data)

# scattering_data is a 7x182x3 array. For each pt (182) in each dataset (7), there are three values recorded: q, I(q), and Ierr
#    print('scattering_data2: ', mcavars.scattering_data)

# Since the number of data points are the same for all files (which was checked in the module filter), then set the scalar value for the module variable, number_of_data_points
    mcavars.number_of_data_points = number_of_data_lines[0]
    log.debug('number_of_data_points: ' +
              str(mcavars.number_of_data_points) + '\n')

# Guinier analysis
#    print('scale factor before Guinier: ', mcavars.scale_factor)

# new mcavars after Guinier fit:  points_used_guinier, chi_squared_guinier, rg_guinier, rg_error_guinier, izero_guinier, izero_error_guinier, q_min_guinier, q_max_guinier, q_rg_min_guinier, q_rg_max_guinier = guinier_analysis.guinier_fit(
# vars used in Guinier fit:  mvars.number_of_contrast_points, mcavars.number_of_data_points, mcavars.scattering_data, mvars.q_rg_limit_guinier, mvars.starting_data_point_guinier, mvars.initial_points_to_use_guinier, mvars.refine_scale_factor_flag, mcavars.scale_factor, mcavars.delta_rho_v, mvars.initial_guess_guinier)
    guinier_analysis.guinier_fit(other_self)

#    print('scale factor after Guinier: ', mcavars.scale_factor)
#    print('points used: ', mcavars.points_used_guinier)
#    print('chi squared: ', mcavars.chi_squared_guinier)
#    print('Rg, Rg_error', mcavars.rg_guinier, mcavars.rg_error_guinier)
#    print('I(0), I(0) error: ', mcavars.izero_guinier,
#          mcavars.izero_error_guinier)
#    print('qmin, qmax: ', mcavars.q_min_guinier, mcavars.q_max_guinier)
#    print('qRgmin, qRgmax: ', mcavars.q_rg_min_guinier, mcavars.q_rg_max_guinier)

    log.debug('\nscale factor after Guinier ' +
              str(mcavars.scale_factor) + '\n')

    pgui('calculating the composite scattering intensities\n')
# Multiple linear regression

# rescale the data with respect to the input concentrations; if the input concentrations are similar at each contrast, this rescaling doesn't matter that much.  But, if one or more are significantly different, rescaling will give a more accurate result.
    mcavars.rescaled_scattering_data = numpy.zeros(
        numpy.shape(mcavars.scattering_data))
    # print('shape of rescaled scattering data array: ',numpy.shape(mcavars.rescaled_scattering_data))
    for i in range(mvars.number_of_contrast_points):
        for j in range(mcavars.number_of_data_points):
            mcavars.rescaled_scattering_data[i][j][0] = mcavars.scattering_data[i][j][0]
            mcavars.rescaled_scattering_data[i][j][1] = mcavars.scale_factor[i] * \
                mcavars.scattering_data[i][j][1]
            mcavars.rescaled_scattering_data[i][j][2] = mcavars.scale_factor[i] * \
                mcavars.scattering_data[i][j][2]

# mean_square_difference is a list of lists containing all mean square differences. composite_scattering_intensity is a numpy array containing I(q) and sigmaI(q) for all the composite scattering functions, and reduced_chi_squared_list is a list of reduced_chi_squared values for each q
    mean_square_difference, composite_scattering_intensity, calculated_scattering_data, reduced_chi_squared_list = do_multiple_linear_regression(
        mvars.delta_rho, mcavars.rescaled_scattering_data, mcavars.number_of_data_points, mvars.number_of_contrast_points)

#    for i in range(mcavars.number_of_data_points):
#        print('i, calculated scattering data final: ', i, calculated_scattering_data[i])


# Output for the user

# screen output

    pgui('results written to output file: %s' %
         (mcavars.multi_component_analysis_path+mvars.output_file_name))
    pgui('component 1 composite scattering intensity written to output file: %s' %
         (mcavars.multi_component_analysis_path+mcavars.composite_intensity_file_name[0]))
    pgui('component 2 composite scattering intensity written to output file: %s' %
         (mcavars.multi_component_analysis_path+mcavars.composite_intensity_file_name[2]))
    pgui('inter-component composite scattering intensity written to output file: %s' %
         (mcavars.multi_component_analysis_path+mcavars.composite_intensity_file_name[1]))
    pgui('rescaled data files written to: %s' %
         (mcavars.multi_component_analysis_path))
    pgui('calculated data files written to: %s' %
         (mcavars.multi_component_analysis_path))
    pgui('-------------------------------')
    pgui('Final Results\n')

    pgui("\nDecomposition: Module for extracting composite scattering functions from contrast variation data (v2 2024)")
    pgui("Run name: " + mvars.run_name)
    pgui("\nGuinier analysis and scale factor calculation")
    pgui("fD2O\tRg\terrRg \tI(0)\terrI(0)\tqmin\tqmax\tqRgmin\tqRgmax\tX^2    #pts  Init Scale  Final Scale  delta-rho^2(10^-20 cm^-4)")
    for j in range(mvars.number_of_contrast_points):
        drho_sq = mcavars.delta_rho_v[j]*mcavars.delta_rho_v[j]
        pgui(str("%.2f" % mvars.fraction_d2o[j])+"\t"+str("%.2f" % mcavars.rg_guinier[j])+"\t"+str("%.2f" % mcavars.rg_error_guinier[j])+"\t"+str("%.3f" % mcavars.izero_guinier[j])+"\t"+str("%.3f" % mcavars.izero_error_guinier[j])+"\t"+str("%.4f" % mcavars.q_min_guinier[j])+"\t"+str("%.4f" % mcavars.q_max_guinier[j])+"\t"+str(
            "%.2f" % mcavars.q_rg_min_guinier[j])+"\t"+str("%.2f" % mcavars.q_rg_max_guinier[j])+"\t"+str("%.2f" % mcavars.chi_squared_guinier[j])+"\t"+str("%2i" % mcavars.points_used_guinier[j])+"\t"+str("%.2f" % mcavars.initial_scale_factor[j])+"\t   "+str("%.2f" % mcavars.scale_factor[j])+"\t\t"+str("%.3f" % drho_sq))

# this output is long; do we want to output it to the screen since the user will have to scroll up quite a bit to see the results of the Guinier analysis?
    
    outstring = ""
    outstring1 = ""
    outstring2 = ""
    pgui("\nDeviations between the composite scattering functions and the contrast variation series") 
    outstring="fD2O:"
    for i in range(mvars.number_of_contrast_points):
        outstring = outstring+"\t"+str("%.2f" %mvars.fraction_d2o[i])
    outstring = outstring+"\n"    
    pgui(outstring)

    outstring1 = "q\t(delta_I(q)/sigma(I(q)))^2"
    tabs = mvars.number_of_contrast_points - 3
    if(tabs > 1):
        for i in range(tabs):
            outstring1 = outstring1+"\t"
    else:
        outstring1 = "q\t(delta_I(q)/sigma(I(q)))^2\t"
    outstring1 = outstring1+"chi^2"
    pgui(outstring1)
    for j in range(mcavars.number_of_data_points):
        outstring2 = str("%.4f" %mcavars.scattering_data[0][j][0])+"\t"
        for i in range(mvars.number_of_contrast_points):
            outstring2 = outstring2+mean_square_difference[j][i]+"\t"
        pgui(outstring2+str("%.2f" %reduced_chi_squared_list[j]))
    

# output to files

# general output file
    mcavars.outfile.write('--------------------------------\n')
    mcavars.outfile.write('Final Results\n')
#    with io.open(output_file_path+output_file_name, 'w') as outfile:
    mcavars.outfile.write("Run name: "+mvars.run_name+"\n")
    mcavars.outfile.write("\nGuinier analysis and scale factor calculation\n")
    mcavars.outfile.write(
        "fD2O\tRg\terrRg \tI(0)\terrI(0)\tqmin\tqmax\tqRgmin\tqRgmax\tX^2    #pts  Init Scale  Final Scale  delta-rho^2(10^-20 cm^-4)\n")
    for j in range(mvars.number_of_contrast_points):
        drho_sq = mcavars.delta_rho_v[j]*mcavars.delta_rho_v[j]
        mcavars.outfile.write(str("%.2f" % mvars.fraction_d2o[j])+"\t"+str("%.2f" % mcavars.rg_guinier[j])+"\t"+str("%.2f" % mcavars.rg_error_guinier[j])+"\t"+str("%.3f" % mcavars.izero_guinier[j])+"\t"+str("%.3f" % mcavars.izero_error_guinier[j])+"\t"+str("%.4f" % mcavars.q_min_guinier[j])+"\t"+str("%.4f" % mcavars.q_max_guinier[j])+"\t"+str(
            "%.2f" % mcavars.q_rg_min_guinier[j])+"\t"+str("%.2f" % mcavars.q_rg_max_guinier[j])+"\t"+str("%.2f" % mcavars.chi_squared_guinier[j])+"\t"+str("%2i" % mcavars.points_used_guinier[j])+"\t"+str("%.2f" % mcavars.initial_scale_factor[j])+"\t   "+str("%.2f" % mcavars.scale_factor[j])+"\t\t"+str("%.3f" % drho_sq) + '\n')
    mcavars.outfile.write(
        "\nDeviations between the composite scattering functions and the contrast variation series\n")

    outstring = ""
    outstring1 = ""
    outstring2 = ""
    outstring = "fD2O:"
    for i in range(mvars.number_of_contrast_points):
        outstring = outstring+"\t"+str("%.2f" % mvars.fraction_d2o[i])
    outstring = outstring+"\n\n"
#       print('outstring: ', outstring)
    mcavars.outfile.write(outstring)

    outstring1 = "q\t(delta_I(q)/sigma(I(q)))^2"
    tabs = mvars.number_of_contrast_points - 3
    if (tabs > 1):
        for i in range(tabs):
            outstring1 = outstring1+"\t"
    else:
        outstring1 = "q\t(delta_I(q)/sigma(I(q)))^2\t\t"
    outstring1 = outstring1+"chi^2\n"
#       print('outstring1: ', outstring1)
    mcavars.outfile.write(outstring1)
    for j in range(mcavars.number_of_data_points):
        outstring2 = str("%.4f" % mcavars.scattering_data[0][j][0])+"\t"
        for i in range(mvars.number_of_contrast_points):
            outstring2 = outstring2+mean_square_difference[j][i]+"\t"
        mcavars.outfile.write(outstring2+str("%.2f" %
                                             reduced_chi_squared_list[j])+"\n")
    mcavars.outfile.close()


# rescaled data files
# NOTE:  The rescaling will either be the initial rescaling based on the inputted concentration values or the refined scale factors.  This is indicated in the output file, along with the scale factor used.  It is also indicated that the scale factor is wrt the 1st data set and that file name is written out.
    for i in range(mvars.number_of_contrast_points):
        #        print('rescaled data file: ', rescaled_data_file_name[i])
        with io.open(mcavars.multi_component_analysis_path+mcavars.rescaled_data_file_name[i], 'w') as rfile:
            rfile.write("# rescaled file for: "+mvars.data_file_name[i]+"\n")
            if (mvars.refine_scale_factor_flag == True):
                rfile.write(
                    "# rescaled based on refined concentrations with respect to file: "+mvars.data_file_name[0]+"\n")
            else:
                rfile.write(
                    "# rescaled based on input concentrations with respect to file: "+mvars.data_file_name[0]+"\n")
            rfile.write("# scale factor: "+str(mcavars.scale_factor[i])+"\n")
            if (error_flag_message[i] != " "):
                rfile.write("# " + error_flag_message[i])
                if (error_flag_message[i][0:19] != "Error values of 0.0"):
                    rfile.write("#Amplitude of Gaussian used to describe S/N vs q: %f\n" %
                                (mvars.signal_to_noise_amplitude[i]))
            rfile.write("#     q           I(q)            sigma(I(q))\n")
            for j in range(mcavars.number_of_data_points):
                rfile.write(
                    str("%.6f" % mcavars.rescaled_scattering_data[i][j][0])+"\t\t")
                rfile.write(
                    str("%.6f" % mcavars.rescaled_scattering_data[i][j][1])+"\t\t")
                rfile.write(
                    str("%.6f" % mcavars.rescaled_scattering_data[i][j][2])+"\n")
        rfile.close()

# calculated data files
# NOTE: calculated data are obtained from the composite scattering intensities
# TODO: plot the calculated data along with the rescaled data and the mean square differences -- or do we want to plot the residuals?
    for i in range(mvars.number_of_contrast_points):
        # print('calculated data file: ', mcavars.calculated_data_file_name[i])
        with io.open(mcavars.multi_component_analysis_path+mcavars.calculated_data_file_name[i], 'w') as cfile:
            cfile.write("# calculated file for: " +
                        mcavars.rescaled_data_file_name[i]+"\n")
            if (mvars.refine_scale_factor_flag == True):
                cfile.write(
                    "# rescaled based on refined concentrations with respect to file: "+mvars.data_file_name[0]+"\n")
            else:
                cfile.write(
                    "# rescaled based on input concentrations with respect to file: "+mvars.data_file_name[0]+"\n")
            cfile.write("# scale factor: "+str(mcavars.scale_factor[i])+"\n")
            cfile.write("#     q           I(q)            sigma(I(q))\n")
            for j in range(mcavars.number_of_data_points):
                cfile.write(
                    str("%.6f" % mcavars.rescaled_scattering_data[0][j][0])+"\t\t")
                cfile.write(
                    str("%.6f" % calculated_scattering_data[j][i])+"\t\t")
                cfile.write("0.000000\n")
        cfile.close()

# write out the composite scattering intensities
# TODO:  This could be done in a loop

    with io.open(mcavars.multi_component_analysis_path+mcavars.composite_intensity_file_name[0], 'w') as f1:
        f1.write("# Component 1 composite scattering intensity\n")
        f1.write("#     q           I(q)         sigma(I(q))\n")
        for i in range(mcavars.number_of_data_points):
            f1.write(str("%.6f" % mcavars.scattering_data[0][i][0])+"\t\t")
            f1.write(
                str("%.6f" % composite_scattering_intensity[0][i][0])+"\t\t")
            f1.write(
                str("%.6f" % composite_scattering_intensity[0][i][1])+"\n")
    f1.close()
    with io.open(mcavars.multi_component_analysis_path+mcavars.composite_intensity_file_name[1], 'w') as f2:
        f2.write("# Inter-component composite scattering intensity\n")
        f2.write("#     q           I(q)         sigma(I(q))\n")
        for i in range(mcavars.number_of_data_points):
            f2.write(str("%.6f" % mcavars.scattering_data[0][i][0])+"\t\t")
            f2.write(
                str("%.6f" % composite_scattering_intensity[1][i][0])+"\t\t")
            f2.write(
                str("%.6f" % composite_scattering_intensity[1][i][1])+"\n")
    f2.close()
    with io.open(mcavars.multi_component_analysis_path+mcavars.composite_intensity_file_name[2], 'w') as f3:
        f3.write("# Component 2 composite scattering intensity\n")
        f3.write("#     q           I(q)         sigma(I(q))\n")
        for i in range(mcavars.number_of_data_points):
            f3.write(str("%.6f" % mcavars.scattering_data[0][i][0])+"\t\t")
            f3.write(
                str("%.6f" % composite_scattering_intensity[2][i][0])+"\t\t")
            f3.write(
                str("%.6f" % composite_scattering_intensity[2][i][1])+"\n")
    f3.close()


#HERE

    #save_data_to_plot_as_json(other_self, square_root_izero, square_root_izero_error,
    #                          square_root_izero_calculated, diff, match_point, match_point_error)

    time.sleep(0.5)

    return


# TODO: this method can be written in terms of mvars and mcavars,just passing the index, q?


def calculate_reduced_chi_squared(delta_rho, rescaled_scattering_data, composite_scattering_intensity, number_of_contrast_points, q):
    '''
    **Get Reduced Chi Squared** calculates the **global** reduced_chi_squared value for all SANS contrast variation datasets at a particular q value.

    Parameters
    ----------

    q: int
        The index of the q value for which the calculation is being performed
    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
     rescaled_scattering_data:  3D array (dimensions = number_of_contrast_points x number_of_data_points x 3)
        q, I(q) and I(q) error for the rescaled data at at each fraction D\ :sub:`2`\ O
    composite_scattering_intensity:  3D float array (dimensions = 3 x number_of_data_points x 2)
        I(q) and I(q) error for each of the 3 composite scattering intensities, :math:`I_{11}`, :math:`I_{12}` and :math:`I_{22}`

    Returns
    -------

    msd_list: float array (dimension = number_of_contrast_points)
        mean square difference at each fraction D\ :sub:`2`\ O value  at the desired q value
    calculated_intensity: float array (dimension = number_of_contrast_points)
        calculated scattering intensity at each fraction D\ :sub:`2`\ O value  at the desired q value
    reduced_chi_squared: float
        global reduced chi-squared at the desired q value

    '''

    delta_rho_xy = numpy.zeros(3)
    nSize = 0
    reduced_chi_squared = 0.0
    calculated_intensity = numpy.zeros(number_of_contrast_points)
    msd_list = []

    for i in range(number_of_contrast_points):
        delta_rho_xy[0] = delta_rho[i][0]*delta_rho[i][0]
        delta_rho_xy[1] = delta_rho[i][0]*delta_rho[i][1]
        delta_rho_xy[2] = delta_rho[i][1]*delta_rho[i][1]
#        intensity = scale_factor[i]*scattering_data[i][q][1]
#        intensity_error = scale_factor[i]*scattering_data[i][q][2]
        intensity = rescaled_scattering_data[i][q][1]
        intensity_error = rescaled_scattering_data[i][q][2]

# TODO:  Don't need to check for zero error here because it is checked in the module filter

        if intensity_error > 0:
            nSize += 1
            w = 1/intensity_error
            calc_intensity = delta_rho_xy[0]*composite_scattering_intensity[0][q][0] + delta_rho_xy[1] * \
                composite_scattering_intensity[1][q][0] + \
                delta_rho_xy[2]*composite_scattering_intensity[2][q][0]
#            D = intensity - delta_rho_xy[0]*composite_scattering_intensity[0][q][0] - delta_rho_xy[1]*composite_scattering_intensity[1][q][0] - delta_rho_xy[2]*composite_scattering_intensity[2][q][0]
            D = intensity - calc_intensity
            if (w*D*w*D < 9.0):
                MSD = str("%.2f" % (w*D*w*D))
            else:
                MSD = str("%.2f" % (w*D*w*D))+"*"
            reduced_chi_squared = reduced_chi_squared + w*w*D*D
            msd_list.append(MSD)
            calculated_intensity[i] = calc_intensity

    if nSize <= 3:
        reduced_chi_squared = 0.0
    else:
        reduced_chi_squared = reduced_chi_squared / (float(nSize) - 3.0)

#    print('q, nsize, msd list, reduced chi squared: ', q, nSize, msd_list, reduced_chi_squared)
#    print('q, calculated intensity: ', q, calculated_intensity)
#    print('q, nsize, reduced chi squared: ', q, nSize, reduced_chi_squared)
#    print('msd list: ', msd_list)
    return msd_list, reduced_chi_squared, calculated_intensity


# TODO: this method can be written in terms of mvars and mcavars
def do_multiple_linear_regression(delta_rho, rescaled_scattering_data, number_of_data_points, number_of_contrast_points):
    r'''

    **Do Multiple Linear Regression** performs a multiple linear regression **at each q value** to solve the set equations defined at each contrast:

    :math:`I_{exp}=\Delta \rho_1^2 I_{11} + \Delta \rho_1 \Delta \rho_2 I_{12} + \Delta \rho_2^2 I_{22}`

    where :math:`I_{exp}` is the **rescaled** measured scattering data and :math:`\Delta \rho_1` and :math:`\Delta \rho_2` are the contrasts of components 1 and 2, respectively. 

    The matrix inversion method used is as described, for example, in Data Analysis and Error Reduction for the Physical Sciences, by Philip R. Bevington and D. Keith Robinson, McGraw-Hill, New York, NY 2003, 1992, 1969 (Chapter 7).  

    Parameters
    ----------

    number_of_contrast_points:  int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    number_of_data_points: int
        number of data points at each fraction D\ :sub:`2`\ O  (all files have the same number of data points)
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )
    rescaled_scattering_data:  3D array (dimensions = number_of_contrast_points x number_of_data_points x 3)
        q, I(q) and I(q) error for the rescaled data at at each fraction D\ :sub:`2`\ O

    Returns
    -------

    mean_square_difference:  2D float array (dimensions = number_of_data_points x number_of_contrast_points)
        a complete list (all q values, all datasets) of mean square differences
    reduced_chi_squared_list: float array (dimension = number_of_data_points)
        global reduced chi-squared at all q values
    composite_scattering_intensity:  3D float array (dimensions = 3 x number_of_data_points x 2)
        I(q) and I(q) error for each of the 3 composite scattering intensities, :math:`I_{11}`, :math:`I_{12}` and :math:`I_{22}`
    calculated_scattering_data:  2D float array (dimensions = number_of_data_points x number_of_contrast_points)
        calculated I(q) values from the above :math:`I_{exp}` equation at each fraction D\ :sub:`2`\ O

    '''

#    composite_scattering_intensity = numpy.empty((3,number_of_data_points,2))
    composite_scattering_intensity = numpy.zeros((3, number_of_data_points, 2))
    calculated_scattering_data = numpy.zeros(
        (number_of_data_points, number_of_contrast_points))
    mean_square_difference = []
    reduced_chi_squared_list = []
#    print('composite_scattering_intensity: ', composite_scattering_intensity)

# big loop over each q
    for q in range(number_of_data_points):
        X = numpy.zeros((3, 3))
        Xi = numpy.zeros((3, 3))
        Y = numpy.zeros((3, 1))
        I = numpy.zeros((3, 1))
        delta_rho_xy = numpy.zeros(3)
        mean_square_difference.append([])

# for a given q value, build the X and Y matrices at each contrast
# NOTE:  these delta_rho variables seem to be calculated again in calculate reduced chi squared
        for i in range(number_of_contrast_points):
            delta_rho_xy[0] = delta_rho[i][0]*delta_rho[i][0]
            delta_rho_xy[1] = delta_rho[i][0]*delta_rho[i][1]
            delta_rho_xy[2] = delta_rho[i][1]*delta_rho[i][1]
#            intensity = scale_factor[i]*scattering_data[i][q][1]
#            intensity_error = scale_factor[i]*scattering_data[i][q][2]
            intensity = rescaled_scattering_data[i][q][1]
            intensity_error = rescaled_scattering_data[i][q][2]

#  The weight is 1/errI(q)**2
            if intensity_error > 0:
                w = 1/intensity_error
# Set up the multiple linear regression matrix per Bevington (chapter 7):  X*I = Y
# Y is a column matrix of order 3 (since we are solving for 3 unknowns)
# Y[k] = sum over k[(1/IexpErr**2)*Iexp*delta_rho_xy[k]]  (eqs. 7.3 and 7.14)
# I is a column matrix of order 3 that contains the solution (I11, I12 and I22)
# X is a 3x3 matrix containing the delta_rho_xy parameters
# X[k][l]  = sum over k and l[(1/IexpErr**2)*delta_rho_xy[k]*delta_rho_xy[l]] (eq. 7.15)
# invert the X matrix to obtain Xi (the covariance matrix); then Xi*Y = I

                for k in range(3):
                    Y[k] += w*w*delta_rho_xy[k]*intensity
                    for l in range(3):
                        X[k][l] += w*w*delta_rho_xy[k]*delta_rho_xy[l]

# using numpy arrays rather than numpy.matrix, as the latter may be removed according to the numpy manual
        Xi = numpy.linalg.inv(X)
        I = numpy.dot(Xi, Y)
#        Xi = numpy.linalg.inv(numpy.matrix(X))
#        print('covariance matrix: ', Xi)
#        I = Xi*Y
#        print('solution: ', I)

        for j in range(3):
            composite_scattering_intensity[j][q][0] = I[j]

        msd_list, reduced_chi_squared, calculated_intensity = calculate_reduced_chi_squared(
            delta_rho, rescaled_scattering_data, composite_scattering_intensity, number_of_contrast_points, q)

        for m in range(number_of_contrast_points):
            mean_square_difference[q].append(msd_list[m])
            calculated_scattering_data[q][m] = calculated_intensity[m]
        reduced_chi_squared_list.append(reduced_chi_squared)

# error is the square root of the corresponding diagonal element in the covariance matrix
        for k in range(3):
            composite_scattering_intensity[k][q][1] = numpy.sqrt(
                Xi[k, k])
#            composite_scattering_intensity[k][q][1] = numpy.sqrt(reduced_chi_squared*Xi[k, k])
#        print('q, model scattering data: ', q, composite_scattering_intensity)

#    print('mean_square_difference, composite_scattering_intensity, reduced_chi_squared_list: ', mean_square_difference, composite_scattering_intensity, reduced_chi_squared_list)

#    print('model_scattering data final: ', composite_scattering_intensity)

    return mean_square_difference, composite_scattering_intensity, calculated_scattering_data, reduced_chi_squared_list
