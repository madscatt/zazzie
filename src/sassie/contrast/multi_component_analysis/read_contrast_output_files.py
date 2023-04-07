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
#       MULTI-COMPONENT ANALYSIS
#
#       08/12/2021       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Read Contrast Output Files** contains the methods that reads the output file from
    **Contrast Calculator** to obtain contrast (:math:`\\Delta \\rho`), Mw, I(0) or other information for
    further calculation. It can be used for both experiment planning and data analysis.

    **Inputs:**
    
        input file name, number of contrast points, number of components in the complex
        fraction of D\ :sub:`2`\ O in the solvent for each contrast point
    
    **Outputs:**
    
        contrast, Mw, I(0) for each component at each contrast read from **Contrast Calculator** output files
            
    Called by the **Multi-component Analysis** module.
    
'''

# Right now, this file only contains a method that reads contrast values from the *_contrast.txt output files in order to obtain contrast values at a particular fraction D2O.

import io
import os
import sys
import string
import locale
import numpy
import scipy.stats


def read_contrast_file(other_self):
    '''
    **Read Contrast File** reads the contrast (:math:`\\Delta \\rho`) values for each component from the **Contrast Calculator** output file and returns the contrasts at the fraction D\ :sub:`2`\ O values chosen by the user.    

    Parameters
    ----------
    input_file_name: string
        The name of the contrast calculator output file that contains the desired information
    number_of_contrast_points: int
        The number of solvent conditions with different fraction D\ :sub:`2`\ O values
    number_of_components: int
        The number of components in the molecule with different scattering length densities
    fraction_d2o:   float array (dimension = number_of_contrast_points)
        The fraction D\ :sub:`2`\ O values that define the contrasts

    Returns
    -------
    delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
        The contrast for each component at all fraction D\ :sub:`2`\ O values of interest

    Note
    ----
    - At this time, a table of contrast values is read from the \\*_contrast.txt contrast calculator output file.  Since the contrast values are linear with fraction D\ :sub:`2`\ O in the solvent, contrast (`y`) vs fraction D\ :sub:`2`\ O (`x`) values are fit to a straight line using **scipy.stats.linregress**, which performs an **unweighted fit**.  The resultant equation is used to find the contrast values (`delta_rho`) at the experimental fraction D\ :sub:`2`\ O values that are input by the user. 

    - In the future, additional return values are expected to be implemented:

        **izero:**                      float array (dimension = number_of_contrast points)

        **partial_specific_volume:**    float array (dimension = number_of_components)

        **molecular_weight:**           float array (dimension = number_of_componenets) 

    - While the same code can be used to call **scipy.stats.linregress** when converting Python 2 (scipy version 1.2.1) to Python 3 (scipy version 1.7), the read_contrast_file method can now use the newer Python 3 syntax if desired. The Python 2 syntax still works.  The main difference is that the newer syntax allows the intercept standard error to be captured.  As of now, neither the slope standard error nor the intercept standard error is used further.  

        **Python 2:**

        slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(x,y)

        **Python 3:**

        result = scipy.stats.linregress(x,y)

        slope = result.slope

        intercept = result.intercept

        r_value = result.rvalue

        p_value = result.pvalue

        slope_std_error = result.stderr

        intercept_std_error = result.intercept_stderr  

    '''

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    mvars = other_self.module_variables

    log.debug('in read contrast file')

    with io.open(mvars.input_file_name, 'r') as infile:
        lines = infile.readlines()
    # infile.close()

    # Read the fraction_d2o values of interest and associate them with numpy array x (fraction D2O).
    # Convert to numpy array
    x = numpy.zeros(shape=(mvars.number_of_contrast_points))
    for i in range(mvars.number_of_contrast_points):
        x[i] = mvars.fraction_d2o[i]
#    print ('frac D2O: ', x)

    number_of_lines = len(lines)
    # print ("number of lines ", number_of_lines)

# Find linear equation to calculate desired delta_rho values. This equation can eventually be part of the contrast calculator output file so it can be read directly instead of having to re-calculate it.
    log.debug('\n')
    log.debug('reading contrasts from input file: %s' %
              (mvars.input_file_name))
# First, find out where the desired delta_rho values are in the file.
    for i in range(number_of_lines):
        # print(i, lines[i][0:4])
        if(lines[i][0:4] == u'# fr'):
            starting_line = i
            # print('start here: ', i)

    # print(lines[starting_line])
#    l = string.split(lines[starting_line], ',')
    l = lines[starting_line].split(',')  # for Python 3 compatibility
    # print(l)

# Associate the contrast arrays with a title [protein, RNA, DNA, Molecule 1, etc.] for debugging purposes. This may eventually be a module variable that is input from the GUI, as we may have to ask users which contrast goes first to be in agreement with the way they are inputing partial specific volume, etc.
    title = []
    title.append(u' frac D2O')
    for i in range(mvars.number_of_components+1):
        title.append(l[i+1])
    # print(title)
    # print(title[0])
    # print(title[1])
    # print(title[2])

    number_of_data_points = number_of_lines - starting_line - 1
    # print('number of data points for linear fit: ', number_of_data_points)

    data = numpy.zeros(
        shape=(mvars.number_of_components + 1, number_of_data_points))
    # print('data: ', data)

    for i in range(starting_line+1, number_of_lines):
        # print('i: ', i)
        #        l = string.split(lines[i], '\t')
        l = lines[i].split('\t')  # for Python 3 compatibility
        # print(l)
        # print(l[0])
        # print(l[1])
        # print(l[2])
        for j in range(mvars.number_of_components+1):
            # print('j,line(j): ',j, l[j])
            data[j][i-starting_line-1] = locale.atof(l[j])
            # print(data[j][i-starting_line-1])

    # print('data: ', data)
    # print('data[0]: ', data[0])
    # print('data[1]: ', data[1])
    # print('data[0][0]: ', data[0][0])


# Find the slope and intercept for each component and put them into arrays. (I'm doing it this way in case we end up just reading the values from an equation that is included in the contrast calculator output file. Right now we have to find the linear equation from the table of values we just read.)

    slope_values = numpy.zeros(shape=mvars.number_of_components)
    intercept_values = numpy.zeros(shape=mvars.number_of_components)

# TODO: This part should be a separate method such as "calculate match point" (because it is a loop rather than a single line), perhaps to used in contrast calculator also; for now, it could just be a method in this file.
    for i in range(mvars.number_of_components):
        slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(
            data[0], data[i+1])
        x_intercept = -intercept / slope
        # print ('slope, intercept, x_intercept: ', slope, intercept, x_intercept)
        # Match points are calculated as a check.
        match_point = x_intercept * 100  # calculated as a check
        log.debug('match point %s: %s' % (title[i+1], str(match_point)))
        slope_values[i] = slope
        intercept_values[i] = intercept
    # print('slope_values: ', slope_values)
    # print('intercept_values: ', intercept_values)

# Find the delta_rho values and put them into a delta_rho array.
    mvars.delta_rho = [[] for n in range(mvars.number_of_contrast_points)]
    for i in range(mvars.number_of_contrast_points):
        info = []  # for debugging; if writing to file, make 'title' a global variable to pass back to main program
        for j in range(mvars.number_of_components):
            result = slope_values[j]*x[i] + intercept_values[j]
            this_delta_rho = round(result, 3)
#            this_delta_rho = '%.3f' % (slope_values[j]*x[i] + intercept_values[j]) #this gives a unicode output.
# TODO: If these values are going to populate the GUI, they need to be converted to string.  In the GUI mimic, we use a nested_float_array in the format: '-3.2, -5.7; 1.6, 0.26; 0.031, -1.74'.
            mvars.delta_rho[i].append(this_delta_rho)
            info.append(str(title[j+1]))
        log.debug('delta_rho for fraction D2O %s %s: %s' %
                  (str(mvars.fraction_d2o[i]), info, str(mvars.delta_rho[i])))

    log.debug(vars(mvars))
    return
