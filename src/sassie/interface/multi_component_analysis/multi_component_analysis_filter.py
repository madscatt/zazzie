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
#       MULTI-COMPONENT ANALYSIS FILTER
#
#       08/30/2021       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    MULTI-COMPONENT ANALYSIS FILTER is the method that checks the inputs for
    the MULTI-COMPONENT ANALYSIS module that were not previously checked by
    INPUT FILTER, which only checks for valid string, float, integer, boolean, etc.
        
    Called from GUI_MIMIC_MULTI_COMPONENT_ANALYSIS
    Calls INPUT_FILTER
          
    INPUTS:
        run_name
        output_file_name
        input_file_name
        read_from_file
        number_of_contrast_points
        number_of_components
        stoichiometry_flag
        match_point_flag
        stuhrmann_parallel_axis_flag
        decomposition_flag
        fraction_d2o
        initial_matchpoint_guess
        izero
        izero_error
        concentration
        concentration_error
        partial_specific_volume
        molecular_weight
        radius_of_gyration
        radius_of_gyration_error
        delta_rho
            
    OUTPUTS:
        error string 
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import locale
import string
import sassie.interface.input_filter as input_filter


#this version uses the new input filter that recognizes a nested_float_array
'''
 ### pseudo code stub to filter boolean list so at least ONE and only ONE value is TRUE 
 bool IsExactlyOneBooleanTrue( bool *boolAry, int size )
    {
      bool areAnyTrue = false;
      bool areTwoTrue = false;
      for(int i = 0; (!areTwoTrue) && (i < size); i++) {
        areTwoTrue = (areAnyTrue && boolAry[i]);
        areAnyTrue |= boolAry[i];
      }
      return ((areAnyTrue) && (!areTwoTrue));
    }

https://stackoverflow.com/questions/14888174/how-do-i-determine-if-exactly-one-boolean-is-true-without-type-conversion

'''

def check_multi_component_analysis(variables, **kwargs):

    run_name = variables['run_name'][0]
    output_file_name = variables['output_file_name'][0]
    input_file_name = variables['input_file_name'][0]
    read_from_file = variables['read_from_file'][0]
    number_of_contrast_points = variables['number_of_contrast_points'][0]
    fraction_d2o = variables['fraction_d2o'][0]
    stoichiometry_flag = variables['stoichiometry_flag'][0]
    match_point_flag = variables['match_point_flag'][0] 
    stuhrmann_parallel_axis_flag = variables['stuhrmann_parallel_axis_flag'][0]
    decomposition_flag = variables['decomposition_flag'][0]

#check run_name
    error = []
    error = input_filter.check_name(run_name)
    if(error != []):
        return error

#check read_from_file (Will an error be raised prior to this point?  YES! The program will crash in the input filter, so this test isn't needed. 
    if(read_from_file != True and read_from_file != False):
        error.append('read_from_file must be True or False')
        return error

#check that input file exists
    if(read_from_file == True):
        error = input_filter.check_file_exists(input_file_name)
        if(len(error) > 0):
            error.append('input file is not readable or does not exist')
            return error

#check that at least one method is chosen.
    if(not stoichiometry_flag and not match_point_flag and not stuhrmann_parallel_axis_flag and not decomposition_flag):
        error.append('at least one method must be selected')
        return error

#TODO:  Do we need to check if only one flag is True? Will more than one flag be allowed to be True at a time to execute methods sequentially? Right now, this is not the way the program is envisioned.  Like MulCh, you would choose a method and execute it.  Then choose another an execute, etc.  If the method is selected by a drop down menu, then we don't need to check to make sure that at least one flag is chosen and that only one flag is True.
        
#check the multi_component analysis variables depending on the method.
#    print('len frac d2o: ', len(fraction_d2o))
#variables common to all methods:  fraction_d2o
#check if length of fraction D2O = number of contrasts
    if(len(fraction_d2o) != number_of_contrast_points):
        error.append('fraction D2O must have %i values' %(number_of_contrast_points))
        return error
#check that fraction D2O is between 0 and 1
    for i in range(number_of_contrast_points):
        if(fraction_d2o[i] < 0 or fraction_d2o[i] > 1):
            error.append('fraction D2O[%i] must be between 0 and 1' %(i))
            return error

    if(stoichiometry_flag == True):
        number_of_components = variables['number_of_components'][0]
        izero = variables['izero'][0]
        concentration = variables['concentration'][0]
        partial_specific_volume = variables['partial_specific_volume'][0]
        delta_rho = variables['delta_rho'][0]

#check if number of contrast points is >= the number of components
        if(number_of_contrast_points < number_of_components):
            error.append('number of contrasts must be >= number of components')
            return error
#check if length of izero, concentration and delta_rho = number of contrasts
        if(len(izero) != number_of_contrast_points):
            error.append('I(0) must have %i values' %(number_of_contrast_points))
            return error
        if(len(concentration) != number_of_contrast_points):
            error.append('concentration must have %i values' %(number_of_contrast_points))
            return error
#If read_from_file is True, the values are read from a contrast calculator output file, which is done in the unpack variables method in the main program. The way the program is written now, this executes after these checks are performed.  If we execute the reading of the values from the file to fill the value into the GUI before executing the main program, then we can perform this whether read_from_file is True or False and the if statement can be removed.
        if(read_from_file == False):
            if(len(delta_rho) != number_of_contrast_points):
                error.append('delta rho must have %i sets of values' %(number_of_contrast_points))
                return error
#check if length delta_rho[i] = number of components 
            for i in range(number_of_contrast_points):
                if(len(delta_rho[i]) != number_of_components):
                    error.append('delta rho[%i] must have %i values' %(i,number_of_components))
                    return error
#check if length of partial specific volume = number of components
        if(len(partial_specific_volume) != number_of_components):
            error.append('partial_specific_volume must have %i values' %(number_of_components))
            return error


    elif(match_point_flag == True):
        izero = variables['izero'][0]
        izero_error = variables['izero_error'][0]
        concentration = variables['concentration'][0]
        concentration_error = variables['concentration_error'][0]
        initial_match_point_guess = variables['initial_match_point_guess'][0]
#        print('init matchpoint guess: ', initial_match_point_guess)
        
#check if length of izero, izero_error, concentration and concentration_error = number of contrasts
        if(len(izero) != number_of_contrast_points):
            error.append('I(0) must have %i values' %(number_of_contrast_points))
            return error
        if(len(izero_error) != number_of_contrast_points):
            error.append('I(0) error must have %i values' %(number_of_contrast_points))
            return error
        if(len(concentration) != number_of_contrast_points):
            error.append('concentration must have %i values' %(number_of_contrast_points))
            return error
        if(len(concentration_error) != number_of_contrast_points):
            error.append('concentration error must have %i values' %(number_of_contrast_points))
            return error
#check if izero_error and concentration_error are non-zero since a weighted fit is performed
        for i in range(number_of_contrast_points):
                if(izero_error[i] == 0.0):
                    error.append('I(0) error[%i] cannot equal zero' %(i))
                    return error
                if(concentration_error[i] == 0.0):
                    error.append('concentration error[%i] cannot equal zero' %(i))
                    return error
#check that initial match point guess is between 0 and 1
        if(initial_match_point_guess < 0 or initial_match_point_guess > 1):
            error.append('initial match point guess must be between 0 and 1')
            return error


    elif(stuhrmann_parallel_axis_flag == True):
        number_of_components = variables['number_of_components'][0]
        molecular_weight = variables['molecular_weight'][0]
        partial_specific_volume = variables['partial_specific_volume'][0]
        radius_of_gyration = variables['radius_of_gyration'][0]
        radius_of_gyration_error = variables['radius_of_gyration_error'][0]
        delta_rho = variables['delta_rho'][0]

#check if number of contrast points is >= the number of components
        if(number_of_contrast_points < number_of_components):
            error.append('number of contrasts must be >= number of components')
            return error
#delta_rho_checks if the values are input by hand
        if(read_from_file == False):
            if(len(delta_rho) != number_of_contrast_points):
                error.append('delta rho must have %i sets of values' %(number_of_contrast_points))
                return error
#check if length delta_rho[i] = number of components. 
            for i in range(number_of_contrast_points):
                if(len(delta_rho[i]) != number_of_components):
                    error.append('delta rho[%i] must have %i values' %(i,number_of_components))
                    return error
#check if length of partial specific volume = number of components
        if(len(partial_specific_volume) != number_of_components):
            error.append('partial_specific_volume must have %i values' %(number_of_components))
            return error
#check if length of molecular weight = number of components
        if(len(molecular_weight) != number_of_components):
            error.append('Mw must have %i values' %(number_of_components))
            return error
#check if length of izero, concentration and delta_rho = number of contrasts
        if(len(radius_of_gyration) != number_of_contrast_points):
            error.append('Rg must have %i values' %(number_of_contrast_points))
            return error
        if(len(radius_of_gyration_error) != number_of_contrast_points):
            error.append('Rg error must have %i values' %(number_of_contrast_points))
            return error
#check if radius_of_gyration_error is non-zero since a weighted fit is performed
        for i in range(number_of_contrast_points):
                if(radius_of_gyration_error[i] == 0.0):
                    error.append('Rg error[%i] cannot equal zero' %(i))
                    return error

    return error


