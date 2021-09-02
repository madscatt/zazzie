# -*- coding: utf-8 -*-
'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import locale
import string
import sassie.interface.input_filter_sasmol as input_filter



#check the variables; the variables checked will depend on the method used
#don't worry about following SASSIE 2.0 format using "self"
#this version doesn't check the multi_component_analysis_variables to see if they are floating numbers
#need to rewrite the input filter to allow for a nested set of lists such as delta_rho


def check_multi_component_analysis(variables, multi_component_analysis_variables, **kwargs):

    run_name = variables['run_name'][0]
    output_file_name = variables['output_file_name'][0]
    input_file_name = variables['input_file_name'][0]
    read_from_file = variables['read_from_file'][0]
    number_of_contrast_points = variables['number_of_contrast_points'][0]
    number_of_components = variables['number_of_components'][0]
    stoichiometry_flag = variables['stoichiometry_flag'][0]
    match_point_flag = variables['match_point_flag'][0] 
    stuhrmann_parallel_axis_flag = variables['stuhrmann_parallel_axis_flag'][0]
    decomposition_flag = variables['decomposition_flag'][0]

#check run_name
    error = []
    error = input_filter.check_name(run_name)
    if(error != []):
        return error

#check read_from_file (Will an error be raised prior to this point?  YES! The program will crash earlier on, so this test isn't needed.)
#is the first line below better?
#   if(read_from_file not True and read_from_file not False):
    if(read_from_file != True and read_from_file != False):
        error.append('read_from_file must be True or False')
        return error

#check that input file exists
    if(read_from_file == True):
        error = input_filter.check_file_exists(input_file_name)
        if(len(error) > 0):
            error.append('input file is not readable or does not exist')
            return error

#check if number of contrast points is >= the number of components
    if(number_of_contrast_points < number_of_components):
        error.append('number of contrasts must be >= number of components')
        return error
        
#check if at least one method has been chosen. NOTE: This check isn't needed since multi_component_analysis_variables won't be defined if one of these flags aren't true and there will be an error before we get to this point.  How do we capture this prior error so the program doesn't just crash? Perhaps make multi_component_analysis_variables blank to start?
    if(not stoichiometry_flag and not match_point_flag and not stuhrmann_parallel_axis_flag and not decomposition_flag):
        error.append('at least one method must be selected')
        return error
        
#check the multi_component analysis variables depending on the method (will become if, elif). 

    if(stoichiometry_flag == True):
        fraction_d2o = multi_component_analysis_variables[0]
        partial_specific_volume = multi_component_analysis_variables[1]
        izero = multi_component_analysis_variables[2]
        concentration = multi_component_analysis_variables[3]
        delta_rho = multi_component_analysis_variables[4]

#check if length of fraction_d2o, izero, concentration and delta_rho = number of contrasts
        if(len(fraction_d2o) != number_of_contrast_points):
            error.append('fraction_d2o must have %i values' %(number_of_contrast_points))
        if(len(izero) != number_of_contrast_points):
            error.append('izero must have %i values' %(number_of_contrast_points))
        if(len(concentration) != number_of_contrast_points):
            error.append('concentration must have %i values' %(number_of_contrast_points))
        if(len(delta_rho) != number_of_contrast_points):
            error.append('delta_rho must have %i sets of values' %(number_of_contrast_points))

#check if length delta_rho[i] = number of components. This check is only needed if read_from_file is False since delta_rho will be initialized properly if it is being read from a contrast output file.
        if(read_from_file == False):
            for i in range(number_of_contrast_points):
                if(len(delta_rho[i]) != number_of_components):
                    error.append('delta rho[%i] must have %i values' %(i,number_of_components))
                    return error

#check if length of partial specific volume = number of components

        if(len(partial_specific_volume) != number_of_components):
            error.append('partial_specific_volume must have %i values' %(number_of_components))
            return error
            
    return error


