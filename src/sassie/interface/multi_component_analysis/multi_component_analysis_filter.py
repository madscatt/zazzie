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
#       04/03/2023       --      python 3 coding        :   Joseph E. Curtis
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
"""
    **Multi-component Analysis Filter** is the method that checks the inputs for
    the **Multi-component Analysis** module that were not previously checked by
    **Input Filter**, which only checks for valid string, float, integer, boolean, etc.

    **Inputs:**
    
        The inputs checked depend on which method is used.
     
        The method is chosen based on which of the following flags is True:
    
        - match_point_flag: Match Point Analysis
        - stuhrmann_parallel_axis_flag: Sturhmann and Parallel Axis Theorem Analyses
        - decomposition_flag: Decomposition Analysis
        - stoichiometry_flag: Stoichiometry Analysis
        
        Only one flag at a time can be True.
    
    **Outputs:**
    
        Error string

    Called from **Gui Mimic Multi-component Analysis**
    
    Calls **Input Filter**

"""

import sassie.interface.input_filter as input_filter

# this version uses the new input filter that recognizes a nested_float_array


def check_multi_component_analysis(variables, **kwargs):
    """
    Method to check the **Multi-component Analysis** variables.
        
    Parameters
    ----------
    
        run_name: string
            run name
        output_file_name: string
            user-specified output file name
        stoichiometry_flag: boolean
            flag to determine if stoichiometry analysis is being used
        match_point_flag: boolean
            flag to determine if match point analysis is being used 
        stuhrmann_parallel_axis_flag: boolean
            flag to determine if Stuhrmann and Parallel Axis methods are being used
        decomposition_flag: boolean
            flag to determine if decomposition analysis is being used
        number_of_contrast_points:  int
            The number of solvent conditions with different fraction D\ :sub:`2`\ O values
        fraction_d2o:   float array (dimension = number_of_contrast_points)
            The fraction D\ :sub:`2`\ O values that define the contrasts
        izero:  float array (dimension = number_of_contrast_points)
            I(0) value at each contrast in cm\ :sup:`-1`\ 
        izero_error:  float array (dimension = number_of_contrast_points)
            I(0) error value at each contrast
        concentration:  float array (dimension = number_of_contrast_points)
            concentration at each contrast in mg/mL
        concentration_error:  float array (dimension = number_of_contrast_points)
            concentration error at each contrast
        initial_match_point_guess:  float
            The fraction D\ :sub:`2`\ O value to be used as initial match point guess    
        partial_specific_volume: float array (dimension = number of components)
            partial specific volume of each component
        molecular_weight: float array (dimension = number of components)
            molecular_weight of each component
        radius_of_gyration: float array (dimension = number_of_contrast_points)
            radius of gyration at each contrast in Angstroms
        radius_of_gyration_error: float array (dimension = number_of_contrast_points)
            radius of gyration error at each contrast in Angstroms
        delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
            The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )        

    Returns
    -------
    
        error: string
            The error message generated when a check fails. If there are no failures, the error is blank.
    
    """
    #### INPUT VARIABLES FOR ALL METHODS

    match_point_flag = variables["match_point_flag"][0]
    stuhrmann_parallel_axis_flag = variables["stuhrmann_parallel_axis_flag"][0]
    stoichiometry_flag = variables["stoichiometry_flag"][0]
    decomposition_flag = variables["decomposition_flag"][0]

    run_name = variables["run_name"][0]
    output_file_name = variables["output_file_name"][0]
    #input_file_name = variables["input_file_name"][0]
    read_from_contrast_calculator_output_file = variables[
        "read_from_contrast_calculator_output_file"
    ][0]

    number_of_contrast_points = variables["number_of_contrast_points"][0]
    fraction_d2o = variables["fraction_d2o"][0]

    # define empty error list to return

    error = []

    # check flags to make sure they are boolean type

    if type(match_point_flag) is not bool:
        error.append("match_point_flag is not a boolean")
        return error
    elif type(stuhrmann_parallel_axis_flag) is not bool:
        error.append("stuhrmann_parallel_axis_flag is not a boolean")
        return error
    elif type(stoichiometry_flag) is not bool:
        error.append("stoichiometry_flag is not a boolean")
        return error
    elif type(decomposition_flag) is not bool:
        error.append("decomposition_flag is not a boolean")
        return error

    flag_list = [
        match_point_flag,
        stuhrmann_parallel_axis_flag,
        stoichiometry_flag,
        decomposition_flag,
    ]

    # count the number of True values in the list
    number_of_true_values = flag_list.count(True)

    if number_of_true_values != 1:
        error.append(
            "one method value must be True and only one method flag can be True\n"
        )
        error.append("match_point_flag = " + str(match_point_flag) + "\n")
        error.append(
            "stuhrmann_parallel_axis_flag = " + str(stuhrmann_parallel_axis_flag) + "\n"
        )
        error.append("stoichiometry_flag = " + str(stoichiometry_flag) + "\n")
        error.append("decomposition_flag = " + str(decomposition_flag) + "\n")
        return error

    # check run_name

    error = input_filter.check_name(run_name)
    if error != []:
        return error

    # check read_from_contrast_calculator_output_file (Will an error be raised prior to this point?  YES! The program will crash in the input filter, so this test isn't needed.
    if type(read_from_contrast_calculator_output_file) is not bool:
        error.append("read_from_contrast_calculator_output_file must be True or False")
        return error

    # check that input file exists
    if read_from_contrast_calculator_output_file:
        error = input_filter.check_file_exists(input_file_name)
        if len(error) > 0:
            error.append(
                "read_from_contrast_calculator_output_file file is not readable or does not exist"
            )
            return error

    if len(fraction_d2o) != number_of_contrast_points:
        error.append("fraction D2O must have %i values" % (number_of_contrast_points))
        return error
    # check that fraction D2O is between 0 and 1
    for i in range(number_of_contrast_points):
        if fraction_d2o[i] < 0 or fraction_d2o[i] > 1:
            error.append("fraction D2O[%i] must be between 0 and 1" % (i))
            return error

    if match_point_flag:
        izero = variables["izero"][0]
        izero_error = variables["izero_error"][0]
        concentration = variables["concentration"][0]
        concentration_error = variables["concentration_error"][0]
        initial_match_point_guess = variables["initial_match_point_guess"][0]

        # check if length of izero, izero_error, concentration and concentration_error = number of contrasts
        if len(izero) != number_of_contrast_points:
            error.append("I(0) must have %i values" % (number_of_contrast_points))
            return error
        if len(izero_error) != number_of_contrast_points:
            error.append("I(0) error must have %i values" % (number_of_contrast_points))
            return error
        if len(concentration) != number_of_contrast_points:
            error.append(
                "concentration must have %i values" % (number_of_contrast_points)
            )
            return error
        if len(concentration_error) != number_of_contrast_points:
            error.append(
                "concentration error must have %i values" % (number_of_contrast_points)
            )
            return error
        # check if izero_error and concentration_error are non-zero since a weighted fit is performed
        for i in range(number_of_contrast_points):
            if izero_error[i] == 0.0:
                error.append("I(0) error[%i] cannot equal zero" % (i))
                return error
            if concentration_error[i] == 0.0:
                error.append("concentration error[%i] cannot equal zero" % (i))
                return error
        # check that initial match point guess is between 0 and 1
        if initial_match_point_guess < 0 or initial_match_point_guess > 1:
            error.append("initial match point guess must be between 0 and 1")
            return error

    elif stuhrmann_parallel_axis_flag:
        number_of_components = variables["number_of_components"][0]
        molecular_weight = variables["molecular_weight"][0]
        partial_specific_volume = variables["partial_specific_volume"][0]
        radius_of_gyration = variables["radius_of_gyration"][0]
        radius_of_gyration_error = variables["radius_of_gyration_error"][0]
        delta_rho = variables["delta_rho"][0]

        # check if number of contrast points is >= the number of components
        if number_of_contrast_points < number_of_components:
            error.append("number of contrasts must be >= number of components")
            return error
        # delta_rho_checks if the values are input by hand

        #### TODO: read_from_sascalc_output_file = False
        #### TODO: read_from_contrast_calculator_output_file:

        read_from_contrast_calculator_output_file = variables["read_from_contrast_calculator_output_file"][0]
        read_from_sascalc_output_file = variables["read_from_sascalc_output_file"][0]

        if not read_from_contrast_calculator_output_file:
            if len(delta_rho) != number_of_contrast_points:
                error.append(
                    "delta rho must have %i sets of values"
                    % (number_of_contrast_points)
                )
                return error
            # check if length delta_rho[i] = number of components.
            for i in range(number_of_contrast_points):
                if len(delta_rho[i]) != number_of_components:
                    error.append(
                        "delta rho[%i] must have %i values" % (i, number_of_components)
                    )
                    return error
        # check if length of partial specific volume = number of components
        if len(partial_specific_volume) != number_of_components:
            error.append(
                "partial_specific_volume must have %i values" % (number_of_components)
            )
            return error
        # check if length of molecular weight = number of components
        if len(molecular_weight) != number_of_components:
            error.append("Mw must have %i values" % (number_of_components))
            return error
        # check if length of izero, concentration and delta_rho = number of contrasts
        if len(radius_of_gyration) != number_of_contrast_points:
            error.append("Rg must have %i values" % (number_of_contrast_points))
            return error
        if len(radius_of_gyration_error) != number_of_contrast_points:
            error.append("Rg error must have %i values" % (number_of_contrast_points))
            return error
        # check if radius_of_gyration_error is non-zero since a weighted fit is performed
        for i in range(number_of_contrast_points):
            if radius_of_gyration_error[i] == 0.0:
                error.append("Rg error[%i] cannot equal zero" % (i))
                return error

    elif stoichiometry_flag:
        number_of_components = variables["number_of_components"][0]
        izero = variables["izero"][0]
        concentration = variables["concentration"][0]
        partial_specific_volume = variables["partial_specific_volume"][0]
        delta_rho = variables["delta_rho"][0]

        # check if number of contrast points is >= the number of components
        if number_of_contrast_points < number_of_components:
            error.append("number of contrasts must be >= number of components")
            return error
        # check if length of izero, concentration and delta_rho = number of contrasts
        if len(izero) != number_of_contrast_points:
            error.append("I(0) must have %i values" % (number_of_contrast_points))
            return error
        if len(concentration) != number_of_contrast_points:
            error.append(
                "concentration must have %i values" % (number_of_contrast_points)
            )
            return error
        # If read_from_file is True, the values are read from a contrast calculator output file, which is done in the unpack variables method in the main program. The way the program is written now, this executes after these checks are performed.  If we execute the reading of the values from the file to fill the value into the GUI before executing the main program, then we can perform this whether read_from_file is True or False and the if statement can be removed.
        read_from_contrast_calculator_output_file = variables["read_from_contrast_calculator_output_file"][0]
        if not read_from_contrast_calculator_output_file:
            if len(delta_rho) != number_of_contrast_points:
                error.append(
                    "delta rho must have %i sets of values"
                    % (number_of_contrast_points)
                )
                return error
            # check if length delta_rho[i] = number of components
            for i in range(number_of_contrast_points):
                if len(delta_rho[i]) != number_of_components:
                    error.append(
                        "delta rho[%i] must have %i values" % (i, number_of_components)
                    )
                    return error
        # check if length of partial specific volume = number of components
        if len(partial_specific_volume) != number_of_components:
            error.append(
                "partial_specific_volume must have %i values" % (number_of_components)
            )
            return error

    return error


if __name__ == "__main__":
    variables = {}

    variables["match_point_flag"] = (True, "boolean")
    variables["match_point_flag"] = (False, "boolean")
    #variables["stuhrmann_parallel_axis_flag"] = (True, "boolean")
    variables["stuhrmann_parallel_axis_flag"] = (False, "boolean")
    variables["stoichiometry_flag"] = (True, "boolean")
    #variables["stoichiometry_flag"] = (False, "boolean")
    variables["decomposition_flag"] = (False, "boolean")

    variables["run_name"] = ("run_0", "string")
    variables["output_file_name"] = "test.out" "string"
    variables["input_file_name"] = ("input_test.txt", "string")
    variables["read_from_contrast_calculator_output_file"] = (False, "boolean")

    if variables["match_point_flag"][0]:

        variables["number_of_contrast_points"] = (7, "int")
        variables["fraction_d2o"] = ([1.0, 0.9, 0.8, 0.4, 0.2, 0.1, 0.0], "int_array")

        variables["izero"] = ([0.537, 0.332, 0.19, 0.0745, 0.223, 0.352, 0.541], 'float_array')
        variables["izero_error"] = ([0.001, 0.002, 0.001, 0.002, 0.002, 0.002, 0.003], 'float_array')
        variables["concentration"] = ([11.9, 11.9, 11.9, 26.9, 11.9, 11.9, 11.9], 'float_array')
        variables["concentration_error"] = ([0.6, 0.6, 0.6, 1.3, 0.6, 0.6, 0.6], 'float_array')
        variables["initial_match_point_guess"] = (0.5, 'float')

    elif variables["stuhrmann_parallel_axis_flag"][0]:

        variables["read_from_sascalc_output_file"] = (False, 'boolean')

        variables["number_of_contrast_points"] = (7, "int")
        variables["fraction_d2o"] = ([1.0, 0.9, 0.8, 0.4, 0.2, 0.1, 0.0], "int_array")

        variables["number_of_components"] = (2, 'int')
        variables["partial_specific_volume"] = ([0.73, 0.73], 'float_array')
        variables["molecular_weight"] = ([50.7, 11.7], 'float_array')
        variables['delta_rho'] = ([[-3.34, 0.41], [-2.78, 0.98], [-2.21, 1.54], [0.055, 3.8], [1.18, 4.93], [1.75, 5.49], [2.31, 6.06]], 'nested_float_array')
        variables['radius_of_gyration'] = ([25.11, 24.16, 23.03, 23.4, 28.22, 28.33, 28.85], 'float_array')
        variables['radius_of_gyration_error'] = ([0.09, 0.14, 0.2, 0.7, 0.29, 0.19, 0.12], 'float_array')

    elif variables["stoichiometry_flag"][0]:

        variables['number_of_contrast_points'] = (3, 'int')
        variables['fraction_d2o'] = ([0.99, 0.12, 0.41], 'float_array')
        variables['concentration'] = ([3.7, 3.6, 3.1], 'float_array')
        variables['concentration_error'] = ([0.18, 0.18, 0.18], 'float_array')
        variables['number_of_components'] = (2, 'int')
        variables['partial_specific_volume'] = ([0.745, 0.903], 'float_array')
        variables['delta_rho'] = ([[-3.2, -5.7], [1.6, 0.26], [0.031, -1.74]], 'nested_float_array')
        variables['izero'] = ([8.4, 0.6, 0.17], 'float_array')
        variables['izero_error'] = ([0.2, 0.04, 0.01], 'float_array')

    elif variables["decomposition_flag"][0]:
        pass

    error = check_multi_component_analysis(variables)
    if (len(error)) == 0:
        print("NO ERRORS FOUND")
    else:
        print("len(__main__ error) = " + str(len(error)))
        from pprint import pprint

        pprint(error)
