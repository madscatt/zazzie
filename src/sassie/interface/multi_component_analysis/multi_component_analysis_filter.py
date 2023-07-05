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
#       08/30/2021       --      initial coding               :  Susan Krueger
#       04/03/2023       --      python 3 coding              :  Joseph E. Curtis
#       05/23/2023       --      added decomposition variables:  Susan Krueger
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

    Calls **Input Filter** and **Read Data File**

"""

import numpy
import sassie.interface.input_filter as input_filter
import sassie.contrast.multi_component_analysis.read_data_file as read_data_file
#import input_filter_new as input_filter
#import read_data_file as read_data_file

# this version uses the new input filter that recognizes a nested_float_array and a string array


def check_data_file(data_file_name):
    """
    Method to check that the number of data points and q binning are the same for all data files in an array.

    """

#    print("in check data file")
#    print("data file name: ", data_file_name)
    error = []
    scattering_data = []
    number_of_data_lines = []
    number_of_files = len(data_file_name)
#    print("number of files: ", number_of_files)

    for item in data_file_name:
        #        print("item: ", item)
        #        print("data file: ", item)
        q = []
        i = []
        ierr = []

        numpts, q, i, ierr = read_data_file.read_file(item)
#        print("numpts: ", numpts)
#        print("q: ", q)
#        print("i: ", i)
#        print("ierr: ", ierr)
        data = numpy.zeros((numpts, 3))
        for j in range(numpts):
            data[j][0] = q[j]
            data[j][1] = i[j]
            data[j][2] = ierr[j]
#        print(data)
        number_of_data_lines.append(numpts)
        scattering_data.append(data)

#    print("number of data lines: ", number_of_data_lines)
#    print("scattering_data1: ", scattering_data)

# checks that the number of data points is the same for all data sets
    for i in range(1, number_of_files):
        if number_of_data_lines[0] != number_of_data_lines[i]:
            error.append(
                "The number of data points, " + str(number_of_data_lines[i]) + ", in " + data_file_name[i] + " differs from " + str(
                    number_of_data_lines[0]) + " in the reference data file, " + data_file_name[0]
            )
            return error

# this checks that binning is identical for all data sets
    for i in range(number_of_files):
        for j in range(number_of_data_lines[0]):
            if scattering_data[0][j][0] != scattering_data[i][j][0]:
                error.append(
                    "The binning for point " +
                    str(j+1) + " in dataset " + str(i+1) +
                    " differs from the reference"
                )
                return error

#    print("error: ", error)
    return error


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
        read_from_contrast_calculator_output_file: boolean
            flag to determine if the contrasts (:math:`\Delta \rho`) values are read from a contrast calculator output file
        contrast_calculator_output_file_name: string
            user-specified contrast calculator output file name
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
        data_file_name: string array (dimension = number_of_contrast_points)
            contrast variation data file name at each fraction D\ :sub:`2`\ O 
        q_rg_limit_guinier: float array (dimension = number_of_contrast_points)
            qR\ :sub:`g`\  limit for the Guinier analysis at each fraction D\ :sub:`2`\ O
        starting_data_point_guinier: int array (dimension = number_of_contrast_points)
            index of the starting data point for the Guinier fit at each fraction D\ :sub:`2`\ O  (index of the first data point = 1)
        initial_points_to_use_guinier: int array (dimension = number_of_contrast_points)
            number of data points to use initially for the Guinier fit at each fraction D\ :sub:`2`\ O  (the final number of points used depends on the qR\ :sub:`g`\  limit)
        refine_scale_factor_flag: boolean
            Indicates whether the scale factor at each fraction D\ :sub:`2`\ O  will be adjusted based on the I(0) values 
        delta_rho:  2D float array (dimensions = number_of_contrast_points x number_of_components)
            The contrast for each component at all fraction D\ :sub:`2`\ O values of interest in 10\ :sup:`10`\ cm\ :sup:`-2`\  (10 :sup:`-6`\ A\ :sup:`-2`\ )

    Returns
    -------

        error: string
            The error message generated when a check fails. If there are no failures, the error is blank.

    """
    # INPUT VARIABLES FOR ALL METHODS

    match_point_flag = variables["match_point_flag"][0]
    stuhrmann_parallel_axis_flag = variables["stuhrmann_parallel_axis_flag"][0]
    stoichiometry_flag = variables["stoichiometry_flag"][0]
    decomposition_flag = variables["decomposition_flag"][0]

    run_name = variables["run_name"][0]
    output_file_name = variables["output_file_name"][0]
    read_from_contrast_calculator_output_file = variables[
        "read_from_contrast_calculator_output_file"
    ][0]
    if read_from_contrast_calculator_output_file:
        contrast_calculator_output_file_name = variables[
            "contrast_calculator_output_file_name"
        ][0]

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

    # check to make sure at least one -- but only one -- flag is True
    #print("flag list: ", flag_list)
    # count the number of True values in the list
    number_of_true_values = flag_list.count(True)
    # print("number of True values: ", number_of_true_values)

    if number_of_true_values != 1:
        error.append(
            "one method value must be True and only one method flag can be True"
        )
        error.append("match_point_flag = " + str(match_point_flag))
        error.append(
            "stuhrmann_parallel_axis_flag = " +
            str(stuhrmann_parallel_axis_flag)
        )
        error.append("stoichiometry_flag = " + str(stoichiometry_flag))
        error.append("decomposition_flag = " + str(decomposition_flag))
        return error

    # these variables are only defined in the Gui mimic if at least one of the above flags is True, so set them after flags are checked
    number_of_contrast_points = variables["number_of_contrast_points"][0]
    fraction_d2o = variables["fraction_d2o"][0]

    # check run_name

    error = input_filter.check_name(run_name)
    if error != []:
        return error

    # check read_from_contrast_calculator_output_file (Will an error be raised prior to this point?  YES! The program will crash in the input filter, so this test isn"t needed.)
    if type(read_from_contrast_calculator_output_file) is not bool:
        error.append(
            "read_from_contrast_calculator_output_file must be True or False")
        return error

    # check that contrast calculator output file exists
    if read_from_contrast_calculator_output_file:
        error = input_filter.check_file_exists(
            contrast_calculator_output_file_name)
        if len(error) > 0:
            error.append(
                contrast_calculator_output_file_name + " is not readable or does not exist"
            )
            return error

    if len(fraction_d2o) != number_of_contrast_points:
        error.append("fraction D2O must have %i values" %
                     (number_of_contrast_points))
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

        # TODO: read_from_contrast_calculator_output_file = True: add option to read izero from a contrast calculator output file.  Errors must be nonzero so we need a way to handle this.

        # check if length of izero, izero_error, concentration and concentration_error = number of contrasts
        if len(izero) != number_of_contrast_points:
            error.append("I(0) must have %i values" %
                         (number_of_contrast_points))
            return error
        if len(izero_error) != number_of_contrast_points:
            error.append("I(0) error must have %i values" %
                         (number_of_contrast_points))
            return error
        if len(concentration) != number_of_contrast_points:
            error.append(
                "concentration must have %i values" % (
                    number_of_contrast_points)
            )
            return error
        if len(concentration_error) != number_of_contrast_points:
            error.append(
                "concentration error must have %i values" % (
                    number_of_contrast_points)
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
        read_from_sascalc_output_file = variables["read_from_sascalc_output_file"][0]
        if read_from_sascalc_output_file:
            sascalc_output_file_name = variables[
                "sascalc_output_file_name"
            ][0]

        # check if number of contrast points is >= 3 to solve for R1, R2 and D
        if number_of_contrast_points < 3:
            error.append("number of contrasts must be >= 3")
            return error

        # TODO: read_from_sascalc_output_file = True: add option to read radius_of_gyration from a sascalc output file.  Errors must be nonzero so we need a way to handle this.
        # TODO: read_from_contrast_calculator_output_file = True: add option to read delta_rho, molecular_weight and partial_specific_volume from a contrast_calculator_output file.

        # NOTE: If either of the above flags is True, the relevant values are read from a contrast calculator or sascalc output file. If we execute the reading of the values from the file and fill them into the GUI before executing the main program, then we can perform the checks on delta_rho whether the flag is True or False and the if statement below can be removed. As of 5/2023, we aren't checking any flags before doing the molecular weight or partial specific volume checks.

        # delta_rho_checks if the values are input by hand
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
                        "delta rho[%i] must have %i values" % (
                            i, number_of_components)
                    )
                    return error
        # check if length of partial specific volume = number of components
        if len(partial_specific_volume) != number_of_components:
            error.append(
                "partial_specific_volume must have %i values" % (
                    number_of_components)
            )
            return error
        # check if length of molecular weight = number of components
        if len(molecular_weight) != number_of_components:
            error.append("Mw must have %i values" % (number_of_components))
            return error
        # check if length of radius_of_gyration and radius_of_gyration_error = number of contrasts
        if len(radius_of_gyration) != number_of_contrast_points:
            error.append("Rg must have %i values" %
                         (number_of_contrast_points))
            return error
        if len(radius_of_gyration_error) != number_of_contrast_points:
            error.append("Rg error must have %i values" %
                         (number_of_contrast_points))
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

        # check if number of contrast points is >= the number of components to solve for the Mw of each component
        if number_of_contrast_points < number_of_components:
            error.append("number of contrasts must be >= number of components (%i)" % (
                number_of_components))
            return error

        # TODO: read_from_contrast_calculator_output_file = True: add option to read izero from a contrast calculator output file.

        # check if length of izero, concentration and delta_rho = number of contrasts
        if len(izero) != number_of_contrast_points:
            error.append("I(0) must have %i values" %
                         (number_of_contrast_points))
            return error
        if len(concentration) != number_of_contrast_points:
            error.append(
                "concentration must have %i values" % (
                    number_of_contrast_points)
            )
            return error

        # TODO: read_from_contrast_calculator_output_file = True: add option to read delta_rho and partial_specific_volume from a contrast_calculator_output file.

        # delta_rho_checks if the values are input by hand
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
                        "delta rho[%i] must have %i values" % (
                            i, number_of_components)
                    )
                    return error
        # check if length of partial specific volume = number of components
        if len(partial_specific_volume) != number_of_components:
            error.append(
                "partial_specific_volume must have %i values" % (
                    number_of_components)
            )
            return error

    elif decomposition_flag:
        number_of_components = variables["number_of_components"][0]
        molecular_weight = variables["molecular_weight"][0]
        partial_specific_volume = variables["partial_specific_volume"][0]
        concentration = variables["concentration"][0]
        # concentration_error is not currently used
        # concentration_error = variables["concentration_error"][0]
        delta_rho = variables["delta_rho"][0]
        data_file_name = variables["data_file_name"][0]
        q_rg_limit_guinier = variables["q_rg_limit_guinier"][0]
        starting_data_point_guinier = variables["starting_data_point_guinier"][0]
        initial_points_to_use_guinier = variables["initial_points_to_use_guinier"][0]
        refine_scale_factor_flag = variables["refine_scale_factor_flag"][0]

        # check if number of contrast points is >= 3 to solve for I11, I12 and I22
        if number_of_contrast_points < 3:
            error.append("number of contrasts must be >= 3")
            return error

        # TODO: read_from_contrast_calculator_output_file = True: add option to read delta_rho, molecular_weight and partial_specific_volume from a contrast_calculator_output file.

        # delta_rho_checks if the values are input by hand
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
                        "delta rho[%i] must have %i values" % (
                            i, number_of_components)
                    )
                    return error
        # check if length of partial specific volume = number of components
        if len(partial_specific_volume) != number_of_components:
            error.append(
                "partial_specific_volume must have %i values" % (
                    number_of_components)
            )
            return error
        # check if length of molecular weight = number of components
        if len(molecular_weight) != number_of_components:
            error.append("Mw must have %i values" % (number_of_components)
                         )
            return error
        # check if length of concentration  = number of contrasts
        if len(concentration) != number_of_contrast_points:
            error.append("concentration must have %i values" % (number_of_contrast_points)
                         )
            return error
        # check if concentration_error is non-zero
        # concentration_error is not currently used
        # for i in range(number_of_contrast_points):
        #    if concentration_error[i] == 0.0:
        #        error.append("concentration error[%i] cannot equal zero" % (i))
        #        return error

        # check if length of data file name = number of contrasts
        if len(data_file_name) != number_of_contrast_points:
            error.append(
                "data file name array must have %i values" % (
                    number_of_contrast_points)
            )
            return error
        # check if data files exist
        for i in range(number_of_contrast_points):
            error = input_filter.check_file_exists(data_file_name[i])
            if len(error) > 0:
                error.append(
                    data_file_name[i] + " is not readable or does not exist"
                )
                return error
        # check if number of points and q binning are the same for all data files
        error = check_data_file(data_file_name)
#        print("error after data file check: ", error)
        if len(error) > 0:
            #            print("error = ", error)
            return error

        # check if length of starting data points and points to use for the Guinier fit = number of contrasts
        if len(starting_data_point_guinier) != number_of_contrast_points:
            error.append("starting data point for Guinier analysis must have %i values" % (number_of_contrast_points)
                         )
            return error

        # check if length of starting data points and points to use for the Guinier fit = number of contrasts
        if len(starting_data_point_guinier) != number_of_contrast_points:
            error.append("starting data points for Guinier analysis much have %i values" % (number_of_contrast_points)
                         )
            return error
        if len(initial_points_to_use_guinier) != number_of_contrast_points:
            error.append("initial number of data points to use for Guinier analysis much have %i values" % (number_of_contrast_points)
                         )
            return error

       # check if qRg limit for the Guinier analysis is > 0
        if q_rg_limit_guinier <= 0.0:
            error.append(
                "qRg limit for Guinier analysis must be greater than zero")
            return error

        # check if refine scale factor flag is boolean
        if type(refine_scale_factor_flag) is not bool:
            error.append("refine_scale_factor_flag is not a boolean")
            return error

    return error


if __name__ == "__main__":
    variables = {}

    # variables["match_point_flag"] = (True, "boolean")
    variables["match_point_flag"] = (False, "boolean")
    # variables["stuhrmann_parallel_axis_flag"] = (True, "boolean")
    variables["stuhrmann_parallel_axis_flag"] = (False, "boolean")
    # variables["stoichiometry_flag"] = (True, "boolean")
    variables["stoichiometry_flag"] = (False, "boolean")
    variables["decomposition_flag"] = (True, "boolean")
    # variables["decomposition_flag"] = (False, "boolean")

    variables["run_name"] = ("run_0", "string")
    variables["output_file_name"] = ("test.out" "string")
    variables["contrast_calculator_output_file_name"] = (
        "./input_contrast.txt", "string")
    variables["read_from_contrast_calculator_output_file"] = (False, "boolean")

    if variables["match_point_flag"][0]:

        variables["number_of_contrast_points"] = (7, "int")
        variables["fraction_d2o"] = (
            [1.0, 0.9, 0.8, 0.4, 0.2, 0.1, 0.0], "int_array")

        variables["izero"] = ([0.537, 0.332, 0.19, 0.0745,
                              0.223, 0.352, 0.541], "float_array")
        variables["izero_error"] = (
            [0.001, 0.002, 0.001, 0.002, 0.002, 0.002, 0.003], "float_array")
        variables["concentration"] = (
            [11.9, 11.9, 11.9, 26.9, 11.9, 11.9, 11.9], "float_array")
        variables["concentration_error"] = (
            [0.6, 0.6, 0.6, 1.3, 0.6, 0.6, 0.6], "float_array")
        variables["initial_match_point_guess"] = (0.5, "float")

    elif variables["stuhrmann_parallel_axis_flag"][0]:

        variables["read_from_sascalc_output_file"] = (False, "boolean")

        variables["number_of_contrast_points"] = (7, "int")
        variables["fraction_d2o"] = (
            [1.0, 0.9, 0.8, 0.4, 0.2, 0.1, 0.0], "int_array")

        variables["number_of_components"] = (2, "int")
        variables["partial_specific_volume"] = ([0.73, 0.73], "float_array")
        variables["molecular_weight"] = ([50.7, 11.7], "float_array")
        variables["delta_rho"] = ([[-3.34, 0.41], [-2.78, 0.98], [-2.21, 1.54], [0.055, 3.8], [
            1.18, 4.93], [1.75, 5.49], [2.31, 6.06]], "nested_float_array")
        variables["radius_of_gyration"] = (
            [25.11, 24.16, 23.03, 23.4, 28.22, 28.33, 28.85], "float_array")
        variables["radius_of_gyration_error"] = (
            [0.09, 0.14, 0.2, 0.7, 0.29, 0.19, 0.12], "float_array")

    elif variables["stoichiometry_flag"][0]:

        variables["number_of_contrast_points"] = (3, "int")
        variables["fraction_d2o"] = ([0.99, 0.12, 0.41], "float_array")
        variables["concentration"] = ([3.7, 3.6, 3.1], "float_array")
        variables["concentration_error"] = ([0.18, 0.18, 0.18], "float_array")
        variables["number_of_components"] = (2, "int")
        variables["partial_specific_volume"] = ([0.745, 0.903], "float_array")
        variables["delta_rho"] = (
            [[-3.2, -5.7], [1.6, 0.26], [0.031, -1.74]], "nested_float_array")
        variables["izero"] = ([8.4, 0.6, 0.17], "float_array")
        variables["izero_error"] = ([0.2, 0.04, 0.01], "float_array")

    elif variables["decomposition_flag"][0]:

        variables["number_of_contrast_points"] = (7, "int")
        variables["fraction_d2o"] = (
            [0.0, 0.1, 0.2, 0.4, 0.8, 0.9, 1.0], "int_array")
        variables["concentration"] = (
            [11.9, 11.9, 11.9, 26.9, 11.9, 11.9, 11.9], "float_array")
        # variables["concentration_error"] = ([0.6, 0.6, 0.6, 1.3, 0.6, 0.6, 0.6], "float_array")
        variables["data_file_name"] = (
            ["./0.dat", "./10.dat", "./20.dat", "./40.dat", "./80.dat", "./90.dat", "./100.dat"], "string_array")
        variables["number_of_components"] = (2, "int")
        variables["partial_specific_volume"] = ([0.73, 0.73], "float_array")
        variables["molecular_weight"] = ([50.7, 11.7], "float_array")
        variables["delta_rho"] = ([[2.31, 6.06], [1.75, 5.49], [1.18, 4.93], [
            0.055, 3.8], [-2.21, 1.54], [-2.78, 0.98], [-3.34, 0.41]], "nested_float_array")
        variables["q_rg_limit_guinier"] = (1.3, "float")
        variables["starting_data_point_guinier"] = (
            [1, 1, 1, 1, 1, 1, 1], "int_array")
        variables["initial_points_to_use_guinier"] = (
            [6, 6, 6, 6, 6, 6, 6], "int_array")
        variables["refine_scale_factor_flag"] = (True, "boolean")

    error = check_multi_component_analysis(variables)
    if (len(error)) == 0:
        print("NO ERRORS FOUND")
    else:
        print("len(__main__ error) = " + str(len(error)))
        from pprint import pprint

        pprint(error)
