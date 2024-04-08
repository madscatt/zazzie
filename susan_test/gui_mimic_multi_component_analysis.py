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
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#       MULTI-COMPONENT ANALYSIS
#
#       08/09/2021       --      initial coding         :   Susan Krueger
#       03/13/2023       --      python 3               :   Joseph E. Curtis
#       03/27/2024       --   use PAI-VN and RSV params :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
"""
    **Gui Mimic Multi-component Analysis** is the driver method to run the
    **Multi-compoment Analysis** module.

    Requires **Input Filter**, **Multi-component Analysis Filter** and **Multi-component Analysis**.

"""

import os
import shutil
import time

# import sassie.contrast.multi_component_analysis.multi_component_analysis as multi_component_analysis
# import sassie.interface.input_filter as input_filter
# import sassie.interface.multi_component_analysis.multi_component_analysis_filter as multi_component_analysis_filter
import multi_component_analysis as multi_component_analysis
import input_filter_new as input_filter
import multi_component_analysis_filter as multi_component_analysis_filter
import multiprocessing


def user_variables(self, **kwargs):
    """
    Method for the user to enter the input variables.

    Note: if data is read from file (either from constrast calculator or sascalc)
    the data are used to "autofill" the values in the GUI to be then checked and/or
    altered by the user.

    """

    #### user input ####
    #### user input ####
    #### user input ####

    # METHOD SELECTION FLAGS
    # ONLY ONE OF THE FOLLOWING FOUR SHOULD BE TRUE

    self.path = "./"
#    self.match_point_flag = True
    self.match_point_flag = False
#    self.stuhrmann_parallel_axis_flag = True
    self.stuhrmann_parallel_axis_flag = False
    self.decomposition_flag = True
#    self.decomposition_flag = False
#    self.stoichiometry_flag = True
    self.stoichiometry_flag = False

    # INPUT VARIABLES FOR ALL METHODS

    self.run_name = "run_0"
    self.output_file_name = "general_output_file.out"

#    self.read_from_contrast_calculator_output_file = True
    self.read_from_contrast_calculator_output_file = False

    if self.match_point_flag or self.stuhrmann_parallel_axis_flag or self.decomposition_flag:
        self.number_of_contrast_points = "4"
        self.fraction_d2o = "0.0, 0.2, 0.85, 1.0"
        if self.read_from_contrast_calculator_output_file:
            self.contrast_calculator_output_file_name = os.path.join(
                self.path, "sample_contrast_calculator_output_file.txt"
            )
    elif self.stoichiometry_flag:
        self.number_of_contrast_points = "3"
        self.fraction_d2o = "0.99, 0.12, 0.41"
        if self.read_from_contrast_calculator_output_file:
            self.contrast_calculator_output_file_name = os.path.join(
                self.path, "sample_contrast_calculator_output_file1.txt"
            )  # note that this file is different from the one above since it is for RSV and not PAI-VN

    # MATCH POINT ANALYSIS VARIABLES

    if self.match_point_flag:
        self.concentration = "7.7, 7.7, 7.7, 7.7"
        self.concentration_error = "0.4, 0.4, 0.4, 0.4"
        # False is the default for the intial_match_point_guess_flag; an initial guess will be made by fitting a polynomial to I(0)/c vs fraction D2O
        # NOTE: allow the user to change to True as an advanced option in case the automated attempt fails to find a reasonable value
        self.initial_match_point_guess_flag = False
        # initial_match_point_guess is only used if intial_match_point_guess_flag = True
        if self.initial_match_point_guess_flag:
            self.initial_match_point_guess = "0.5"

        if self.read_from_contrast_calculator_output_file:
            # TODO NEED TO READ THESE IN; VALUES ARE HERE AS A PLACEHOLDER
            self.izero = "0.85, 0.534, 0.013, 0.095"
            self.izero_error = "0.01, 0.044, 0.003, 0.002"

        else:
            self.izero = "0.85, 0.534, 0.013, 0.095"
            self.izero_error = "0.01, 0.044, 0.003, 0.002"

    # STUHRMANN PARALLEL AXIS ANALYSIS VARIABLES

    elif self.stuhrmann_parallel_axis_flag:
        self.number_of_components = "2"
        self.component_name = "vn, pai"
        # intial guess for Stuhrmann equation coefficients
        # NOTE: user shouldn't have to change the initial guess; make it an advanced option just in case?
        self.initial_guess_stuhrmann = "1.0, 1.0, 1.0"
        self.read_from_sascalc_output_file = False

        if self.read_from_contrast_calculator_output_file:
            # TODO NEED TO READ THESE IN; VALUES ARE HERE AS A PLACEHOLDER
            self.partial_specific_volume = "0.73, 0.73"  # 1 value for each component
            self.molecular_weight = "14.3, 44.2"  # kDa; 1 value for each component
            # 2 values for each contrast since there are 2 components.
            self.delta_rho = "2.551, 5.104; 1.383, 3.928; -2.415, 0.109; -3.292, -0.773"
        else:
            self.partial_specific_volume = "0.73, 0.73"  # 1 value for each component
            self.molecular_weight = "14.3, 44.2"  # kDa; 1 value for each component
            # 2 values for each contrast point since there are 2 components.
            self.delta_rho = "2.551, 5.104; 1.383, 3.928; -2.415, 0.109; -3.292, -0.773"

        if self.read_from_sascalc_output_file:
            # TODO NEED TO READ THESE IN; VALUES ARE HERE AS A PLACEHOLDER
            self.sascalc_output_file_name = os.path.join(
                self.path, "sample_sascalc_output_file.txt"
            )
            # 1 value for each contrast point
            # these are the Rg values for the starting model structure
            self.radius_of_gyration = "31.1, 27.2, 42.1, 57.8"
            # 1 value for each contrast point
            self.radius_of_gyration_error = "0.1, 0.2, 1.0, 0.5"

        else:
            # 1 value for each contrast point
            # these are the experimental Rg values
            self.radius_of_gyration = "25.45, 24.95, 28.0, 31.34"
            # 1 value for each contrast point
            self.radius_of_gyration_error = "0.07, 0.09, 3.0, 0.4"

    # DECOMPOSITION ANALYSIS VARIABLES

    elif self.decomposition_flag:
        # NOTE: The path needs to be part of the file name so that the check as to whether each file exists will look for the files in the right place.
        self.data_file_name = self.path+"0p.dat, "+self.path + \
            "20p.dat, " + self.path+"85p1.dat, "+self.path+"100p1.dat"
#        print('data file name: ', self.data_file_name)
        self.concentration = "7.7, 7.7, 7.7, 7.7"
        self.concentration_error = "0.4, 0.4, 0.4, 0.4"
        self.number_of_components = "2"
        self.component_name = "vn, pai"
        self.q_rg_limit_guinier = "1.3"
        self.starting_data_point_guinier = "1, 1, 1, 1"
#        self.initial_points_to_use_guinier = "4, 4, 4, 3"  #for model data
        self.initial_points_to_use_guinier = "6, 6, 6, 6"
        self.refine_scale_factor_flag = True
        # self.refine_scale_factor_flag = False
        # intial guess for Guinier equation coefficients
        # NOTE: user shouldn't have to change the initial guess; make it an advanced option just in case?
        self.initial_guess_guinier = "1.0, 1.0"
        # amplitude for the Gaussian that describes the S/N of typical SANS data; only used if there are no errors specified in the data file
        # the values below are the defaults; they can be changed by the user to vary as a function of contrast; only used for model data; treat like an advanced option?
#        self.sn_amplitude = "50.0, 50.0, 50.0, 50.0"  # default
        self.sn_amplitude = "400.0, 300.0, 10.0, 100.0"  # PAI-VN CV series

        if self.read_from_contrast_calculator_output_file:
            # TODO NEED TO READ THESE IN; VALUES ARE HERE AS A PLACEHOLDER

            self.partial_specific_volume = "0.73, 0.73"  # 1 value for each component
            self.molecular_weight = "14.3, 44.2"  # kDa; 1 value for each component
            # 2 values for each contrast since there are 2 components.
            self.delta_rho = "2.551, 5.104; 1.383, 3.928; -2.415, 0.109; -3.292, -0.773"

        else:
            self.partial_specific_volume = "0.73, 0.73"  # 1 value for each component
            self.molecular_weight = "14.3, 44.2"  # kDa; 1 value for each component
            # 2 values for each contrast since there are 2 components.
            self.delta_rho = "2.551, 5.104; 1.383, 3.928; -2.415, 0.109; -3.292, -0.773"

    # STOICHIOMETRY ANALYSIS VARIABLES

    elif self.stoichiometry_flag:
        self.concentration = "3.1, 3.6, 3.1"
        self.concentration_error = "0.18, 0.18, 0.18"
        self.number_of_components = "2"
        self.component_name = "rsv, ps80"

        if self.read_from_contrast_calculator_output_file:
            # TODO NEED TO READ THESE IN; VALUES ARE HERE AS A PLACEHOLDER
            self.partial_specific_volume = "0.745, 0.903"  # 1 value for each component
            # 2 values for each contrast point since there are 2 components.
            self.delta_rho = "-3.2, -5.7; 1.6, 0.26; 0.031, -1.74"
            self.izero = "8.1, 0.6, 0.17"  # 1 value for each contrast point
            self.izero_error = "0.2, 0.04, 0.01"  # 1 value for each contrast point

        else:
            self.partial_specific_volume = "0.745, 0.903"  # 1 value for each component
            # 2 values for each contrast point since there are 2 components.
            self.delta_rho = "-3.2, -5.7; 1.6, 0.26; 0.031, -1.74"
            self.izero = "8.4, 0.6, 0.17"  # 1 value for each contrast point
            self.izero_error = "0.2, 0.04, 0.01"  # 1 value for each contrast point

    #### end user input ####
    #### end user input ####
    #### end user input ####


def test_variables(self, paths):
    """

    Method that defines the variables that will be used to test the module as well as its input filter. Variables are defined outside the gui_mimic_multi_component_analysis class so that they can be used by these other programs.

    Users of gui_mimic as a driver script to run this module should not edit the values below as they are used for development tests.

    """

# TODO: test variables need to be updated to include all variables as in the user variables above
    pdb_data_path = paths["pdb_data_path"]
    dcd_data_path = paths["dcd_data_path"]
    module_data_path = paths["module_data_path"]
    other_data_path = paths["other_data_path"]

    self.run_name = "run_0"
    self.path = ""
    self.number_of_contrast_points = "3"
    self.number_of_components = "2"
    self.read_from_contrast_calculator_output_file = False
    if self.read_from_contrast_calculator_output_file:
        self.contrast_calculator_output_file_name = os.path.join(
            other_data_path, "input_contrast.txt")
    self.output_file_name = os.path.join(other_data_path, "99_12_41.out")
    self.match_point_flag = False
    self.stuhrmann_parallel_axis_flag = False
    self.decomposition_flag = False
    self.stoichiometry_flag = True
    self.fraction_d2o = "0.99, 0.12, 0.41"
    self.concentration = "3.7, 3.6, 3.1"
    self.concentration_error = "0.4, 0.4, 0.3"
    self.izero = "11.8, 0.6, 0.17"
    self.izero_error = "0.1, 0.04, 0.03"
    self.partial_specific_volume = "0.745, 0.903"
    self.delta_rho = "-3.2,-5.7; 1.6, 0.26; 0.031, -1.74"

    self.precision = 3  # not needed? What is the default?


def run_module(self, **kwargs):
    """
    Method to run the module and/or its input filter.
    Only the module input filter is run if kwargs is: test_filter=True
    The method is defined outside the class so that it can be used by other programs such as test_module and test_module_filter.

    """

    svariables = {}

    # INPUT VARIABLES FOR ALL METHODS

    svariables["run_name"] = (self.run_name, "string")
    svariables["path"] = (self.path, "string")
    svariables["output_file_name"] = (self.output_file_name, "string")

    svariables["match_point_flag"] = (self.match_point_flag, "boolean")
    svariables["stoichiometry_flag"] = (self.stoichiometry_flag, "boolean")
    svariables["stuhrmann_parallel_axis_flag"] = (
        self.stuhrmann_parallel_axis_flag,
        "boolean",
    )
    svariables["decomposition_flag"] = (self.decomposition_flag, "boolean")

    svariables["read_from_contrast_calculator_output_file"] = (
        self.read_from_contrast_calculator_output_file,
        "boolean",
    )

    if self.read_from_contrast_calculator_output_file:
        svariables["contrast_calculator_output_file_name"] = (
            self.contrast_calculator_output_file_name,
            "string",
        )

    # MATCH POINT ANALYSIS VARIABLES

    if self.match_point_flag:
        svariables["number_of_contrast_points"] = (
            self.number_of_contrast_points, "int")
        svariables["fraction_d2o"] = (self.fraction_d2o, "float_array")
        svariables["initial_match_point_guess_flag"] = (
            self.initial_match_point_guess_flag,
            "boolean",
        )
        if self.initial_match_point_guess_flag:
            svariables["initial_match_point_guess"] = (
                self.initial_match_point_guess,
                "float",
            )
        svariables["concentration"] = (self.concentration, "float_array")
        svariables["concentration_error"] = (
            self.concentration_error, "float_array")

        svariables["izero"] = (self.izero, "float_array")
        svariables["izero_error"] = (self.izero_error, "float_array")

    # STUHRMANN PARALLEL AXIS ANALYSIS VARIABLES

    elif self.stuhrmann_parallel_axis_flag:
        svariables["number_of_contrast_points"] = (
            self.number_of_contrast_points, "int")
        svariables["fraction_d2o"] = (self.fraction_d2o, "float_array")
        svariables["initial_guess_stuhrmann"] = (
            self.initial_guess_stuhrmann, "float_array")
        svariables["number_of_components"] = (self.number_of_components, "int")
        svariables["component_name"] = (self.component_name, "string_array")
# already defined above
#        svariables["read_from_contrast_calculator_output_file"] = (
#            self.read_from_contrast_calculator_output_file,
#            "boolean",
#        )
        svariables["read_from_sascalc_output_file"] = (
            self.read_from_sascalc_output_file,
            "boolean",
        )

        svariables["partial_specific_volume"] = (
            self.partial_specific_volume,
            "float_array",
        )
        svariables["molecular_weight"] = (self.molecular_weight, "float_array")
        svariables["delta_rho"] = (self.delta_rho, "nested_float_array")

        if self.read_from_sascalc_output_file:
            svariables["sascalc_output_file_name"] = (
                self.sascalc_output_file_name,
                "string",
            )

        svariables["radius_of_gyration"] = (
            self.radius_of_gyration, "float_array")
        svariables["radius_of_gyration_error"] = (
            self.radius_of_gyration_error,
            "float_array",
        )

    # DECOMPOSITION ANALYSIS VARIABLES

    elif self.decomposition_flag:
        svariables["number_of_contrast_points"] = (
            self.number_of_contrast_points, "int")
        svariables["fraction_d2o"] = (self.fraction_d2o, "float_array")
        svariables["data_file_name"] = (self.data_file_name, "string_array")
        svariables["concentration"] = (self.concentration, "float_array")
        svariables["concentration_error"] = (
            self.concentration_error, "float_array")
        svariables["number_of_components"] = (self.number_of_components, "int")
        svariables["component_name"] = (self.component_name, "string_array")
        svariables["q_rg_limit_guinier"] = (self.q_rg_limit_guinier, "float")
        svariables["starting_data_point_guinier"] = (
            self.starting_data_point_guinier, "int_array")
        svariables["initial_points_to_use_guinier"] = (
            self.initial_points_to_use_guinier, "int_array")
        svariables["refine_scale_factor_flag"] = (
            self.refine_scale_factor_flag, "boolean")
        svariables["initial_guess_guinier"] = (
            self.initial_guess_guinier, "float_array")
        svariables["sn_amplitude"] = (self.sn_amplitude, "float_array")
# already defined above
#        svariables["read_from_contrast_calculator_output_file"] = (
#            self.read_from_contrast_calculator_output_file,
#            "boolean",
#        )

        svariables["partial_specific_volume"] = (
            self.partial_specific_volume,
            "float_array",
        )
        svariables["molecular_weight"] = (self.molecular_weight, "float_array")
        svariables["delta_rho"] = (self.delta_rho, "nested_float_array")

    # STOICHIOMETRY ANALYSIS VARIABLES

    elif self.stoichiometry_flag:
        svariables["number_of_contrast_points"] = (
            self.number_of_contrast_points, "int")
        svariables["fraction_d2o"] = (self.fraction_d2o, "float_array")
        svariables["concentration"] = (self.concentration, "float_array")
        svariables["concentration_error"] = (
            self.concentration_error, "float_array")
        svariables["number_of_components"] = (self.number_of_components, "int")
        svariables["component_name"] = (self.component_name, "string_array")
# already defined above
#        svariables["read_from_contrast_calculator_output_file"] = (
#            self.read_from_contrast_calculator_output_file,
#            "boolean",
#        )

        svariables["partial_specific_volume"] = (
            self.partial_specific_volume,
            "float_array",
        )
        svariables["delta_rho"] = (self.delta_rho, "nested_float_array")
        svariables["izero"] = (self.izero, "float_array")
        svariables["izero_error"] = (self.izero_error, "float_array")

    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print("error = ", error)
        return error

#    print(self.variables)
    # import sys; sys.exit()

    try:
        if kwargs["file_check"]:
            error = multi_component_analysis_filter.check_multi_component_analysis(
                self.variables
            )
    except:
        error = multi_component_analysis_filter.check_multi_component_analysis(
            self.variables, no_file_check="true"
        )

    if len(error) > 0:
        print("error = ", error)
        return error

    try:
        if kwargs["test_filter"]:
            return error
    except:
        #        print("not testing")
        pass

    run_name = self.variables["run_name"][0]

    # I don't want the directory tree to be deleted since I want to be able to put more than one
    # output file in the directory for scenarios involving different contrasts.
    # But, then the sassie_json and sassie_log files don't get removed even if an output
    # file is overwritten.
    # This will eventually be handled at the GenApp level, so I am leaving it as is here.

    if os.path.exists(os.path.join(run_name, self.module)):
        shutil.rmtree(os.path.join(run_name, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_multi_component_analysis = multi_component_analysis.multi_component_analysis()
    this_multi_component_analysis.main(self.variables, txtQueue)


class gui_mimic_multi_component_analysis:
    """
    gui_mimic class contains the name of the module
    """

    module = "multi_component_analysis"

    def __init__(self, test, paths):
        if not test:
            user_variables(self)
        else:
            test_variables(self, paths)

        # run_module(self, test_filter=True)
        run_module(self)


if __name__ == "__main__":
    test = False  # option to run with test variables not implemented
    paths = None

    # We are thinking of defining the install path so the gui mimic can be run from anywhere as long as it is called from that particular python
    # That way, the test files will always be available to the user.
    if test:
        pdb_data_path = (
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "..",
                "..",
                "data",
                "pdb_common",
            )
            + os.path.sep
        )
        dcd_data_path = (
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "..",
                "..",
                "data",
                "dcd_common",
            )
            + os.path.sep
        )
        module_data_path = (
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "..",
                "..",
                "data",
                "interface",
                "align",
            )
            + os.path.sep
        )

        paths = {
            "pdb_data_path": pdb_data_path,
            "dcd_data_path": dcd_data_path,
            "module_data_path": module_data_path,
        }

    start = time.time()
    run_gui = gui_mimic_multi_component_analysis(test, paths)
    print("time used: ", time.time() - start)
