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
#       GUI MIMIC RG CM DISTANCE CALCULATOR
#
#       08/20/2024       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
"""
    **Gui Mimic Rg CM Distance Calculator** is the driver method to run the
    **Rg CM Distance Calculator** module.

    Requires **Input Filter**, **Rg CM Distance Calculator Filter** and **Rg CM Distance Calculator**.

"""

import sys
import os
import shutil
import time
import multiprocessing
import sassie.contrast.rg_cm_distance_calculator.rg_cm_distance_calculator as rg_cm_distance_calculator
import sassie.interface.input_filter as input_filter
import sassie.interface.rg_cm_distance_calculator.rg_cm_distance_calculator_filter as rg_cm_distance_calculator_filter
# import rg_cm_distance_calculator as rg_cm_distance_calculator
# import rg_cm_distance_calculator_filter as rg_cm_distance_calculator_filter
# import input_filter as input_filter


def user_variables(self, **kwargs):
    """
    Method for the user to enter the input variables.

    """

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.run_name = 'run_0'
    self.path = './'
#    self.path = ('./test/')
    self.pdb_file_name = self.path+'pai_vn_start.pdb'
#    self.trajectory_file_name = self.path+'any_old_file.pdb'
#    self.pdb_file_name = self.path+'some_file.pdb'
#    self.trajectory_file_name = self.path+'pai_vn_start.pdb'
    self.trajectory_file_name = self.path+'pai_vn_2_frames.pdb'
#   self.trajectory_file_name = self.path+'pai_vn_20_frames.dcd'
#    self.trajectory_file_name = self.path+'best_all.dcd'
#    self.trajectory_file_name = self.path+'struct_k6_dimer.dcd'  # this should fail
    self.number_of_components = '2'
    self.component_name = 'VN, PAI'
    self.basis_string = 'segname VN1, segname PAI1'
#    self.basis_string = 'segname = VN1, segname PAI1'  # this should fail

    self.testflag = False

    ### END USER INPUT ###
    ### END USER INPUT ###
    ### END USER INPUT ###


def test_variables(self, paths):
    """

    Method that defines the variables that will be used to test the module as well as its input filter. Variables are defined outside the gui_mimic_rg_cm_distance_calculator class so that they can be used by these other programs.

    Users of gui_mimic as a driver script to run this module should not edit the values below as they are used for development tests.

    """

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    other_data_path = paths['other_data_path']
    module_data_path = paths['module_data_path']

    self.run_name = 'run_0'
    self.path = './'
    self.pdb_file_name = os.path.join((self.path), 'pai_vn_start.pdb')
#    self.trajectory_file_name = os.path.join((self.path), 'pai_vn_2_frames.pdb')
    self.trajectory_file_name = os.path.join((self.path), 'best_all.dcd')
    self.number_of_components = "2"  # hardwired at this time
    self.component_name = 'VN, PAI'
    self.basis_string = 'segname VN1, segname PAI1'

    self.testflag = True
    self.precision = 3  # is this needed?


def run_module(self, **kwargs):
    """
    Method to run the module and/or its input filter.
    Only the module input filter is run if kwargs is: test_filter=True
    The method is defined outside the class so that it can be used by other programs such as test_module and test_module_filter.

    """

    svariables = {}

    svariables['run_name'] = (self.run_name, 'string')
    svariables['pdb_file_name'] = (self.pdb_file_name, 'string')
    svariables['trajectory_file_name'] = (self.trajectory_file_name, 'string')
    svariables["number_of_components"] = (self.number_of_components, "int")
    svariables["component_name"] = (self.component_name, "string_array")
    svariables["basis_string"] = (self.basis_string, "string_array")
    svariables['path'] = (self.path, 'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)
    if (len(error) > 0):
        print('error = ', error)
        if not (self.testflag):
            sys.exit()
        return error

    print(self.variables)
    # sys.exit()

    try:
        if kwargs['file_check']:
            error = rg_cm_distance_calculator_filter.check_rg_cm_distance_calculator(
                self.variables)
    except:
        error = rg_cm_distance_calculator_filter.check_rg_cm_distance_calculator(
            self.variables, no_file_check="true")

    if (len(error) > 0):
        print('error = ', error)
        if not (self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        print('not testing')
        pass

    run_name = self.variables['run_name'][0]

#    if os.path.exists(os.path.join(run_name, self.module)):
#        shutil.rmtree(os.path.join(run_name, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_rg_cm_distance_calculator = rg_cm_distance_calculator.rg_cm_distance_calculator()
    this_rg_cm_distance_calculator.main(self.variables, txtQueue)


class gui_mimic_rg_cm_distance_calculator():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'rg_cm_distance_calculator'

    def __init__(self, test, paths):

        if not test:
            user_variables(self)
        else:
            test_variables(self, paths)

#        run_module(self, test_filter=True)
        run_module(self)


if __name__ == '__main__':

    test = False  # option to run with test variables not implemented in 1.0.
    paths = None


# We are thinking of defining the install path so the gui mimic can be run from anywhere as long as it is called from that particular python
# That way, the test files will always be available to the user.
    if test:
        pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
        dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
        other_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'other_common') + os.path.sep
        module_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'interface', 'extract_utilities') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_rg_cm_distance_calculator(test, paths)
    print("time used: ", time.time() - start)
