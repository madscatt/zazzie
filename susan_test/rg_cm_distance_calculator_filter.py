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
#       RG CM DISTANCE CALCULATOR FILTER
#
#       08/20/2024       --      initial coding         :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
"""
    **Rg CM Distance Calculator Filter** is the method that checks the inputs for
    the **Rg CM Distance Calculator** module that were not previously checked by
    **Input Filter**, which only checks for valid string, float, integer, boolean, etc.

    **Inputs:**

        reference PDB file name, trajectory file name, component names, basis string with VMD syntax defining selection basis for each component


    **Outputs:**

        Error string

    Called from **Gui Mimic Rg CM Distance Calculator**

    Requires: **sasmol**, **Input Filter**, **Basis to Python Filter**, **Basis to Python**

"""

import sasmol.system as system
# import input_filter as input_filter
# import basis_to_python_filter as basis_to_python_filter
# import basis_to_python as basis_to_python
import sassie.interface.input_filter as input_filter
import sassie.interface.basis_to_python_filter as basis_to_python_filter
import sassie.util.basis_to_python as basis_to_python


def check_rg_cm_distance_calculator(variables, **kwargs):
    """
    Method to check the **Rg CM Distance Calculator** variables.

    Calls **Input Filter**, **Basis to Python** and **Basis to Python Filter**

    Parameters
    ----------

        Parameters
        ----------

        run_name: string
            run name
        pdb_file_name: string
            name of the reference PDB file
        trajectory_file_name: string
            name of the trajectory file (PDB or DCD)
        number_of_components: int
            number of components in the molecule
        component_name: string array (dimension = number_of_components) 
            names of the components in the molecule
        basis_string: string array (dimension = number_of_components)
            basis string for each component in the molecule

    Returns
    -------

        error: string
            The error message generated when a check fails. If there are no failures, the error is blank.

    """

    run_name = variables['run_name'][0]
    pdb_file_name = variables['pdb_file_name'][0]
    trajectory_file_name = variables['trajectory_file_name'][0]
    path = variables['path'][0]
    number_of_components = variables['number_of_components'][0]
    component_name = variables['component_name'][0]
    basis_string = variables['basis_string'][0]
#    print('run_name: ', run_name)
#    print('pdb_file_name: ', pdb_file_name)
#    print('trajectory_file_name: ', trajectory_file_name)
#    print('path: ', path)
#    print('number_of_components: ', number_of_components)
#    print('component_name: ', component_name)
#    print('basis_string: ', basis_string)

    error = []

# check the run_name
    error = input_filter.check_name(run_name)
    if (len(error)):
        return error

# check the path and permissions
    if 'no_file_check' not in kwargs:
        ev, rv, wv = input_filter.check_permissions(path)
        if (not ev or not rv or not wv):
            error.append('permission error in input file path ' +
                         path + '  [code = ' + str(ev) + str(rv) + str(wv) + ']')
            if (ev == False):
                error.append('path does not exist')
            elif (rv == False):
                error.append('read permission not allowed')
            elif (wv == False):
                error.append('write permission not allowed')
            return error

# check if files end in pdb or dcd
    if ((trajectory_file_name[-3:] != 'dcd' and trajectory_file_name[-3:] != 'pdb') or pdb_file_name[-3:] != 'pdb'):
        error.append('input file names must end with ".pdb" or ".dcd" ')
        return error

# check if reference pdb file exists
    error = input_filter.check_file_exists(pdb_file_name)
    if (len(error)):
        error.append('input pdb file, '+pdb_file_name+', does not exist')
        return error
# check if reference pdb file is a valid pdb file and if it can be read
    ev, value = input_filter.check_pdb_dcd(pdb_file_name, 'pdb')
# if file doesn't exist (ev = 0), error is returned from check_file_exists above.  So, tests won't get to this if stmt.
#       if(ev == 0):
#            error.append('check input pdb file: ' + pdb_file_name)
#            return error
    if (value == 0):
        error.append('input pdb file, ' + pdb_file_name +
                     ', is not a valid pdb file')
        return error

# the try/except below isn't needed?  If the file can't be read successfully, the previous check would have returned an error.
# kept this check just in case
    try:
        m1 = system.Molecule(0)
        m1.read_pdb(pdb_file_name)
        number_of_frames = m1.number_of_frames()
#        print('> found ' + str(number_of_frames) + ' frames in PDB file')
    except:
        error.append('could not read frames in pdb file ', + pdb_file_name)
        return error

# check if trajectory pdb or dcd file exists
    error = input_filter.check_file_exists(trajectory_file_name)
    if (len(error)):
        error.append('input trajectory file, ' +
                     trajectory_file_name+', does not exist')
        return error

# check if trajectory pdb or dcd file is a valid pdb or dcd file and whether it can be read
    ev, value = input_filter.check_pdb_dcd(trajectory_file_name, 'dcd')
# if file doesn't exist (ev = 0), error is returned from check_file_exists above.  So, tests won't get to this if stmt.
#    if(ev == 0):
#        error.append('check input trajectory filename : ' +
#                         trajectory_filename)
#        return error
    if (value == 0):
        ev, value = input_filter.check_pdb_dcd(trajectory_file_name, 'pdb')
        if (value == 0):
            error.append('input trajectory file, ' +
                         trajectory_file_name + ', is not a valid pdb file')
            return error
        input_file_type = 'pdb'
    else:
        input_file_type = 'dcd'

# check if the reference pdb file and trajectory file are compatible
    if (input_file_type == 'dcd'):
        value = input_filter.certify_pdb_dcd(
            pdb_file_name, trajectory_file_name)
        if (value == 0):
            error.append('input pdb file ' + pdb_file_name + ' and trajectory file ' +
                         trajectory_file_name + ', are not compatiable')
            return error
    elif (input_file_type == 'pdb'):
        # corrected error when calling input_filter.certify_pdb_pdb:
        # method returns both ev and value -- not just value
        ev, value = input_filter.certify_pdb_pdb(
            pdb_file_name, trajectory_file_name)
        if (value == 0):
            error.append('input pdb file ' + pdb_file_name + ' and trajectory file ' +
                         trajectory_file_name + ', are not compatiable')
            return error

# the try/except below isn't needed?  If the trajectory file can't be read successfully, the previous check would have returned an error.
# kept this check just in case
    try:
        m1 = system.Molecule(0)
        if (input_file_type == 'dcd'):
            read_trajectory_file = m1.open_dcd_read(trajectory_file_name)
            number_of_frames = read_trajectory_file[2]
#            print('> found ' + str(number_of_frames) +
#                  ' frames in the trajectory DCD file')
        elif (input_file_type == 'pdb'):
            m1.read_pdb(trajectory_file_name)
            number_of_frames = m1.number_of_frames()
#            print('> found ' + str(number_of_frames) +
#                  ' frames in the trajectory PDB file')
    except:
        if (input_file_type == 'pdb'):
            error.append(
                'could not read frames in the trajectory PDB file ', + trajectory_file_name)
        elif (input_file_type == 'dcd'):
            error.append(
                'could not read frames in the trajectory DCD file ', + trajectory_file_name)
        return error

# check that component name has length equal to number of components
    if len(component_name) != number_of_components:
        error.append('number of components (' + str(number_of_components) +
                     ') does not match length of component names (' + str(len(component_name)) + ')')
        return error
# check that the basis string has length equal to number of components
    if len(basis_string) != number_of_components:
        error.append('number of components (' + str(number_of_components) +
                     ') does not match length of basis strings (' + str(len(basis_string)) + ')')
        return error

# check the basis string

    for i in range(len(basis_string)):
        #        print('checking basis: ' + basis_string[i])
        error = basis_to_python_filter.check_basis_syntax(basis_string[i])
#        print('len(error): ', len(error))
        if (len(error)):
            return error
        try:
            this_python_basis = basis_to_python.parse_basis(basis_string[i])
#            print('python_basis: ', this_python_basis)
        except:
            error.append(
                'unable to convert input string #' + basis_string[i] + '# to python readable string')
            return error

    return error


if __name__ == "__main__":
    variables = {}

    variables["path"] = ("./", "string")
    variables["run_name"] = ("paivn_0", "string")
    variables["pdb_file_name"] = ("./pai_vn_start.pdb", "string")
#    variables["pdb_file_name"] = ("./any_old_file.pdb", "string")
#    variables["trajectory_file_name"] = ("./test.pdb", "string")
#    variables["trajectory_file_name"] = ("./pai_vn_20_frames.dcd", "string")
#    variables["trajectory_file_name"] = ("./pai_vn_20_frames.pdb", "string")
#    variables["trajectory_file_name"] = ("./pai_vn_2_frames.pdb", "string")
    variables["trajectory_file_name"] = ("./best_all.dcd", "string")
    variables["number_of_components"] = (2, "int")
    variables["component_name"] = (["VN", "PAI"], "string_array")
    variables["basis_string"] = (
        ["segname VN1", "segname PAI1"], "string_array")

    error = check_rg_cm_distance_calculator(variables)
    if (len(error)) == 0:
        print("NO ERRORS FOUND")
    else:
        print("len(__main__ error) = " + str(len(error)))
        from pprint import pprint

        pprint(error)
