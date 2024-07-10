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
import os
import string
import locale
import numpy
import io

import sasmol.system as system


def read_structure_files(other_self):

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    log.debug('in read_structure_files')
    mvars = other_self.module_variables
    avars = other_self.asaxs_variables

    pdb_file_name = mvars.pdb_file_name

    mol = system.Molecule(pdb_file_name)

    print("natoms = ", mol.natoms())

    return


def read_crysol_files(other_self):

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    log.debug('in read_structure_files')
    mvars = other_self.module_variables

    crysol_file_name = mvars.crysol_file_name

# Open the file for reading
    with open(crysol_file_name, 'r') as file:
        # Initialize an empty list to store the data
        data = []

        # Read the file line by line
        for line in file:
            # Split the line into components based on whitespace
            components = line.split()

            # Check if the line has the correct number of components (5)
            if len(components) == 5:
                # Convert the first component to float and name it 's' after dividing by 2*PI
                s = float(components[0])/(2.0*numpy.pi)
                # Convert the remaining components to float and store them as intensities
                intensities = [float(component) for component in components[1:]]
                # Append the 's' and intensities as a tuple to the data list
                data.append((s, *intensities))

    # Example: Print the first few rows to verify
    for row in data[:5]:
        print(f"s: {row[0]}, intensity_1: {row[1]}, intensity_2: {row[2]}, intensity_3: {row[3]}, intensity_4: {row[4]}")


    # use q*2Pi, intensity_1 (total scattering with hydration), intensity_2 (vacuum scattering)

