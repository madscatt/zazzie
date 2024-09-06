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
#       RG CM DISTANCE CALCULATOR
#
#       08/20/2024       --      initial coding               :  Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Rg CM Distance Calculator** is the module that is used to calculate, from a PDB or trajectory (DCD) file, the radius of gyration of each component in a two-component system as well as the distance between their center of mass.

    **Inputs:**

        reference PDB file name, trajectory file name, component names, basis string with VMD-like syntax defining selection basis for each component

    **Outputs:**

        radii of gyration for the two components, distance between their centers of mass

    Called from **Gui Mimic Rg CM Distance Calculator**

    Calls **Basis to Python**

    Requires **sasmol**, **numpy**, **module_utilities**, **sasconfig**, **basis_to_python**

    TODO: Need boilerplate line that shows the flow, i.e.,  module utilities, setup logging, unpack variables, run main, etc. This will be in all modules.

'''

import os
import io
import time
import numpy
import sasmol.system as system
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.basis_to_python as basis_to_python
# import basis_to_python as basis_to_python

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'rg_cm_distance_calculator'


class module_variables():
    ''' Module variables class'''

    def __init__(self, parent=None):
        self.app = app


class rg_cm_distance_calculator_variables():
    ''' Rg CM distance calculator variables class'''

    def __init__(self, parent=None):
        pass


class rg_cm_distance_calculator():
    ''' Base class containing methods that control the program execution, unpacking of the module variables, intialization of the Rg CM distance calculator variables and clean up after execution.'''

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):
        '''The main method that calls all of the other methods.'''

        self.module_variables = module_variables()

        self.rg_cm_distance_calculator_variables = rg_cm_distance_calculator_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.calculate_rg_cm_distance()

        self.epilogue()

        return

    def unpack_variables(self, variables):
        '''Method to unpack variables passed from the GUI.'''

        log = self.log
        mvars = self.module_variables
        log.debug('in unpack_variables')

# Unpack the variables that are used for all methods.
        mvars.run_name = variables['run_name'][0]
        mvars.pdb_file_name = variables['pdb_file_name'][0]
        mvars.trajectory_file_name = variables['trajectory_file_name'][0]
#        mvars.path = variables['path'][0]
        mvars.number_of_components = variables['number_of_components'][0]
        mvars.component_name = variables['component_name'][0]
        mvars.basis_string = variables['basis_string'][0]

#        print(vars(mvars))

        log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        Method to initialize Rg CM Distance Calculator variables and write initial information to an output file.

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

        output_file_path: string
            sub-path where output file will be written: run_name + \'rg_cm_distance_calculator'
        input_file_type: string
            type of trajectory file (PDB or DCD)
        output_file_name: string
            name of the output file (depends on the trajectory file name)
        outfile: string
            output file name (with full path): path + output_file_name
        python_basis: string array (dimension = number_of_components)
            basis string for each component in the molecule in python-readable format

        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        rgcmdvars = self.rg_cm_distance_calculator_variables

# output file path
        if (mvars.run_name[-1] == '/'):
            log.debug('run_name(1) = %s' % (mvars.run_name))
            rgcmdvars.output_file_path = mvars.run_name + 'rg_cm_distance_calculator/'
            log.debug('output_file_path = %s' %
                      (rgcmdvars.output_file_path))
        else:
            log.debug('run_name(2) = %s' % (mvars.run_name))
            rgcmdvars.output_file_path = mvars.run_name + '/rg_cm_distance_calculator/'
            log.debug('output_file_path = %s' %
                      (rgcmdvars.output_file_path))


# get the output file name
# first, strip the input path from the file name since it may be different from the output file path
        numslash = mvars.trajectory_file_name.count('/')
        if numslash == 0:
            stripped_trajectory_file_name = mvars.trajectory_file_name
        else:
            groups = mvars.trajectory_file_name.split('/')
            log.debug('groups = %s' % (groups))
#            print('groups: ', groups)
            stripped_trajectory_file_name = ('.'.join(groups[numslash:]))

#        print('stripped file name: ', stripped_trajectory_file_name)

# output file name is then defined using the stripped dcd file name (output path is added later)
# find the number of '.' in each stripped file name and locate the last one so that more characters can be added to the filename at that location

        numdots = stripped_trajectory_file_name.count('.')
#        print('numdots: ', numdots)
        if numdots == 0:
            if mvars.trajectory_file_name[-3:] == 'dcd':
                rgcmdvars.input_file_type = 'dcd'
                rgcmdvars.output_file_name = stripped_trajectory_file_name + '_rgdist_dcd.txt'
            elif mvars.trajectory_file_name[-3:] == 'pdb':
                rgcmdvars.input_file_type = 'pdb'
                rgcmdvars.output_file_name = stripped_trajectory_file_name + '_rgdist_pdb.txt'
            log.debug('output file name(1) = %s' %
                      (rgcmdvars.output_file_name))
            log.debug('input_file_type = %s' %
                      (rgcmdvars.input_file_type))
#            print('input_file_type, output_file_name(1): ',
#                  rgcmdvars.input_file_type, rgcmdvars.output_file_name)
        else:
            groups = stripped_trajectory_file_name.split('.')
            log.debug('groups = %s' % (groups))
#            print('groups: ', groups)
            if mvars.trajectory_file_name[-3:] == 'dcd':
                rgcmdvars.input_file_type = 'dcd'
                rgcmdvars.output_file_name = (
                    '.'.join(groups[:numdots]) + '_rgdist_dcd.txt')
            elif mvars.trajectory_file_name[-3:] == 'pdb':
                rgcmdvars.input_file_type = 'pdb'
                rgcmdvars.output_file_name = (
                    '.'.join(groups[:numdots]) + '_rgdist_pdb.txt')
            log.debug('output file name(2) = %s' %
                      (rgcmdvars.output_file_name))
            log.debug('input_file_type = %s' %
                      (rgcmdvars.input_file_type))
#            print('input_file_type, output_file_name(2): ',
#                  rgcmdvars.input_file_type, rgcmdvars.output_file_name)

        # check for existence of output file path and create if necessary

        direxist = os.path.exists(rgcmdvars.output_file_path)
        if (direxist == 0):
            os.system('mkdir -p ' + rgcmdvars.output_file_path)

# convert basis to python readable format
# The basis string should be correct, since it was checked in the module filter; however, it is checked again here just in case
        rgcmdvars.python_basis = []
        error = []
        for i in range(len(mvars.basis_string)):
            try:
                this_python_basis = basis_to_python.parse_basis(
                    mvars.basis_string[i])
            except:
                error.append(
                    'unable to convert input string ' + mvars.basis_string[i] + ' to python readable string')
                pgui('error: ', error)
                return error
            rgcmdvars.python_basis.append(this_python_basis)

        log.debug('python basis = %s' % (rgcmdvars.python_basis))
#        print('python_basis: ', rgcmdvars.python_basis)
#        print('length of python_basis: ', len(rgcmdvars.python_basis))

# open the general output file for writing
        rgcmdvars.outfile = io.open(
            rgcmdvars.output_file_path+rgcmdvars.output_file_name, 'w')

        # write the input variables in the output file

        rgcmdvars.outfile.write('input variables: ' + repr(vars(mvars)) + '\n')

        log.debug(vars(mvars))
        log.debug(vars(rgcmdvars))

        return

    def do_the_calculations(self, molecule, structure_number):
        '''
        **Do the Calculations** calculates the radius of gyration of each component in a two-component system as well as the distance between their center of mass for a single structure in a trajectory.

        Called from **Calculate Rg CM Distance**

        Parameters
        ----------

        python_basis: string array (dimension = number_of_components)
            basis string for each component in the molecule in python-readable format
        outfile: string
            output file name (with full path): path + output_file_name
        molecule: system.Molecule
            instance of the molecule class containing the coordinates of the current structure
        structure_number: int
            number of the current structure in the trajectory

        Returns
        -------

        radius_of_gyration_molecule1: float
            radius of gyration of the first component
        radius_of_gyration_molecule2: float
            radius of gyration of the second component
        distance: float
            distance between the centers of mass of the two components

        '''

        #   rgcmdvars used: python_basis, outfile

        log = self.log
        log.debug('performing calculations')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        rgcmdvars = self.rg_cm_distance_calculator_variables

        error = []
        error, mask1 = molecule.get_subset_mask(rgcmdvars.python_basis[0])
        if error:
            pgui('error, mask: ', error, mask1)
            return error
        error, mask2 = molecule.get_subset_mask(rgcmdvars.python_basis[1])
        if error:
            pgui('error, mask2: ', error, mask2)
            return error

        molecule1 = system.Molecule(1)
        molecule2 = system.Molecule(2)

        error = molecule.copy_molecule_using_mask(molecule1, mask1, 0)
        if error:
            pgui('error: ', error)
            return error
        error = molecule.copy_molecule_using_mask(molecule2, mask2, 0)
        if error:
            pgui('error: ', error)
            return error

        center_of_mass_molecule1 = molecule1.calculate_center_of_mass(0)
        center_of_mass_molecule2 = molecule2.calculate_center_of_mass(0)

        radius_of_gyration_molecule1 = molecule1.calculate_radius_of_gyration(
            0)
        radius_of_gyration_molecule2 = molecule2.calculate_radius_of_gyration(
            0)

        log.debug('center_of_mass_molecule1: %s' % (center_of_mass_molecule1))
        log.debug('center_of_mass_molecule2: %s' % (center_of_mass_molecule2))
#       print('center_of_mass_molecule1: ', center_of_mass_molecule1)
#       print('center_of_mass_molecule2: ', center_of_mass_molecule2)

        distance = numpy.sqrt(
            numpy.sum((center_of_mass_molecule1 - center_of_mass_molecule2)**2))

        pgui('distance between centers of mass of %s, %s: %0.3f' %
             (mvars.component_name[0], mvars.component_name[1], distance))
#       print('distance = ', distance)

        pgui('radius of gyration of %s: %0.3f' %
             (mvars.component_name[0], radius_of_gyration_molecule1))
        pgui('radius of gyration of %s: %0.3f' %
             (mvars.component_name[1], radius_of_gyration_molecule2))
#       print('radius_of_gyration_molecule1: ', radius_of_gyration_molecule1)
#       print('radius_of_gyration_molecule2: ', radius_of_gyration_molecule2)

        rgcmdvars.outfile.write(str("%6i" % structure_number) + '\t\t' + str("%.3f" % distance) + '\t\t\t' +
                                str("%.3f" % radius_of_gyration_molecule1) + '\t\t\t' + str("%.3f" % radius_of_gyration_molecule2) + '\n')

        return

    def calculate_rg_cm_distance(self):
        '''
        **Calculate Rg CM Distance** reads the input PDB/DCD files, calculates the radius of gyration of each component in a two-component system as well as the distance between their center of mass, and writes the output file.

        Called from **Rg CM Distance Calculator** main program

        Calls **Do the Calculations**

        Parameters
        ----------

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
        output_file_path: string
            sub-path where output file will be written: run_name + \'rg_cm_distance_calculator'
        input_file_type: string
            type of trajectory file (PDB or DCD)
        output_file_name: string
            name of the output file (depends on the trajectory file name)
        outfile: string
            output file name (with full path): path + output_file_name
        python_basis: string array (dimension = number_of_components)
            basis string for each component in the molecule in python-readable format 

        Returns
        -------

        structure_number: int
            number of the current structure in the trajectory
        molecule: system.Molecule
            instance of the molecule class containing the coordinates of the current structure
        radius_of_gyration_molecule1: float
            radius of gyration of the first component
        radius_of_gyration_molecule2: float
            radius of gyration of the second component
        distance: float
            distance between the centers of mass of the two components

        '''

        # mvars used: pdb_file_name, trajectory_file_name, component_name, basis_string
        # rgcmdvars used: input_file_type, output_file_path, outfile, python_basis

        log = self.log
        log.debug('in calculate_rg_cm_distance')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        rgcmdvars = self.rg_cm_distance_calculator_variables

        molecule = system.Molecule(0)

        rgcmdvars.outfile.write('#\n')
        rgcmdvars.outfile.write('# Input PDB file: ' +
                                mvars.pdb_file_name+'\n')
        rgcmdvars.outfile.write('# Input DCD file: ' +
                                mvars.trajectory_file_name+'\n')
        rgcmdvars.outfile.write(
            '# Input basis string (' + mvars.component_name[0] + '): ' + mvars.basis_string[0]+'\n')
        rgcmdvars.outfile.write(
            '# Input basis string (' + mvars.component_name[1] + '): ' + mvars.basis_string[1]+'\n')
        rgcmdvars.outfile.write(
            '# Python basis string (' + mvars.component_name[0] + '): ' + rgcmdvars.python_basis[0]+'\n')
        rgcmdvars.outfile.write(
            '# Python basis string (' + mvars.component_name[1] + '): ' + rgcmdvars.python_basis[1]+'\n')
        rgcmdvars.outfile.write('#\n')
        rgcmdvars.outfile.write('#Structure ' + '\t' + 'CM Distance' +
                                '\t\t' + 'Rg (' + mvars.component_name[0] + ')\t\t\t' + 'Rg (' + mvars.component_name[1] + ')\n')
        rgcmdvars.outfile.write('#\n')

        if (rgcmdvars.input_file_type == 'pdb'):
            molecule.read_pdb(mvars.trajectory_file_name)
            natoms = molecule.natoms()
#           print('natoms = ', natoms)
            coor = numpy.zeros((1, natoms, 3), numpy.float32)
            number_of_frames_pdb = molecule.number_of_frames()
            pgui(">> input trajectory file is a pdb file")
#            pgui('number of pdb frames: %d' % (number_of_frames_pdb))
            pgui('\n There are %d files to process\n' %
                 (number_of_frames_pdb))
# create another instance of the molecule class to write the coordinates in each frame to a temporary PDB file
            molecule1 = system.Molecule(1)
            molecule1.read_pdb(mvars.trajectory_file_name)
# loop over files and calculate Rg and distance
            number_of_frames_processed = 0
# this is to mimic the behavior for reading a dcd file
            for i in range(number_of_frames_pdb):
                pgui('.\n')
            for i in range(number_of_frames_pdb):
                #                print('coords = ', molecule.coor()[i])
                number_of_frames_processed += 1
                structure_number = i+1
                log.debug('structure_number: %d' % (structure_number))
                temporary_file_name = rgcmdvars.output_file_path + 'temporary_file.pdb'
                pgui('extracting coordinates from frame: %d' %
                     (structure_number))
                coor[0, :, :] = molecule.coor()[i]
                molecule1.setCoor(coor)
                molecule1.write_pdb(temporary_file_name, 0, 'w')
                time.sleep(0.5)
#                print(temporary_file_name+'\n')
                molecule1.read_pdb(temporary_file_name)
                self.do_the_calculations(molecule1, structure_number)
                os.system('rm -f ' + temporary_file_name)
                fraction_done = (float(i+1)/float(number_of_frames_pdb))
                progress_string = 'COMPLETED ' + \
                    str(i+1)+' of '+str(number_of_frames_pdb)+' : ' + \
                    str(fraction_done*100.0)+' % done'
                pgui('\n%s \n' % progress_string)
            pgui('Processed %s pdb frames\n' % (number_of_frames_processed))

        elif (rgcmdvars.input_file_type == 'dcd'):
            molecule.read_pdb(mvars.pdb_file_name)
            read_trajectory_file = molecule.open_dcd_read(
                mvars.trajectory_file_name)
            number_of_frames_dcd = read_trajectory_file[2]
            pgui(">> input trajectory file is a dcd file")
#            print('number of dcd frames: %d' % (number_of_frames_dcd))
            pgui('\n There are %d files to process\n' %
                 (number_of_frames_dcd))
# loop over files and calculate Rg and distance
            number_of_frames_processed = 0
            for i in range(number_of_frames_dcd):
                number_of_frames_processed += 1
                structure_number = i+1
                log.debug('structure_number: %d' % (structure_number))
                temporary_file_name = rgcmdvars.output_file_path + 'temporary_file.pdb'
                pgui('extracting coordinates from frame: %d' %
                     (structure_number))
                molecule.read_dcd_step(read_trajectory_file, i)
                molecule.write_pdb(temporary_file_name, 0, 'w')
                time.sleep(0.5)
 #               print(temporary_file_name+'\n')
                molecule.read_pdb(temporary_file_name)
                self.do_the_calculations(molecule, structure_number)
                os.system('rm -f ' + temporary_file_name)
                fraction_done = (float(i+1)/float(number_of_frames_dcd))
                progress_string = 'COMPLETED ' + \
                    str(i+1)+' of '+str(number_of_frames_dcd)+' : ' + \
                    str(fraction_done*100.0)+' % done'
                pgui('\n%s \n' % progress_string)
            molecule.close_dcd_read(read_trajectory_file[0])
            pgui('Processed %s dcd frames\n' % (number_of_frames_processed))

        rgcmdvars.outfile.close()
        return

    def epilogue(self):
        '''Method to print out results and to move results to appropriate places.'''

        log = self.log
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        st = ''.join(['=' for x in range(60)])
        pgui('\n%s \n\n' % (st))
        pgui('RG CM DISTANCE CALCULATOR IS DONE')

        time.sleep(1.0)

        return
