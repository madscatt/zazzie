"""
    SASSIE: Copyright (C) 2011-2015 Joseph E. Curtis, Ph.D.

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
"""

from __future__ import division

import sys
import os
import random
import logging
import numpy
import string
import time
import subprocess

import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import sassie.util.basis_to_python as basis_to_python
import sassie.simulate.torsion_angle_monte_carlo.ooverlap as overlap
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.torsion as torsion
import sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.double_stranded_nucleic.double_stranded_nucleic_torsion as double_stranded_nucleic_torsion
import sassie.simulate.torsion_angle_monte_carlo.group_psf as group_psf

#       MONTE_CARLO
#
#       09/26/2005      --      gag-dihedral search             :       jc
#       07/14/2008      --      file management changes         :       jc
#       01/12/2011      --      added sasmol support            :       jc
#       01/15/2011      --      generalized rotation basis      :       jc
#       11/28/2014      --      group rotation basis            :       jc
#
# LC     1         2         3         4         5         6         7
# LC567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                     *      **
# $Id: monte_carlo.py 3040 2016-03-01 20:05:14Z schowell $

"""
    MONTE_CARLO is the module that contains the functions
    that are used to generate ensembles of structures by varying
    torision/dihedral angles.

    This module is called from Monte Carlo module from the main
    GUI through the pyqt_monte_carlo.py script.

    This module calls to C / Python extension modules to speed up
    calculations (see: overlap.c).

    REFERENCE: (for CHARMM force-field parameters)

    J. A. D. MacKerell et al.
    Journal of Physical Chemistry B,  102  3586-3616  (1998)

    B. R. Brooks et al.
    Journal of Computational Chemistry  4  187--217  (1983)

"""

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'monte_carlo'


class module_variables():

    def __init__(self, parent=None):
        self.app = app


class mcvars():

    def __init__(self, parent=None):
        pass


class simvars():

    def __init__(self, parent=None):
        pass


class simulation():

    def __init__(self, parent=None):
        pass


    def main(self, input_variables, txtOutput):
        """
        main method to manage simulation
        """

        self.mvars = module_variables()

        self.mcvars = mcvars()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.monte_carlo()

        self.epilogue()

        return


    def unpack_variables(self, variables):
        """
        method to extract variables into system wide class instance
        """

        mvars = self.mvars
        self.log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.dcdfile = variables['dcdfile'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.psffile = variables['psffile'][0]
        mvars.psf_flag = variables['psf_flag'][0]
        mvars.number_of_flexible_regions = variables[
            'number_of_flexible_regions'][0]
        mvars.basis_string_array = variables['basis_string_array'][0]
        mvars.delta_theta_array = variables['delta_theta_array'][0]
        mvars.rotation_type_array = variables['rotation_type_array'][0]
        mvars.rotation_direction_array = variables[
            'rotation_direction_array'][0]
        mvars.overlap_basis = variables['overlap_basis'][0]
        mvars.post_basis_string_array = variables['post_basis_string_array'][0]
        mvars.temperature = variables['temperature'][0]
        mvars.trial_steps = variables['trial_steps'][0]
        mvars.goback = variables['goback'][0]
        mvars.directed_mc = variables['directed_mc'][0]
        mvars.low_rg_cutoff = variables['low_rg_cutoff'][0]
        mvars.high_rg_cutoff = variables['high_rg_cutoff'][0]
        mvars.z_flag = variables['z_flag'][0]
        mvars.z_cutoff = variables['z_cutoff'][0]
        mvars.constraint_flag = variables['constraint_flag'][0]
        mvars.constraint_file = variables['constraint_file'][0]

        mvars.nonbondflag = variables['nonbondflag'][0]
        mvars.seed = variables['seed'][0]

        if mvars.psf_flag:

            mvars.max_steps = variables['max_steps'][0]
            mvars.energy_convergence = variables['energy_convergence'][0]
            mvars.step_size = variables['step_size'][0]

            mvars.toppath = os.path.join(sasconfig.__bin_path__, 'toppar')
            mvars.topfile = os.path.join(
                mvars.toppath, 'top_all27_prot_na.inp')
            mvars.parmfile = os.path.join(mvars.toppath,
                                          'par_all27_prot_na.inp')

        return


    def setup_group_simulation(self, group_pdb_name, group_psf_name):

        log = self.logprint_gui
        mvars = self.mvars

        log.debug('setting up group simulation object')

        mvars.energy_convergence = (mvars.energy_convergence *
                                    simtk.unit.kilojoule /
                                    simtk.unit.mole)  # the default value

        sim_variables = simvars()
        sim_variables.pdbfile = group_pdb_name
        sim_variables.psffile = group_psf_name

        sim_variables.topfile = mvars.topfile
        sim_variables.parmfile = mvars.parmfile

        sim_variables.integrator = simtk.openmm.LangevinIntegrator(
            mvars.temperature * simtk.unit.kelvin, 1 / simtk.unit.picosecond,
            mvars.step_size * simtk.unit.picoseconds)

        open_mm_simulation = open_mm.initialize_system(sim_variables)

        return open_mm_simulation


    def setup_groups(self):
        """
        method to create composite group molecules and collect group variables
        based on rotation sampling types
        """

        log = self.log
        mvars = self.mvars
        full_molecule = self.full_molecule
        pgui = self.run_utils.print_gui

        log.debug('in setup_groups')

        self.group_molecules = []
        self.group_masks = []
        self.group_flexible_molecules = []
        self.group_flexible_masks = []
        self.group_post_molecules = []
        self.group_post_masks = []
        self.group_variables = []

        self.group_simulation = []
        self.group_flexible_simulation = []
        self.group_post_simulation = []

        log.info('setting up %d groups' % mvars.number_of_flexible_regions)

        if mvars.psf_flag:
            mvars.psf_file_data = open(mvars.psffile).readlines()
            mvars.psf_data = group_psf.psf_data()
            group_psf.parse_psf_file(mvars.psf_file_data, mvars.psf_data)

        for group_number in xrange(mvars.number_of_flexible_regions):
            frame = 0
            this_group_rotation = mvars.rotation_type_array[group_number]

            """ now create the composite group molecule from input string """

            flexible_basis = mvars.basis_string_array_python[group_number]
            post_basis = mvars.post_basis_string_array_python[group_number]

            basis_string = '(' + flexible_basis + ') or (' + post_basis + ')'

            flexible_basis_string = '(' + flexible_basis + ')'
            post_basis_string = '(' + post_basis + ')'

            error, mask = full_molecule.get_subset_mask(basis_string)
            if DEBUG:
                assert not error, error

            this_group_molecule = sasmol.SasMol(0)
            this_group_flexible_molecule = sasmol.SasMol(1)
            this_group_post_molecule = sasmol.SasMol(2)

            error = full_molecule.copy_molecule_using_mask(this_group_molecule,
                                                           mask, frame)

            error, flexible_mask = this_group_molecule.get_subset_mask(
                flexible_basis_string)
            error, post_mask = this_group_molecule.get_subset_mask(
                post_basis_string)

            error = this_group_molecule.copy_molecule_using_mask(
                this_group_flexible_molecule, flexible_mask, frame)
            error = this_group_molecule.copy_molecule_using_mask(
                this_group_post_molecule, post_mask, frame)

            self.group_molecules.append(this_group_molecule)
            self.group_masks.append(mask)
            self.group_flexible_molecules.append(this_group_flexible_molecule)
            self.group_flexible_masks.append(flexible_mask)
            self.group_post_molecules.append(this_group_post_molecule)
            self.group_post_masks.append(post_mask)

            if mvars.psf_flag:

                group_pdb_name = os.path.join(self.runpath, "group_%d.pdb"
                                              % group_number)
                group_psf_name = os.path.join(self.runpath, "group_%d.psf"
                                              % group_number)

                group_flexible_pdb_name = os.path.join(self.runpath,
                                                       "group_flexible_%d.pdb"
                                                       % group_number)
                group_flexible_psf_name = os.path.join(self.runpath,
                                                       "group_flexible_%d.psf"
                                                       % group_number)

                group_post_pdb_name = os.path.join(self.runpath,
                                                   "group_post_%d.pdb"
                                                   % group_number)
                group_post_psf_name = os.path.join(self.runpath,
                                                   "group_post_%d.psf"
                                                   % group_number)

                group_psf.get_group_psf(self, group_psf_name, mvars.psf_data,
                                        this_group_molecule)
                group_psf.get_group_psf(self, group_flexible_psf_name,
                                        mvars.psf_data,
                                        this_group_flexible_molecule)
                group_psf.get_group_psf(self, group_post_psf_name,
                                        mvars.psf_data,
                                        this_group_post_molecule)

                index = range(1, this_group_molecule.natoms() + 1)
                this_group_molecule.setIndex(index)

                index = range(1, this_group_flexible_molecule.natoms() + 1)
                this_group_flexible_molecule.setIndex(index)

                index = range(1, this_group_post_molecule.natoms() + 1)
                this_group_post_molecule.setIndex(index)

                this_group_molecule.write_pdb(group_pdb_name, frame, "w")
                this_group_flexible_molecule.write_pdb(group_flexible_pdb_name,
                                                       frame, "w")
                this_group_post_molecule.write_pdb(group_post_pdb_name, frame,
                                                   "w")

                """
                # we removed the openmm bits
                if sasconfig.__openmm__:
                    group_simulation_object = self.setup_group_simulation(
                        group_pdb_name, group_psf_name)
                    group_flexible_simulation_object = self.setup_group_simulation(
                        group_flexible_pdb_name, group_flexible_psf_name)
                    group_post_simulation_object = self.setup_group_simulation(
                        group_post_pdb_name, group_post_psf_name)

                    self.group_simulation.append(group_simulation_object)
                    self.group_flexible_simulation.append(
                        group_flexible_simulation_object)
                    self.group_post_simulation.append(group_post_simulation_object)
                """
            else:
                log.error("ERROR: psf is needed")

            if this_group_rotation == 'double_stranded_nucleic_torsion':
                log.debug('setting up double stranded nucleic nucleic '
                          'rotation group %d' % group_number)

                double_stranded_nucleic_torsion.setup(self, group_number)

                self.group_variables.append(mvars.nvars)


            else:
                log.debug('setting up '+this_group_rotation)
                torsion.setup(self, group_number)

                log.debug('assigning group %d variables from %s' % (group_number, this_group_rotation))

                self.group_variables.append(mvars.pvars)


        # set the overlap basis mask for the full molecule
        overlap_basis = mvars.overlap_basis
        log.debug('creating overlap basis for the full molecule using %s' %
                  overlap_basis)
        if overlap_basis.lower() == 'all':
            overlap_basis_filter = 'not name[i] == None'
            mvars.cutoff = 0.8
        elif overlap_basis.lower() == 'backbone':
            protein_backbone_filter = ('(moltype[i] == "protein" and ('
                                       'name[i] == "N" or name[i] == "CA" or '
                                       'name[i] == "C"))')
            dna_backbone_filter = ('(moltype[i] == "dna" and ('
                                   'name[i] == "O2P" or name[i] == "O1P" or '
                                   'name[i] == "P" or name[i] == "O5\'" or '
                                   'name[i] == "C5\'" or '
                                   'name[i] == "C4\'" or '
                                   'name[i] == "C3\'" or name[i] == "O3\'"))')
            overlap_basis_filter = ('%s or %s' % (protein_backbone_filter,
                                                  dna_backbone_filter))
            mvars.cutoff = 0.8
        elif overlap_basis.lower() == 'heavy':
            overlap_basis_filter = 'not name[i][0] == "H" '
            mvars.cutoff = 0.8
        log.info('overlap_basis_filter = ' + overlap_basis_filter)
        error, self.full_overlap_basis_mask = full_molecule.get_subset_mask(
            overlap_basis_filter)

        return


    def setup_constraints(self):
        """
        method to setup basis masks for filtering structures based on
        selections in the supplied constraint text file
        """

        log = self.log
        mvars = self.mvars
        full_molecule = self.full_molecule

        mvars.filter_flag = 0

        (error, constraint_basis1_array, constraint_basis2_array,
         mvars.distance_array, mvars.type_array
         ) = constraints.read_constraints(full_molecule, mvars.constraint_file,
                                          mvars.filter_flag)

        mvars.mask_a_array = []
        mvars.mask_b_array = []

        for i in xrange(len(distance_array)):
            log.debug('constraint_basis1_array[i] = ' +
                      str(constraint_basis1_array[i]))
            log.debug('constraint_basis2_array[i] = ' +
                      str(constraint_basis2_array[i]))
            log.debug('distance_array[i] = ' + str(distance_array[i]))
            log.debug('type_array[i] = ' + str(type_array[i]))

            error, local_mask_a_array = full_molecule.get_subset_mask(
                constraint_basis1_array[i])
            error, local_mask_b_array = full_molecule.get_subset_mask(
                constraint_basis2_array[i])

            mvars.mask_a_array.append(local_mask_a_array)
            mvars.mask_b_array.append(local_mask_b_array)

        return


    def get_bond_list(self):
        """
        method to get list of atoms that are bonded so that they are
        not used for overlap check

        """

        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui
        full_molecule = self.full_molecule

        mvars.number_of_bonds = 0
        mvars.bonds = []

        for line in mvars.psf_data.nbond:
            this_line = string.split(line)
            if "!NBOND:" not in this_line and len(this_line) > 1:
                try:
                    local_number_of_pairs = int(len(this_line) / 2)
                    j = 0
                    for i in xrange(local_number_of_pairs):
                        mvars.bonds.append([int(this_line[j]) - 1,
                                            int(this_line[j + 1]) - 1])
                        mvars.number_of_bonds += 1
                    j += 2
                except:
                    log.error("ERROR: bond read error")
                    log.error("ERROR: line = " + str(line))
                    log.error("ERROR: len(this_line)/2 = %d" %
                              len(this_line) / 2)
            elif "!NBOND:" in this_line:
                psf_number_of_bonds = int(this_line[0])

        if psf_number_of_bonds != mvars.number_of_bonds:
            log.error("ERROR: PSF file indicates that there are %d bonds" %
                      psf_number_of_bonds)
            log.error("ERROR: get_bond_list() method found %d bonds" %
                      mvars.number_of_bonds)

        log.debug('found %d bonds in complete molecule\n' % mvars.number_of_bonds)

        return


    def initialization(self):
        """
        method to prepare for simulation
        """

        log = self.log
        log.debug('in initialization')
        mvars = self.mvars

        """
        convert input basis string to enable python version
        that is compatiable with eval
        """

        basis_string_array = mvars.basis_string_array
        basis_string_array_python = []

        for basis in basis_string_array:
            basis_string_array_python.append(
                basis_to_python.parse_basis(basis))
        mvars.basis_string_array_python = basis_string_array_python

        post_basis_string_array = mvars.post_basis_string_array
        post_basis_string_array_python = []

        for basis in post_basis_string_array:
            post_basis_string_array_python.append(
                basis_to_python.parse_basis(basis))
        mvars.post_basis_string_array_python = post_basis_string_array_python

        """ the units are J/K """
        kb = 1.380658E-23
        mvars.beta = 1.0 / (mvars.temperature * kb)

        self.full_molecule = sasmol.SasMol(0)

        log.info('opening pdb file: ' + mvars.pdbfile)
        # self.full_molecule.read_pdb(mvars.pdbfile, saspdbrx_topology=True) #SH: TEMPORARY
        self.full_molecule.read_pdb(mvars.pdbfile) #SH: TEMPORARY

        log.info('opening new dcd file to store trajectory: %s' %
                 os.path.join(self.runpath, mvars.dcdfile))
        self.dcdoutfile = self.full_molecule.open_dcd_write(
            os.path.join(self.runpath, mvars.dcdfile))

        cmd = "cp %s %s %s" % (mvars.psffile, mvars.pdbfile, self.runpath)
        return_code = subprocess.call(cmd, shell=True)
        if return_code:
            log.info('FAILED to copy pdb and psf file to run directory')
        else:
            log.info('copied pdb and psf file to run directory')

        self.full_molecule.set_average_vdw()

        new_vdw = []
        for vdw in self.full_molecule.atom_vdw():
            if len(vdw) > 0:
                new_vdw.append(vdw[-1])
            else:
                new_vdw.append(vdw)

        self.full_molecule.setAtom_vdw(new_vdw)

        self.setup_groups()

        if mvars.constraint_flag:
            self.setup_constraints()

        else:
            mvars.mask_a_array = []
            mvars.mask_b_array = []
            mvars.distance_array = []
            mvars.type_array = []

        if mvars.psf_flag:
            self.get_bond_list()

        return


    def evaluate_rg(rg_difference_list, directed_rg_list, accepted_rg_list,
                    this_rg_difference, this_rg, accepted):

        maximum_value = max(rg_difference_list)

        if maximum_value > this_rg_difference:
            index = rg_difference_list.index(maximum_value)
            rg_difference_list[index] = this_rg_difference
            directed_rg_list[index] = this_rg
            accepted_rg_list[index] = accepted

        return


    def update_full_molecule_coor(self, group_number):
        """
        method to update the full molecule coordinates from a certain group
        """

        log = self.log
        mvars = self.mvars
        full_molecule = self.full_molecule

        group_molecule = self.group_molecules[group_number]
        group_mask = self.group_masks[group_number]

        error = full_molecule.set_coor_using_mask(
            group_molecule, 0, group_mask)

        return


    def update_group_coor(self, group_number):
        """
        method to update a given group coordinate from the full molecule
        """

        log = self.log
        mvars = self.mvars
        full_molecule = self.full_molecule

        error, group_coordinates = full_molecule.get_coor_using_mask(
            0, self.group_masks[group_number])
        self.group_molecules[group_number].setCoor(group_coordinates)

        return


    def update_all_group_coor(self):
        """
        method to update all the group coordinates from the full molecule
        """

        log = self.log
        mvars = self.mvars
        full_molecule = self.full_molecule

        for n in xrange(mvars.number_of_flexible_regions):
            error, group_coordinates = full_molecule.get_coor_using_mask(
                0, self.group_masks[n])
            self.group_molecules[n].setCoor(group_coordinates)

        return


    def check_group_overlap(self, group):
        """
        checks all atom overlap over entire molecule
        """

        log = self.log
        mvars = self.mvars
        mcvars = self.mcvars
        pgui = self.run_utils.print_gui
        group_molecule = self.group_molecules[group]
        group_variable = self.group_variables[group]
        overlap_basis_mask = group_variable.overlap_basis_mask
        cutoff = mvars.cutoff

        error, basis_coor = group_molecule.get_coor_using_mask(
            0, overlap_basis_mask)

        flag = overlap.ooverlap(basis_coor[0], float(cutoff))

        if flag == 1:
            mcvars.trial_accepted = False
            mcvars.number_overlap_failures += 1
            log.debug('~~~~ trial step failed: overlap ~~~~')
        else:
            mcvars.trial_accepted = True

        return flag


    def check_all_overlap(self):
        """
        checks all atom overlap over entire molecule
        """
        log = self.log
        mvars = self.mvars
        mcvars = self.mcvars
        pgui = self.run_utils.print_gui
        full_molecule = self.full_molecule
        full_overlap_basis_mask = self.full_overlap_basis_mask
        cutoff = mvars.cutoff

        error, basis_coor = full_molecule.get_coor_using_mask(
            0, full_overlap_basis_mask)

        flag = overlap.ooverlap(basis_coor[0], float(cutoff))

        if flag == 1:
            mcvars.trial_accepted = False
            mcvars.number_overlap_failures += 1
            log.debug('~~~~ trial step failed: overlap ~~~~')
        else:
            mcvars.trial_accepted = True

        return


    def monte_carlo(self):
        """
        method to carry out Monte Carlo simulation
        """

        log = self.log
        mvars = self.mvars
        mcvars = self.mcvars
        pgui = self.run_utils.print_gui
        full_molecule = self.full_molecule

        # start gui output
        pgui("\n%s \n" % ('=' * 60))
        pgui("DATA FROM RUN: %s \n\n" % time.asctime( time.gmtime( time.time() ) ))

        frame = 0

        log.debug('in monte_carlo')
        pgui('>>>> running Monte Carlo simulation for %d steps <<<<\n\n' %
             mvars.trial_steps)

        """ set up random seed """

        if mvars.seed[0] == 1:
            from numpy.random import RandomState
            mvars.seed_object = RandomState(mvars.seed[1])
        else:
            mvars.seed_object = -1

        """ main loop """

        if mvars.directed_mc > 0:
            rg_difference_list = []
            directed_rg_list = []
            accepted_rg_list = []
            rg_list_length = 10  # hardwired

        mcvars.trial_accepted = False
        mcvars.accepted = 0
        mcvars.number_overlap_failures = 0
        mcvars.number_constraint_failures = 0
        mcvars.number_rg_failures = 0
        mcvars.number_z_failures = 0
        mcvars.fail_tally = 0

        for i in xrange(mvars.trial_steps):

            """
            pick a group
            """
            frame = 0

            if DEBUG:
                print i,  # comma prevent printing a line return after 'i'
                if i == mvars.trial_steps - 1:
                    print
                sys.stdout.flush()

            group_number = random.randint(
                0, mvars.number_of_flexible_regions - 1)

            self.update_group_coor(group_number)

            this_group_rotation = mvars.rotation_type_array[group_number]

            if this_group_rotation == 'double_stranded_nucleic_torsion':
                double_stranded_nucleic_torsion.rotate_coarse_grained_beads(
                    self, group_number)
            else:
                torsion.sample(self, group_number)

            self.update_full_molecule_coor(group_number)

            # verify the new structure meets all requiremnts
            if mcvars.trial_accepted:
                self.check_all_overlap()

            if mvars.constraint_flag and mcvars.trial_accepted:
                if 1 == constraints.check_constraints(m1, mask_a_array,
                                                      mask_b_array,
                                                      distance_array,
                                                      type_array):
                    mcvars.trial_accepted = False
                    mcvars.number_constraint_failures += 1
                    log.debug(' ~~~~ trial step failed: constraints ~~~~')

            if mvars.z_flag and mcvars.trial_accepted:
                if not numpy.alltrue(full_molecule.coor()[0, :, 2] >=
                                     mvars.z_cutoff):
                    mcvars.trial_accepted = False
                    mcvars.number_z_failures += 1
                    log.debug(' ~~~~ trial step failed: z_flag ~~~~')

            if mcvars.trial_accepted:
                current_rg = full_molecule.calcrg(frame)
                if not (current_rg > mvars.low_rg_cutoff and
                        current_rg < mvars.high_rg_cutoff):
                    mcvars.trial_accepted = False
                    mcvars.number_rg_failures += 1
                    log.debug(' ~~~~ trial step failed: rg_cutoff ~~~~')

            # process the new structure (store or discard)
            if mcvars.trial_accepted:
                mcvars.accepted += 1
                mcvars.fail_tally = 0
                log.debug(' ~~~~ trial step accepted ~~~~')
                full_molecule.write_dcd_step(self.dcdoutfile, 0,
                                             mcvars.accepted)

                if mvars.directed_mc > 0:
                    if len(rg_difference_list) <= rg_list_length:
                        this_rg_difference = abs(rg_value - directed_mc)
                        rg_difference_list.append(this_rg_difference)
                        directed_rg_list.append(rg_value)
                        accepted_rg_list.append(accepted)
                    else:
                        this_rg_difference = abs(rg_value - directed_mc)
                        evaluate_rg(rg_difference_list, directed_rg_list,
                                    accepted_rg_list, this_rg_difference,
                                    rg_value, accepted)

                # ds_dna track accepted steps
                if this_group_rotation == 'double_stranded_nucleic_torsion':
                    double_stranded_nucleic_torsion.post_process(self,
                                                                 group_number)
                    self.update_full_molecule_coor(group_number)

            else:
                mcvars.fail_tally += 1

                if (mcvars.fail_tally == mvars.goback):

                    mcvars.fail_tally = 0

                    if mcvars.accepted > 0:

                        if mvars.goback == 1:
                            goback_frame = mcvars.accepted

                        elif mvars.seed[0] == 1:
                            ran_num = mvars.seed_object.rand()
                            goback_frame = int(mcvars.accepted * ran_num)

                        elif mvars.directed_mc > 0:
                            local_rg_list_length = len(directed_rg_list)
                            ran_num = random.randrange(0, local_rg_list_length)
                            goback_frame = accepted_rg_list[ran_num] + 1

                        else:
                            goback_frame = random.randint(0, mcvars.accepted)

                    else:
                        goback_frame = 0

                    if goback_frame == 0:
                        log.debug('\nreloading coordinates from original '
                                  'starting structure')
                        full_molecule.read_pdb(mvars.pdbfile, fastread=True,
                                               #saspdbrx_topology=True)
                                               saspdbrx_topology=False)

                    else:
                        log.debug('\nreloading coordinates from a previously '
                                  'accepted structure')
                        full_molecule.read_single_dcd_step(os.path.join(
                            self.runpath, mvars.dcdfile), goback_frame)

                    if ('double_stranded_nucleic_torsion' in
                        mvars.rotation_type_array):
                        double_stranded_nucleic_torsion.step_counter_goback(
                            self, goback_frame)

            # display progress
            fraction_done = (i + 1) * 1.0 / mvars.trial_steps
            report_string = 'STATUS\t%f' % fraction_done
            pgui(report_string)

        full_molecule.close_dcd_write(self.dcdoutfile)

        return


    def epilogue(self):
        """
        method to print out simulation results and to move results
        to appropriate places.
        """

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        if ('double_stranded_nucleic_torsion' in self.mvars.rotation_type_array
            and DEBUG):
            double_stranded_nucleic_torsion.epilogue(self)

        pgui('Configurations and statistics saved in %s directory\n\n' % (self.runpath+os.path.sep))

        pgui('DCD data were written to %s\n\n' % os.path.join(self.runpath, self.mvars.dcdfile))

        fraction_accepted = self.mcvars.accepted / self.mvars.trial_steps
        pgui('accepted %d out of %d : %0.2f%% percent\n' %
             (self.mcvars.accepted, self.mvars.trial_steps, 100.0 * fraction_accepted))
        pgui('overlapped %d out of %d : %0.2f%% percent\n' %
             (self.mcvars.number_overlap_failures, self.mvars.trial_steps, 100.0 *
              self.mcvars.number_overlap_failures / self.mvars.trial_steps))
        self.run_utils.clean_up(log)

        pgui('\nrun json inputs saved to:\n    %s\n' %
             os.path.join(self.runpath, self.parmfile))
        pgui('\nrun log output saved to:\n    %s\n' %
             os.path.join(self.runpath, self.logfile))

        pgui("\n\n")
        pgui("%s \n" % ('=' * 60))
        time.sleep(2)

        return
