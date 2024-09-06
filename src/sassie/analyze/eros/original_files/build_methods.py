from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

"""
    SASSIE: Copyright (C) 2011-2016 Joseph E. Curtis, Ph.D.

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

import sys
import os
import random
import logging
import numpy
import math
import string
import time
import copy

import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import sassie.util.basis_to_python as basis_to_python
import sasmol.sasmath as sasmath
import sasmol.sasutil as sasutil

#
#       BUILD_METHODS 
#
#       09/04/2016      --      initial coding                  :       jc
#
# LC     1         2         3         4         5         6         7
# LC567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                     *      **

"""
    BUILD_METHODS is the module that contains the functions that are
    used to alter and create pdb files for common model building tasks.


"""

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'build_utilities'


class module_variables():

    def __init__(self, parent=None):
        self.app = app

class build_utilities():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):
        """
        main method to manage module
        """

        self.mvars = module_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.pdb_and_fasta()

        self.epilogue()

        return


    def unpack_variables(self, variables):
        """
        method to extract variables into system wide class instance
        """

        mvars = self.mvars
        self.log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        
        mvars.pdb_utilities_flag = variables['pdb_utilities_flag'][0]
        mvars.fasta_utilities_flag = variables['fasta_utilities_flag'][0]
         
        if(mvars.pdb_utilities_flag):
            mvars.input_pdbfile = variables['input_pdbfile'][0]
       
            mvars.renumber_flag = variables['renumber_flag'][0]
            mvars.pdb_constraints_flag = variables['pdb_constraints_flag'][0]
            mvars.modify_fields_flag = variables['modify_fields_flag'][0]
            mvars.translation_rotation_flag = variables['translation_rotation_flag'][0]
            mvars.align_pmi_on_cardinal_axes_flag = variables['align_pmi_on_cardinal_axes_flag'][0]
            mvars.align_pmi_on_axis_flag = variables['align_pmi_on_axis_flag'][0]

            if(mvars.renumber_flag): 
                mvars.renumber_output_filename = variables['renumber_output_filename'][0]
                mvars.renumber_indices_flag = variables['renumber_indices_flag'][0]
        
                if(mvars.renumber_indices_flag):
                    mvars.first_index = variables['first_index'][0]

                mvars.renumber_resids_flag = variables['renumber_resids_flag'][0]
                
                if(mvars.renumber_resids_flag):
                    mvars.first_resid = variables['first_resid'][0]
           
            elif(mvars.pdb_constraints_flag):
                mvars.number_of_constraint_files = variables['number_of_constraint_files'][0]
                mvars.constraint_options = variables['constraint_options'][0]
                mvars.constraint_filenames = variables['constraint_filenames'][0]
                mvars.constraint_fields = variables['constraint_fields'][0]
                mvars.constraint_resets = variables['constraint_resets'][0]
            
            elif(mvars.modify_fields_flag):
                mvars.modify_fields_output_filename = variables['modify_fields_output_filename'][0]
                mvars.number_of_fields_to_modify = variables['number_of_fields_to_modify'][0]
                mvars.field_selections = variables['field_selections'][0]
                mvars.field_options = variables['field_options'][0]
                mvars.field_values = variables['field_values'][0]

            elif(mvars.translation_rotation_flag):
                mvars.translation_rotation_output_filename = variables['translation_rotation_output_filename'][0]
                mvars.pre_center_flag = variables['pre_center_flag'][0] 
                mvars.translation_array = variables['translation_array'][0] 
                mvars.rotation_flag = variables['rotation_flag'][0] 
                mvars.rotation_type = variables['rotation_type'][0] 
                mvars.rotation_axes_order = variables['rotation_axes_order'][0] 
                mvars.rotation_theta = variables['rotation_theta'][0] 
                mvars.user_vector_1 = variables['user_vector_1'][0]
           
            elif(mvars.align_pmi_on_axis_flag):
                mvars.align_pmi_output_filename = variables['align_pmi_output_filename'][0]
                mvars.pmi_eigenvector = variables['pmi_eigenvector'][0]
                mvars.alignment_vector_axis = variables['alignment_vector_axis'][0]
                mvars.user_vector_2 = variables['user_vector_2'][0]
                mvars.settle_on_surface_flag = variables['settle_on_surface_flag'][0]
                if(mvars.settle_on_surface_flag):
                    mvars.surface_plane = variables['surface_plane'][0]
                    mvars.invert_along_axis_flag = variables['invert_along_axis_flag'][0]
            elif(mvars.align_pmi_on_cardinal_axes_flag):
                mvars.align_pmi_output_filename = variables['align_pmi_output_filename'][0]

        elif(mvars.fasta_utilities_flag):
            mvars.fasta_input_option = variables['fasta_input_option'][0]
            mvars.fasta_output_filename  = variables['fasta_output_filename'][0]

            if(mvars.fasta_input_option == 'sequence'):
                mvars.fasta_input_sequence = variables['fasta_input_sequence'][0] 
            elif(mvars.fasta_input_option == 'file'):
                mvars.fasta_input_file = variables['fasta_input_file'][0] 
            
            mvars.fasta_moltype = variables['fasta_moltype'][0]

        mvars.seed = variables['seed'][0]

        return

    def initialization(self):
        """
        method to prepare for data for main method
        """

        log = self.log
        log.debug('in initialization')
        mvars = self.mvars

        if(mvars.pdb_utilities_flag):
            mvars.molecule = sasmol.SasMol(0)
            mvars.molecule.read_pdb(mvars.input_pdbfile)

            if(mvars.pdb_constraints_flag):

                mvars.constraint_options = [x.strip() for x in mvars.constraint_options.split(',')]
                mvars.constraint_fields = [x.strip() for x in mvars.constraint_fields.split(',')]
                mvars.constraint_filenames = [x.strip() for x in mvars.constraint_filenames.split(',')]
                mvars.constraint_resets = [x.strip() for x in mvars.constraint_resets.split(',')]
                
                dum = [] 
                for value in mvars.constraint_resets:
                    dum.append(eval(value))
                mvars.constraint_resets = dum

            if(mvars.modify_fields_flag):

                mvars.field_selections = [x.strip() for x in mvars.field_selections.split(',')]
                mvars.field_options = [x.strip() for x in mvars.field_options.split(',')]
                mvars.field_values = [x.strip() for x in mvars.field_values.split(',')]

                mvars.field_selections_python = []

                for basis in mvars.field_selections:
                    mvars.field_selections_python.append(basis_to_python.parse_basis(basis))

            if(mvars.translation_rotation_flag):

                mvars.translation_array = numpy.array(mvars.translation_array, numpy.float)
                print('translation array = ', mvars.translation_array) 
                if(mvars.rotation_flag):
                    dum_order = []
                    dum_theta = {}
                    for axis in mvars.rotation_axes_order:
                        dum_order.append(axis)
                        if axis == 'x':
                            dum_theta['x'] = mvars.rotation_theta[0] * math.pi/180.0
                        elif axis == 'y':
                            dum_theta['y'] = mvars.rotation_theta[1] * math.pi/180.0
                        elif axis == 'z':
                            dum_theta['z'] = mvars.rotation_theta[2] * math.pi/180.0
                    mvars.rotation_axes_order = dum_order
                    mvars.rotation_theta = dum_theta
                    print('rt = ', mvars.rotation_theta)
                    if(mvars.rotation_type == 'user_vector'):
                        mvars.user_vector_1 = numpy.array(mvars.user_vector_1, numpy.float)
                        print('user_vector_1 = ', mvars.user_vector_1) 

        if mvars.fasta_utilities_flag:

            if(mvars.fasta_input_option == 'file'):
                with open(mvars.fasta_input_file) as fasta_input:
                    mvars.fasta_input_sequence = fasta_input.read().splitlines() 
                all_sequences = sasutil.parse_fasta(mvars.fasta_input_sequence)
                mvars.fasta_input_sequence = all_sequences[0]

    def pdb_and_fasta(self):
        """
        main method of module
        """

        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui

        # start gui output
        pgui("\n%s \n" % ('=' * 60))
        pgui("DATA FROM RUN: %s \n\n" % time.asctime( time.gmtime( time.time() ) ))

        frame = 0

        log.debug('in build_utilities')

        """ set up random seed """

        if mvars.seed[0] == 1:
            from numpy.random import RandomState
            mvars.seed_object = RandomState(mvars.seed[1])
        else:
            mvars.seed_object = -1

        """ main loop """

        if mvars.pdb_utilities_flag:

            if(mvars.renumber_flag):
                if(mvars.renumber_indices_flag):
                    if(mvars.renumber_resids_flag):
                        mvars.molecule.renumber(index=mvars.first_index, resid=mvars.first_resid)
                    else:
                        mvars.molecule.renumber(index=mvars.first_index)
                elif(mvars.renumber_resids_flag):
                    mvars.molecule.renumber(resid=mvars.first_resid)

                mvars.molecule.write_pdb(os.path.join(self.runpath, mvars.renumber_output_filename), frame, 'w')
                mvars.output_pdb_filename = [mvars.renumber_output_filename]

            elif(mvars.pdb_constraints_flag):
                mvars.output_pdb_filename = []
                for i in xrange(mvars.number_of_constraint_files):
                    mvars.molecule.make_constraint_pdb(os.path.join(self.runpath, mvars.constraint_filenames[i]), mvars.constraint_options[i], field=mvars.constraint_fields[i], reset=mvars.constraint_resets[i])
                    mvars.output_pdb_filename.append(mvars.constraint_filenames[i])

            elif(mvars.modify_fields_flag):
                frame = 0 
                
                for i in xrange(mvars.number_of_fields_to_modify):
                    field_selection = mvars.field_selections_python[i]  
                    field_option = mvars.field_options[i]  
                    field_value = mvars.field_values[i] 

                    error, mask = mvars.molecule.get_subset_mask(field_selection)
                    #if DEBUG:
                    #    assert not error, error

                    if field_option == 'record':
                        descriptor = mvars.molecule.atom()
                        error=mvars.molecule.set_descriptor_using_mask(mask,descriptor,field_value)
                    elif field_option == 'altloc':
                        descriptor = mvars.molecule.loc()
                        error=mvars.molecule.set_descriptor_using_mask(mask,descriptor,field_value)
                    elif field_option == 'icode':
                        descriptor = mvars.molecule.rescode()
                        error=mvars.molecule.set_descriptor_using_mask(mask,descriptor,field_value)
                    elif field_option == 'xcoor':
                        for j in xrange(mvars.molecule.natoms()):
                            if mask[j] == 1:
                                mvars.molecule.coor()[frame, j, 0] = field_value
                    elif field_option == 'ycoor':
                        for j in xrange(mvars.molecule.natoms()):
                            if mask[j] == 1:
                                mvars.molecule.coor()[frame, j, 1] = field_value
                    elif field_option == 'zcoor':
                        for j in xrange(mvars.molecule.natoms()):
                            if mask[j] == 1:
                                mvars.molecule.coor()[frame, j, 2] = field_value
                    else:
                        descriptor = eval('mvars.molecule.'+field_option+'()')
                        error=mvars.molecule.set_descriptor_using_mask(mask,descriptor,field_value)

                    #if DEBUG:
                    #    assert not error, error

                mvars.molecule.write_pdb(os.path.join(self.runpath, mvars.modify_fields_output_filename), frame, 'w')
                mvars.output_pdb_filename = [mvars.modify_fields_output_filename]

            elif(mvars.translation_rotation_flag):

                if(mvars.pre_center_flag):
                    mvars.molecule.center(frame)

                mvars.molecule.translate(frame, mvars.translation_array)            

                if(mvars.rotation_flag):
                    if(mvars.rotation_type == 'cardinal'):
                        for axis in mvars.rotation_axes_order:
                            center_of_mass = mvars.molecule.calccom(frame)
                            mvars.molecule.center(frame)
                            theta = mvars.rotation_theta[axis]
                            mvars.molecule.rotate(frame,axis,theta)
                            mvars.molecule.moveto(frame, center_of_mass)
                    elif(mvars.rotation_type == 'user_vector'):
                        for axis in mvars.rotation_axes_order:
                            center_of_mass = mvars.molecule.calccom(frame)
                            mvars.molecule.center(frame)
                            theta = mvars.rotation_theta[axis]
                            ux = mvars.user_vector_1[0]
                            uy = mvars.user_vector_1[1]
                            uz = mvars.user_vector_1[2]
                            mvars.molecule.general_axis_rotate(frame,theta,ux,uy,uz)
                            mvars.molecule.moveto(frame, center_of_mass)
        
                mvars.molecule.write_pdb(os.path.join(self.runpath, mvars.translation_rotation_output_filename), frame, 'w')
                mvars.output_pdb_filename = [mvars.translation_rotation_output_filename]

            elif(mvars.align_pmi_on_cardinal_axes_flag):
                frame = 0
                mvars.molecule.align_pmi_on_cardinal_axes(frame) 
                mvars.molecule.write_pdb(os.path.join(self.runpath,mvars.align_pmi_output_filename), frame, 'w')
                mvars.output_pdb_filename = [mvars.align_pmi_output_filename]

            elif(mvars.align_pmi_on_axis_flag):
                frame = 0
                mvars.molecule.align_pmi_on_axis(frame, mvars.pmi_eigenvector,\
                                  mvars.alignment_vector_axis)

                if(mvars.settle_on_surface_flag):
                    minmax = mvars.molecule.calcminmax()
                    #print('minmax = ', minmax)
                    minimum = minmax[0]
                    maximum = minmax[1]
                    #print('minimum = ', minimum)
                    #print('maximum = ', maximum)
                    if(mvars.alignment_vector_axis == 'x'):
                        mvars.molecule.translate(frame, [-minimum[0], 0.0, 0.0])
                    elif(mvars.alignment_vector_axis == 'y'):
                        mvars.molecule.translate(frame, [0.0, -minimum[1], 0.0])
                    elif(mvars.alignment_vector_axis == 'z'):
                        mvars.molecule.translate(frame, [0.0, 0.0, -minimum[2]])
                    minmax = mvars.molecule.calcminmax()
                    #print('minmax = ', minmax)

                    #mvars.molecule.write_pdb(os.path.join(self.runpath,'moved_'+mvars.align_pmi_output_filename), frame, 'w')
                    if(mvars.invert_along_axis_flag):
                        center_of_mass = mvars.molecule.calccom(frame)
                        mvars.molecule.center(frame)
                        if(mvars.alignment_vector_axis == 'x'):
                            mvars.molecule.rotate(frame, 'y', math.pi)
                        elif(mvars.alignment_vector_axis == 'y'):
                            mvars.molecule.rotate(frame, 'z', math.pi)
                        elif(mvars.alignment_vector_axis == 'z'):
                            mvars.molecule.rotate(frame, 'x', math.pi)
                        mvars.molecule.moveto(frame, center_of_mass)
                        #mvars.molecule.write_pdb(os.path.join(self.runpath,'inverted_moved_'+mvars.align_pmi_output_filename), frame, 'w')

                mvars.molecule.write_pdb(os.path.join(self.runpath,mvars.align_pmi_output_filename), frame, 'w')
                mvars.output_pdb_filename = [mvars.align_pmi_output_filename]
           
        elif mvars.fasta_utilities_flag:
            molecule = sasmol.SasMol(0)
            molecule.setFasta(mvars.fasta_input_sequence)
            molecule.make_backbone_pdb_from_fasta(os.path.join(self.runpath, mvars.fasta_output_filename), mvars.fasta_moltype)
            mvars.output_pdb_filename = [mvars.fasta_output_filename]
 
        return
    

    def epilogue(self):
        """
        method to print out and move results to appropriate places
        """

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui
        mvars = self.mvars

        fraction_done = 1.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)

        pgui('New PDB file(s) saved in %s directory\n' % (self.runpath+os.path.sep)) 
        for name in mvars.output_pdb_filename:
            pgui(name+'\n')
        self.run_utils.clean_up(log) 
        pgui('\nrun json inputs saved to:\n    %s\n' % os.path.join(self.runpath, self.parmfile)) 
        pgui('\nrun log output saved to:\n    %s\n' % os.path.join(self.runpath, self.logfile)) 
        pgui("\n\n") 
        pgui("%s \n" % ('=' * 60)) 
        
        time.sleep(2) 

        return
