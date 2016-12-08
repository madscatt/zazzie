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

import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import sassie.util.basis_to_python as basis_to_python
import sasmol.sasmath as sasmath

#
#       TOOL_METHODS 
#
#       09/04/2016      --      initial coding                  :       jc
#
# LC     1         2         3         4         5         6         7
# LC567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                     *      **

"""
    TOOL_METHODS is the module that contains the functions that are
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
         
        if(mvars.pdb_utilities_flag):
            mvars.pdbfile = variables['pdbfile'][0]
       
        mvars.renumber_flag = variables['renumber_flag'][0]
        
        if(mvars.renumber_flag): 
            mvars.renumber_output_filename = variables['renumber_output_filename'][0]
            mvars.renumber_indices_flag = variables['renumber_indices_flag'][0]
    
            if(mvars.renumber_indices_flag):
                mvars.first_index = variables['first_index'][0]

            mvars.renumber_resids_flag = variables['renumber_resids_flag'][0]
            
            if(mvars.renumber_resids_flag):
                mvars.first_resid = variables['first_resid'][0]
       
        mvars.pdb_constraints_flag = variables['pdb_constraints_flag'][0]
       
        if(mvars.pdb_constraints_flag):
            mvars.number_of_constraint_files = variables['number_of_constraint_files'][0]
            mvars.constraint_options = variables['constraint_options'][0]
            mvars.constraint_filenames = variables['constraint_filenames'][0]
            mvars.constraint_fields = variables['constraint_fields'][0]
            mvars.constraint_resets = variables['constraint_resets'][0]
             
        mvars.translation_rotation_flag = variables['translation_rotation_flag'][0]
        
        if(mvars.translation_rotation_flag):
            mvars.translation_rotation_output_filename = variables['translation_rotation_output_filename'][0]
            mvars.pre_center_flag = variables['pre_center_flag'][0] 
            mvars.translation_array = variables['translation_array'][0] 
            mvars.rotation_axes = variables['rotation_axes'][0] 
            mvars.rotation_order = variables['rotation_order'][0] 
            mvars.rotation_array = variables['rotation_array'][0] 
       
        
        mvars.align_pmi_on_axis_flag = variables['align_pmi_on_axis_flag'][0]
        if(mvars.align_pmi_on_axis_flag):

            mvars.pmi_eigenvector = variables['pmi_eigenvector'][0]
            mvars.alignment_vector_axis = variables['alignment_vector_axis'][0]
            mvars.settle_on_plane = variables['settle_on_plane'][0]
            mvars.plane = variables['plane'][0]

        mvars.fasta_utilities_flag = variables['fasta_utilities_flag'][0]
        

        if(mvars.fasta_utilities_flag):
            mvars.fasta_input_option = variables['fasta_input_option'][0]

            if(mvars.fasta_input_option == 'sequence'):
                mvars.fasta_input_sequence = variables['fasta_input_sequence'][0] 
            elif(mvars.fasta_input_option == 'file'):
                mvars.fasta_input_file = variables['fasta_input_file'][0] 

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
            mvars.molecule.read_pdb(mvars.pdbfile)

            if(mvars.pdb_constraints_flag):

                mvars.constraint_options = [x.strip() for x in mvars.constraint_options.split(',')]
                mvars.constraint_fields = [x.strip() for x in mvars.constraint_fields.split(',')]
                mvars.constraint_filenames = [x.strip() for x in mvars.constraint_filenames.split(',')]
                mvars.constraint_resets = [x.strip() for x in mvars.constraint_resets.split(',')]
                
                dum = [] 
                for value in mvars.constraint_resets:
                    dum.append(eval(value))
                mvars.constraint_resets = dum

            if(mvars.translation_rotation_flag):

                mvars.translation_array = numpy.array(mvars.translation_array, numpy.float)


        if(mvars.fasta_utilities_flag):
            pass


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

        if(mvars.renumber_flag):
            if(mvars.renumber_indices_flag):
                if(mvars.renumber_resids_flag):
                    mvars.molecule.renumber(index=mvars.first_index, resid=mvars.first_resid)
                else:
                    mvars.molecule.renumber(index=mvars.first_index)
            elif(mvars.renumber_resids_flag):
                mvars.molecule.renumber(resid=mvars.first_resid)

            mvars.molecule.write_pdb(os.path.join(self.runpath, mvars.renumber_output_filename), frame, 'w')

        elif(mvars.pdb_constraints_flag):
            for i in xrange(mvars.number_of_constraint_files):
                mvars.molecule.make_constraint_pdb(os.path.join(self.runpath, mvars.constraint_filenames[i]), mvars.constraint_options[i], field=mvars.constraint_fields[i], reset=mvars.constraint_resets[i])

        elif(mvars.translation_rotation_flag):

            if(mvars.pre_center_flag):
                mvars.molecule.center(frame)

            mvars.molecule.translate(frame, mvars.translation_array)            

            ### ROTATION CODE GOES HERE


    
            mvars.molecule.write_pdb(os.path.join(self.runpath, mvars.translation_rotation_output_filename), frame, 'w')

        elif(mvars.align_pmi_on_axis_flag):

            mvars.molecule.center(frame)
            uk, ak, I = mvars.molecule.calcpmi(frame) 
            ''' right handed coordinate frame '''
            ak[0] = -1.0 * ak[0]
            #print('ak = ', ak)
            ''' check if right handed coordinate frame '''
            if numpy.linalg.det([ak]) < 0:
                print('determinant not equal +1')
                print(numpy.linalg.det(ak))
                sys.exit()
            ak = ak[mvars.pmi_eigenvector - 1] 
            #print('ak = ', ak)
            if(mvars.alignment_vector_axis == 'x'):
                axis = numpy.array([1.0, 0.0, 0.0])
            elif(mvars.alignment_vector_axis == 'y'):
                axis = numpy.array([0.0, 1.0, 0.0])
            elif(mvars.alignment_vector_axis == 'z'):
                axis = numpy.array([0.0, 0.0, 1.0])

            #print('axis = ', axis)
            rotvec = numpy.cross(ak, axis)
            sine = numpy.linalg.norm(rotvec)
            rotvec = rotvec/numpy.linalg.norm(rotvec)
            cosine = numpy.dot(ak, axis)
            try:
                theta = math.atan(sine/cosine)
            except:
                print('cosine = 0\nstopping here\n\n')
                sys.exit()
            #print('theta = ', theta)
            #print('rotvec = ', rotvec)
            r1 = rotvec[0] ; r2 = rotvec[1] ; r3 = rotvec[2]
            mvars.molecule.general_axis_rotate(frame, theta, r1,r2,r3) 
            mvars.molecule.write_pdb('junk.pdb', frame, 'w')
            
            #mvars.plane 
            #mvars.settle_on_plane 
        
        return
    

    def epilogue(self):
        """
        method to print out and move results to appropriate places
        """

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        pgui('New PDB file(s) saved in %s directory\n\n' % (self.runpath+os.path.sep)) 
        self.run_utils.clean_up(log) 
        pgui('\nrun json inputs saved to:\n    %s\n' % os.path.join(self.runpath, self.parmfile)) 
        pgui('\nrun log output saved to:\n    %s\n' % os.path.join(self.runpath, self.logfile)) 
        pgui("\n\n") 
        pgui("%s \n" % ('=' * 60)) 
        
        #time.sleep(2) 

        return
