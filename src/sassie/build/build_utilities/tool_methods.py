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
import string
import time

import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import sassie.util.basis_to_python as basis_to_python

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
        mvars.pdbfile = variables['pdbfile'][0]
        
        #mvars.seed = variables['seed'][0]
        mvars.seed = [0]

        return

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

        return

    def initialization(self):
        """
        method to prepare for data for main method
        """

        log = self.log
        log.debug('in initialization')
        mvars = self.mvars


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
        
        time.sleep(2) 

        return
