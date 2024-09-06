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

### TODO: removed non-used imports

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

import eros_util as eros_util

#
#       EROS
#
#       10/22/2018      --      initial coding                  :       jc
#
# LC     1         2         3         4         5         6         7
# LC567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                     *      **

"""
    EROS is a module that contains the functions that are
    used to obtain a maximum entropy solution for an ensemeble
    of scattering profiles to a experimental scattering profile.

    Based on Rozycki, Kim & Hummer, SAXS Ensemble Refinement of ESCRT-III
    CHMP3 Conformational Transitions, Structure 19, 109-116 (2011)

"""

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'eros'


class module_variables():

    def __init__(self, parent=None):
        self.app = app

class eros_variables():

    def __init__(self, parent=None):
        self.app = app

class eros():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):
        """
        main method to manage module
        """

        self.mvars = module_variables()
        self.evars = eros_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.run_eros()

        self.epilogue()

        return

    def unpack_variables(self, variables):
        """
        method to extract variables into system wide class instance
        """

        mvars = self.mvars
        self.log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]

        mvars.iq_data_path = variables['iq_data_path'][0]

        mvars.goal_iq_data_file = variables['goal_iq_data_file'][0]

        mvars.io = variables['io'][0]
       
        mvars.number_of_files_to_use = variables['number_of_files_to_use'][0]
        mvars.number_of_monte_carlo_steps = variables['number_of_monte_carlo_steps'][0]
        mvars.weight_step_size_fraction = variables['weight_step_size_fraction'][0]
        mvars.theta = variables['theta'][0]
        mvars.beta = variables['beta'][0]
        mvars.reduced_x2 = variables['reduced_x2'][0]
        
        mvars.local_bokeh_server = variables['local_bokeh_server'][0]

        mvars.seed = variables['seed'][0]

        return

    def initialization(self):
        """
        method to prepare for data for main method
        """

        log = self.log
        log.debug('in initialization')
        mvars = self.mvars

        return 

    def run_eros(self):
        """
        main method of module
        """

        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui

        # start gui output
        pgui("\n%s \n" % ('=' * 60))
        pgui("DATA FROM RUN: %s \n\n" % time.asctime(time.gmtime(time.time())))

        frame = 0

        log.debug('in eros')

        """ set up random seed """

        if mvars.seed[0] == 1:
            from numpy.random import RandomState
            mvars.seed_object = RandomState(mvars.seed[1])
        else:
            mvars.seed_object = -1

        """ main loop """

        eros_util.main(self)


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

        self.run_utils.clean_up(log)
        pgui('\nrun json inputs saved to:\n    %s\n' %
             os.path.join(self.runpath, self.parmfile))
        pgui('\nrun log output saved to:\n    %s\n' %
             os.path.join(self.runpath, self.logfile))
        pgui("\n\n")
        pgui("%s \n" % ('=' * 60))

        #TODO: uncomment the sleep line when in production
        #time.sleep(2)

        return
