'''
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
'''

from __future__ import division

import sys, os, random, logging, numpy, string, shutil,time

import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities

sys.path.append('./')
import altens_util as altens_util

'''
    ALTENS is the module that calculates alignment tensor extracted from RDCs

    From Professor David Fushman, UMD 

    J. Magn. Reson. 168, 336-345 (2004)
    J. Magn. Reson. 201, 25-33 (2009)
    Prog. Nuc. Magn. Res. Spec. 44, 189-214 (2004)

'''

if sasconfig.__level__ == "DEBUG": DEBUG = True

app = 'altens'

class module_variables():
    def __init__(self, parent = None):
        self.app = app

class altens_input_variables():

    def __init__(self, parent = None):
        pass

class altens():

    def __init__(self, parent = None):
        pass

    def main(self, input_variables, txtOutput):

        '''
        main method to manage simulation
        '''

        self.mvars = module_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.run()

        self.epilogue()

        return

    def unpack_variables(self, variables):
        '''
        method to extract variables into system wise class instance
        '''

        mvars = self.mvars
        self.log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.pdbfile = variables['pdbfile'][0]

        '''
       Current version reads a pdb file of single frame.
       We will use dcd file to treat Ensemble of configurations if David Fushiman is ready.
       mvars.dcdfile = variables['dcdfile'][0]
        '''

        mvars.rdc_input_file = variables['rdc_input_file'][0]
#        mvars.nh_vector_coordinate_file = variables['nh_vector_coordinate_file'][0]

        mvars.residue_list_file = variables['residue_list_file'][0]

        '''
       elif: residue_list will be generated from residue list from input pdb except the N terminal: 
       Note Residue numbering in list file starts from 1 
        '''
        mvars.use_monte_carlo_flag = variables['use_monte_carlo_flag'][0]

        ''' 
        method to MC run for error analysis
        '''

        if mvars.use_monte_carlo_flag:
            mvars.number_of_monte_carlo_steps = variables['number_of_monte_carlo_steps'][0]
            mvars.seed = variables['seed'][0]

#        print "mvars.seed = ", mvars.seed

        self.log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to initialize input variables 
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui
        mvars = self.mvars

        return

    def run(self):
        '''
        method to perform Altens calculation
        '''

        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui

        log.debug('in Altens')

        pgui("\n"+"="*60+" \n")
        pgui("DATA FROM RUN: %s \n\n" %(time.ctime()))
        #pgui('>>> starting Altens\n')

        altens_util.altens_core(self,app)

        return

    def epilogue(self):
        '''
        method to print out computational results and to move output data
        to appropriate places.
        '''

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        self.run_utils.clean_up(log)

        pgui("\n"+"="*60+" \n")
        time.sleep(0.1)

        return

