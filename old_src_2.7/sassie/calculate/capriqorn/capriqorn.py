'''
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
'''

from __future__ import division

import sys, os, random, logging, numpy, string, shutil,time

import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import capriqorn_utils as capriqorn_utils
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'Capriqorn_library'))

'''
    Capriqorn is the module that contains the functions that are used to calculate the neutron or x-ray scattering profile, as well as additional useful output such as p(r) and Vc, based on the user given structure.

'''

if sasconfig.__level__ == "DEBUG": DEBUG = True

app = 'capriqorn'

class module_variables():
    def __init__(self, parent = None):
        self.app = app

class capriqorn_input_variables():
    def __init__(self, parent = None):
        pass

class capriqorn():

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
        mvars.dcdfile = variables['dcdfile'][0]
        if mvars.dcdfile[-3:] == 'dcd':
            mvars.xstfile = variables['xstfile'][0]
        mvars.aliasfile = variables['aliasfile'][0]
        mvars.create_alias_flag = variables['create_alias_flag'][0]
            
        mvars.number_q_values = variables['number_q_values'][0]
        mvars.q_max = variables['q_max'][0]
        
        self.log.debug(vars(mvars))

        return


    def create_alias_file(self):
        pass


        return

    def initialization(self):
        '''
        method to initialize the input variables for capriqorn API
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui
        mvars = self.mvars

        ''' directory preparation '''
        output_folder = os.path.join(mvars.runname,app)
        if os.path.isdir(output_folder):
            shutil.rmtree(output_folder)
            os.mkdir(output_folder)

        ''' initialize capriqorn_input_variables '''
        scvars = capriqorn_input_variables()

        ''' get the sasmol object and save it to capriqorn_input_variables '''
        scvars.runname = mvars.runname
        scvars.pdbfile = mvars.pdbfile
        scvars.dcdfile = mvars.dcdfile
        if scvars.dcdfile[-3:] == 'dcd':
            scvars.xstfile = mvars.xstfile
        scvars.aliasfile = mvars.aliasfile
        scvars.output_folder = output_folder
        mol = sasmol.SasMol(0)
        mol.read_pdb(mvars.pdbfile)

        ''' check if alias.dat file needs to be created '''
        if mvars.create_alias_flag:
            mvars.name = mol.name()
            self.create_alias_file()
        
        try:
            print mvars.dcdfile
            if(mvars.dcdfile.strip()[-3:] == 'dcd'):
                print ">> input file is a DCD file"
                mol.read_dcd(mvars.dcdfile)
                scvars.number_of_frames = mol.number_of_frames()
                self.intype = 'dcd'
            elif(mvars.dcdfile[-3:] =='pdb'):
                print ">> input file is a PDB file"
                mol.read_pdb(mvars.pdbfile)
                scvars.number_of_frames = mol.number_of_frames()
                self.intype = 'pdb'
            elif(mvars.dcdfile.strip()[-6:] == 'crdbox'):
                print ">> input file is a traj file"
                scvars.number_of_frames = capriqorn_utils.get_crdbox_number_of_frames(mvars.dcdfile)
                self.intype = 'traj'
        except:
            message='input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
            message+=' :  stopping here'
            pgui(message)
            sys.exit(1)
        scvars.mol = mol
        scvars.coor= mol.coor()
        scvars.intype = self.intype

        ''' get the arrays of atomic scattering, Q-values, distance values, and save them to capriqorn_input_variables '''
        capriqorn_utils.prepare_capriqorn_inputs(mol, mvars, scvars, output_folder)
        
        self.scvars = scvars

        return

    def run(self):
        '''
        method to perform Capriqorn calculation
        '''

        log = self.log
        mvars = self.mvars
        scvars = self.scvars
        pgui = self.run_utils.print_gui
        output_folder = scvars.output_folder

        frame = 0

        log.debug('in Capriqorn')

        pgui("\n"+"="*60+" \n")
        pgui("DATA FROM RUN: %s \n\n" %(time.ctime()))
        #pgui('>>> starting Capriqorn\n')

        """
        for frame in xrange(scvars.number_of_frames):
            ''' compute and get the results '''
            capriqorn_utils.calculate(frame)

            fraction_done = (frame+1)/float(scvars.number_of_frames)
            time.sleep(0.01)
            pgui('STATUS\t'+str(fraction_done))
        """
        # display progress
        report_string = 'STATUS\t%f' % 0.1
        pgui(report_string)

        capriqorn_utils.calculate(self, mvars, scvars, output_folder)
        
        #pgui('\nProcessed %d DCD frame(s)\n'%scvars.number_of_frames)

        pgui('\nData stored in directory: %s\n'%output_folder)

        return

    def epilogue(self):
        '''
        method to print out simulation results and to move results
        to appropriate places.
        '''

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        self.run_utils.clean_up(log)

        #pgui('\n%s IS DONE\n' % app)
        pgui("\n"+"="*60+" \n")
        
        # display progress
        report_string = 'STATUS\t%f' % 1.0
        pgui(report_string)
        
        time.sleep(0.1)

        return

