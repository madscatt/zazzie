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
import sassie.util.basis_to_python as basis_to_python
#import sassie.calculate.sascalc_utils as sascalc_utils
import sascalc_utils as sascalc_utils
import sassie.calculate.sascalc.sascalc_library.sascalc_lib as sascalc_lib

'''
    SASCALC is the module that contains the functions that are used to calculate the neutron or x-ray scattering profile, as well as additional useful output such as p(r) and Vc, based on the user given structure.

'''

if sasconfig.__level__ == "DEBUG": DEBUG = True

app = 'sascalc'

class module_variables():
    def __init__(self, parent = None):
        self.app = app

class sascalc_input_variables():
    def __init__(self, parent = None):
        pass

class sascalc():

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
        mvars.xon = variables['xon'][0]
        if mvars.xon in ['neutron','neutron_and_xray']:
            mvars.number_contrast_points = variables['number_contrast_points'][0]
            mvars.D2O_percentage_array = variables['D2O_percentage_array'][0]
            mvars.I0_array = variables['I0_array'][0]

            mvars.number_exH_regions = variables['number_exH_regions'][0]
            mvars.exH_basis_string_array = variables['exH_basis_string_array'][0]
            mvars.fraction_exH_array = variables['fraction_exH_array'][0]

            mvars.number_deuteration_regions = variables['number_deuteration_regions'][0]
            mvars.deuterated_basis_string_array = variables['deuterated_basis_string_array'][0]
            mvars.fraction_deuterated_array = variables['fraction_deuterated_array'][0]

        if mvars.xon in ['xray','neutron_and_xray']:
            mvars.xray_number_contrast_points = variables['xray_number_contrast_points'][0]
            mvars.xray_D2O_percentage_array = variables['xray_D2O_percentage_array'][0]
            mvars.xray_I0_array = variables['xray_I0_array'][0]
            
        mvars.number_q_values = variables['number_q_values'][0]
        mvars.q_max = variables['q_max'][0]
        mvars.number_r_values = variables['number_r_values'][0]
        mvars.golden_vector_method_option = variables['golden_vector_method_option'][0]
        if mvars.golden_vector_method_option == 'fixed':
            mvars.number_golden_vectors = variables['number_golden_vectors'][0]
        elif mvars.golden_vector_method_option == 'converge':
            mvars.golden_vector_method_converge_tolerance = variables['golden_vector_method_converge_tolerance'][0]

        mvars.solvent_volume = variables['solvent_volume'][0]
        mvars.VDW_scaling_factor = variables['VDW_scaling_factor'][0]
        mvars.complex_amplitudes = variables['complex_amplitudes'][0]

        self.log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to initialize the input variables for sascalc API
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui
        mvars = self.mvars

        ## code consistent with TAMC
        if mvars.xon in ['neutron','neutron_and_xray']:
            exH_basis_string_array = mvars.exH_basis_string_array
            exH_basis_string_array_python = []
            for basis in exH_basis_string_array:
                exH_basis_string_array_python.append(basis_to_python.parse_basis(basis))
            mvars.exH_basis_string_array_python = exH_basis_string_array_python

            deuterated_basis_string_array = mvars.deuterated_basis_string_array
            deuterated_basis_string_array_python = []
            for basis in deuterated_basis_string_array:
                deuterated_basis_string_array_python.append(basis_to_python.parse_basis(basis))
            mvars.deuterated_basis_string_array_python = deuterated_basis_string_array_python

        ''' directory preparation '''
        output_folder = os.path.join(mvars.runname,app)
        if os.path.isdir(output_folder):
            shutil.rmtree(output_folder)
            os.mkdir(output_folder)

        ''' initialize sascalc_input_variables '''
        scvars = sascalc_input_variables()

        ''' get the sasmol object and save it to sascalc_input_variables '''
        scvars.runname = mvars.runname
        scvars.pdbfile = mvars.pdbfile
        scvars.dcdfile = mvars.dcdfile
        scvars.output_folder = output_folder
        mol = sasmol.SasMol(0)
        mol.read_pdb(mvars.pdbfile)
        
        try:
            if(mvars.dcdfile[-3:] == 'dcd'):
                print ">> input file is a DCD file"
                mol.read_dcd(str(mvars.dcdfile))
                self.intype = 'dcd'
            elif(mvars.dcdfile[-3:] =='pdb'):
                print ">> input file is a PDB file"
                mol.read_pdb(mvars.dcdfile)
                nf = mol.number_of_frames()
                self.intype = 'pdb'
        except:
            message='input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
            message+=' :  stopping here'
            pgui(message)
            sys.exit(1)
        scvars.mol = mol
        scvars.coor= mol.coor()
        scvars.number_of_frames = mol.number_of_frames()
        scvars.number_of_atoms = mol.natoms()
        if 'frames_per_batch' in vars(mvars):
            scvars.frames_per_batch = mvars.frames_per_batch
        else:
            max_gpu_memory = sasconfig.__gpu_memory__*0.2*1024*1024 ## NOTE currently allows only 20% of GPU memory to be used
            memory_per_frame = mol.number_of_frames() * mol.natoms() * 3 * 8 ## NOTE current data precision in coordinate array is float64
            scvars.frames_per_batch = max(int(max_gpu_memory/memory_per_frame), 1)
            #scvars.frames_per_batch = int(max(mol.number_of_frames()/100, 1)) ##@NOTE to ZHL default to 100 batches
            #scvars.frames_per_batch = max(mol.number_of_frames()/1, 1) ##@NOTE to ZHL default to 100 batches

        ''' get the arrays of atomic scattering, Q-values, distance values, and save them to sascalc_input_variables '''
        sascalc_utils.prepare_sascalc_inputs(mol, mvars, scvars, output_folder)

        ''' save golden_vector variables to sascalc_input_variables '''
        scvars.golden_vector_method_option =  mvars.golden_vector_method_option
        if mvars.golden_vector_method_option == 'fixed':
            scvars.number_golden_vectors = mvars.number_golden_vectors;
        if mvars.golden_vector_method_option == 'converge':
            scvars.golden_vector_method_converge_tolerance = mvars.golden_vector_method_converge_tolerance
        ''' save xray/neutron flag information to sascalc_input_variables '''
        scvars.xon = mvars.xon

        self.scvars = scvars

        return

    def run(self):
        '''
        method to perform sascalc calculation
        '''

        log = self.log
        mvars = self.mvars
        scvars = self.scvars
        mol = scvars.mol
        pgui = self.run_utils.print_gui

        frame = 0
        frames_per_batch = scvars.frames_per_batch

        log.debug('in sascalc')

        pgui("\n"+"="*60+" \n")
        pgui("DATA FROM RUN: %s \n\n" %(time.asctime( time.gmtime( time.time() ) ) ))
        pgui('>>> starting sascalc\n')


        print "initializing sascalc object..."
        ''' create the sascalc object '''

        sascalc_object = sascalc_lib.SasCalc(mol, mvars, scvars)

        sascalc_coor = numpy.transpose(mol.coor(),axes=(0,2,1))

        print "computing first batch..."
        ''' compute and get the results '''
        results = sascalc_object.calculate(sascalc_coor, 2)
        ''' save output '''
        results.save('first batch')

        print "computing second batch..."
        ''' compute and get the results '''
        results = sascalc_object.calculate(sascalc_coor[2:], 1)
        ''' save output '''
        results.save('second batch')


        """
        ''' compute and get the results '''
        print "running sascalc "
        results = sascalc_object.calculate(sascalc_coor, 1)
        ''' save output '''
        print "saving result"
        results.save('batch')
        """


        """
        for frame in xrange(scvars.number_of_frames):

            ''' compute and get the results '''
            results = sascalc_object.calculate(mol.coor()[frame], 1)

            ''' save output and report progress '''
            results.save('frame_%d'%frame)

            ''' update fraction done '''
            fraction_done = (frame+1)/float(scvars.number_of_frames)
            time.sleep(0.01)
            pgui('STATUS\t'+str(fraction_done))
        """

        
        pgui('\nProcessed %d DCD frame(s)\n'%scvars.number_of_frames)

        ''' clean up '''
        #sascalc_object.clean()

        return

    def epilogue(self):
        '''
        method to print out simulation results and to move results
        to appropriate places.
        '''
        scvars = self.scvars

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        print 'generating the final h5 file...'
        #sascalc_utils.reoganize_h5(scvars.output_folder)

        self.run_utils.clean_up(log)

        #pgui('\n%s IS DONE\n' % app)
        pgui("\n"+"="*60+" \n")
        time.sleep(0.1)

        return

