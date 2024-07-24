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
import time
import platform
import numpy
import io
import json

import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig

import sassie.analyze.asaxs.asaxs_methods.read_files as read_files

'''
        ASAXS is the module that carries out a calculation of
        anamolous small-angle X-ray scattering from user-supplied
        structure files.

        This module is called from ASAXS from the main
        GUI through the gui_mimic_asaxs.py script.


        REFERENCE:

        Pinfield, V. J., & Scott, D. J. (2014). 
        Anomalous Small Angle X-Ray Scattering Simulations: Proof of Concept for Distance Measurements for 
        Nanoparticle-Labelled Biomacromolecules in Solution. 
        PLOS ONE, 9(4), e95664. 
        https://doi.org/10.1371/journal.pone.0095664

'''

DEBUG = True

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'asaxs'


class module_variables():

    def __init__(self, parent=None):
        self.app = app


class asaxs_variables():

    def __init__(self, parent=None):
        pass


class asaxs():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.module_variables = module_variables()

        self.asaxs_variables = asaxs_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.asaxs()

        self.epilogue()

        return

    def unpack_variables(self, variables):
        '''
        method to extract variables into system wide class instance
        '''

        log = self.log
        mvars = self.module_variables
        log.debug('in unpack_variables')

        mvars.run_name = variables['run_name'][0]
        mvars.pdb_file_name = variables['pdb_file_name'][0]

        mvars.dcd_file_flag = variables['dcd_file_flag'][0]

        mvars.dcd_file_name = None

        if mvars.dcd_file_flag:
            mvars.dcd_file_name = variables['dcd_file_name'][0]

        mvars.gaussian_noise_fraction = variables['gaussian_noise_fraction'][0]

        mvars.number_of_q_values = variables['number_of_q_values'][0]
        mvars.maximum_q_value = variables['maximum_q_value'][0]

        mvars.number_of_d_values = variables['number_of_d_values'][0]
        mvars.maximum_d_value = variables['maximum_d_value'][0]

        mvars.nanocluster_flag = variables['nanocluster_flag'][0]

        if mvars.nanocluster_flag:
        # choice 1
            mvars.nanocluster_atom_name = variables['nanocluster_atom_name'][0]
            mvars.nanocluster_radius = variables['nanocluster_radius'][0]

        mvars.central_energy_value = variables['central_energy_value'][0]
        mvars.energy_range = variables['energy_range'][0]
        mvars.number_of_energy_values = variables['number_of_energy_values'][0]

        mvars.number_of_d_values_for_pair_distribution = svariables['number_of_d_values_for_pair_distribution'][0]
        mvars.maximum_d_value_for_pair_distribution = svariables['maximum_d_value_for_pair_distribution'][0]

        mvars.alternative_scattering_calculator_flag = svariables['alternative_scattering_calculator_flag'][0]

        mvars.crysol_file_name = None
        mvars.sascalc_file_name = None
        mvars.experimental_data_file_name = None

        if mvars.alternative_scattering_calculator_flag:
        # choice 1
            mvars.crysol_flag = svariables['crysol_file_flag'][0]
            mvars.sascalc_file_flag = svariables['sascalc_file_flag'][0]
            mvars.experimental_data_flag = svariables['experimental_data_flag'][0]


        if mvars.crysol_file_flag:
            mvars.crysol_file_name = svariables['crysol_file_name'][0]
        elif mvars.sascalc_file_flag:
            mvars.sascalc_file_name = svariables['sascalc_file_name'][0]
        elif mvars.sascalc_file_flag:
            mvars.experimental_data_file_name = svariables['experimenal_data_file_name'][0] 

        mvars.plot_flag = variables['plot_flag'][0]

        log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to prepare for asaxs 
        '''

        log = self.log
        log.debug('in initialization')
        
        pgui = self.run_utils.print_gui
        mvars = self.module_variables
        avars = self.asaxs_variables

        avars.interpath = mvars.run_name + '/asaxs/'
        direxist = os.path.exists(avars.interpath)
        if(direxist == 0):
            os.system('mkdir -p ' + avars.interpath)

        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in range(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))
        
        return

    def save_data_to_plot_as_json(self):

        mvars = self.module_variables
        avars = self.asaxs_variables

        # open the json output file for writing
#        json_outfile = io.open(avars.interpath + mvars.ofile[:-3] + 'json', 'w')
       
#        data_dict = {
#            'q': [],
#            'iq': [],
#            'iq_error': [],
#            'original_q': [],
#            'original_iq': [],
#            'original_iq_error': [],
#            'signal_to_noise_cutoff_value': avars.cutval
#        }

        #for i in range(len(divars.io_tally)):
        ## Append new values to the lists in the dictionary
        #    data_dict['q'].append(divars.io_tally[i][0])
        #    data_dict['iq'].append(divars.io_tally[i][1])
        #    data_dict['iq_error'].append(divars.io_tally[i][2])
#
#        for i in range(len(divars.odata)):
#            data_dict['original_q'].append(divars.odata[i][0])
#            data_dict['original_iq'].append(divars.odata[i][1])
#            data_dict['original_iq_error'].append(divars.odata[i][2])
#
#        json_data = json.dumps(data_dict)

#        json_outfile.write(json_data)
#        json_outfile.close()

        return


    def asaxs(self):
        '''
        ASAXS is the function to read in variables from GUI input and
        carries out analysis and fitting of ASAXS data.

        INPUT:  variable descriptions:

            run_name:		    project name

            pdb_file_name:      name of reference PDB file

            dcd_file_flag:      boolean to indicate use of DCD trajectory file
            
            dcd_file_name:      (optional) name of DCD trajectory file

            gaussian_noise_fraction: amount of noise to represent experimental uncertainty

            number_of_q_values: number of q-values used in scattering intensity calculations

            maximum_q_value: maximum q-value used in scattering intensity calculations (1/Angstrom)

            number_of_d_values: number of D-values used in scattering intensity calculations 

            maximum_D_value: maximum D-value used in scattering intensity calculations (Angstrom)

            nanocluster_flag:   boolean to indicate modeling single atom in PDB file as a cluster

            nanocluster_atom_name: atom name to define scattering and cluster positions

            nanocluster_radius: user supplied radius (Angstrom)

            central_energy_value: energy near absorption edge (eV)

            energy_range: range of energies to model (+/- energy_range / 2)

            number_of_energy_values: number of energy values to use

            number_of_d_values_for_pair_distribution: number of D-values to consider in final inversion 

            maximum_d_value_for_pair_distribution: maximum D-value to use in final inversion (Angstrom)

            alternative_scattering_calculator_flag: boolean to indicate using calcuated or measured scattering data 

            crysol_file_flag: boolean to use a Crysol file 
            crysol_file_name: (optional) name of Crysol file

            sascalc_file_flag: boolean to use a SasCalc file
            sascalc_file_name: (optional) name of SasCalc file

            experimental_data_file_flag: boolean to use an experimental data file 
            experimental_data_file_name: (optional) name of experimental data file

            plot_flag: (optional) create plot data to save to file 


        OUTPUT:

                file is stored in "run_name"/asaxs directory

                ofile:                  output filename

        '''

        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in asaxs')

        mvars = self.module_variables
        avars = self.asaxs_variables

        read_files.read_structure_files(self)        

        if(mvars.crysol_file_flag):
            read_files.read_crysol_files(self)        

        #outfile2 = open(divars.interpath + mvars.ofile, 'w')

        #self.save_data_to_plot_as_json()

        time.sleep(0.1)

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        time.sleep(1.0)

        log = self.log
        mvars = self.module_variables
        avars = self.asaxs_variables
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        #pgui("\nAsaxs data were written to %s\n" %
        #     ('./' + avars.interpath + mvars.ofile))

        pgui('\nASAXS IS DONE\n')

        st = ''.join(['=' for x in range(60)])
        pgui("\n%s \n" % (st))

        time.sleep(1.0)

        ''' display progress '''

        fraction_done = 1
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        time.sleep(1.0)

        self.run_utils.clean_up(log)

        return
