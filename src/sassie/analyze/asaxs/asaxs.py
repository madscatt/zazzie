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

import sassie.analyze.asaxs.read_files as read_files

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

#        mvars.expdata = variables['expdata'][0]
#        mvars.ofile = variables['ofile'][0]
#        mvars.io = variables['io'][0]
#        mvars.ioe = variables['ioe'][0]
#        mvars.dq = variables['dq'][0]
#        mvars.maxpoints = variables['maxpoints'][0]
#        mvars.plotflag = variables['plotflag'][0]

        log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to prepare for asaxs 
        '''

#mvars:    run_name, expdata, ofile, io, ioe, dq, maxpoints, plotflag

        log = self.log
        log.debug('in initialization')
        
        #readfile = self.readfile

        pgui = self.run_utils.print_gui
        mvars = self.module_variables
        avars = self.asaxs_variables

        avars.interpath = mvars.run_name + '/asaxs/'
        direxist = os.path.exists(avars.interpath)
        if(direxist == 0):
            os.system('mkdir -p ' + avars.interpath)

        # ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in range(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))
        
       # pgui("Input file used : %s\n" % (mvars.expdata))

        #data_file = open(mvars.expdata, 'r').readlines()
#       #data_file=open(mvars.expdata,'r').read().splitlines()

    #   #following line replaces control M (\r) with \n

        #data_file = [x.replace("\r", "\n") for x in data_file]

        #readfile(data_file)
#
#        divars.odata = []
#        flag = 0
#        divars.cut = []
#        for i in range(divars.nval):
#            divars.odata.append([divars.x[i], divars.y[i], divars.z[i]])
#            if(divars.y[i] / divars.z[i] < 2.0 and flag == 0):
#                divars.cut.append([divars.x[i - 1], divars.y[i - 1], mvars.io])
#                divars.cutval = divars.x[i - 1]
#                flag = 1
#            elif((i == divars.nval - 1) and flag == 0):
#                divars.cut.append([divars.x[i - 1], divars.y[i - 1], mvars.io])
#                divars.cutval = divars.x[i - 1]
#                flag = 1

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

                #expdata:        input NCNR data file (*.sub)
                #io:             I(0)
                #ioe:            Error in I(0)
                #dq:             new delta q
                #maxpoints:      number of new points

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
