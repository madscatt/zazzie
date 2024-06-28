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


if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'asaxs'


class module_variables():

    def __init__(self, parent=None):
        self.app = app


class data_interpolation_variables():

    def __init__(self, parent=None):
        pass


class asaxs():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.module_variables = module_variables()

        #self.data_interpolation_variables = data_interpolation_variables()

        #self.run_utils = module_utilities.run_utils(app, txtOutput)

        #self.run_utils.setup_logging(self)

        #self.log.debug('in main')

        #self.unpack_variables(input_variables)

        #self.run_utils.general_setup(self)

        #self.initialization()

        #self.interpolate()

        #self.epilogue()

        return

    def unpack_variables(self, variables):
        '''
        method to extract variables into system wise class instance
        '''

        log = self.log
        mvars = self.module_variables
        log.debug('in unpack_variables')

        mvars.run_name = variables['run_name'][0]
        mvars.expdata = variables['expdata'][0]
        mvars.ofile = variables['ofile'][0]
        mvars.io = variables['io'][0]
        mvars.ioe = variables['ioe'][0]
        mvars.dq = variables['dq'][0]
        mvars.maxpoints = variables['maxpoints'][0]
        mvars.plotflag = variables['plotflag'][0]

        log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to prepare for data interpolation
        '''

#mvars:    run_name, expdata, ofile, io, ioe, dq, maxpoints, plotflag

        log = self.log
        log.debug('in initialization')
        readfile = self.readfile
        pgui = self.run_utils.print_gui
        mvars = self.module_variables
        divars = self.data_interpolation_variables

        divars.interpath = mvars.run_name + '/data_interpolation/'
        direxist = os.path.exists(divars.interpath)
        if(direxist == 0):
            os.system('mkdir -p ' + divars.interpath)

        # ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in range(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))
        pgui("Input file used : %s\n" % (mvars.expdata))

        data_file = open(mvars.expdata, 'r').readlines()
#       data_file=open(mvars.expdata,'r').read().splitlines()

    #   following line replaces control M (\r) with \n

        data_file = [x.replace("\r", "\n") for x in data_file]

        readfile(data_file)

        divars.odata = []
        flag = 0
        divars.cut = []
        for i in range(divars.nval):
            divars.odata.append([divars.x[i], divars.y[i], divars.z[i]])
            if(divars.y[i] / divars.z[i] < 2.0 and flag == 0):
                divars.cut.append([divars.x[i - 1], divars.y[i - 1], mvars.io])
                divars.cutval = divars.x[i - 1]
                flag = 1
            elif((i == divars.nval - 1) and flag == 0):
                divars.cut.append([divars.x[i - 1], divars.y[i - 1], mvars.io])
                divars.cutval = divars.x[i - 1]
                flag = 1

        return

    def save_data_to_plot_as_json(self):

        mvars = self.module_variables
        divars = self.data_interpolation_variables

        # open the json output file for writing
        json_outfile = io.open(divars.interpath + mvars.ofile[:-3] + 'json', 'w')
       
        data_dict = {
            'q': [],
            'iq': [],
            'iq_error': [],
            'original_q': [],
            'original_iq': [],
            'original_iq_error': [],
            'signal_to_noise_cutoff_value': divars.cutval
        }

        for i in range(len(divars.io_tally)):
        # Append new values to the lists in the dictionary
            data_dict['q'].append(divars.io_tally[i][0])
            data_dict['iq'].append(divars.io_tally[i][1])
            data_dict['iq_error'].append(divars.io_tally[i][2])

        for i in range(len(divars.odata)):
            data_dict['original_q'].append(divars.odata[i][0])
            data_dict['original_iq'].append(divars.odata[i][1])
            data_dict['original_iq_error'].append(divars.odata[i][2])

        json_data = json.dumps(data_dict)

        json_outfile.write(json_data)
        json_outfile.close()

        return


    def interpolate(self):
        '''
        INTERPOLATE is the function to read in variables from GUI input and
        calculate an approximate data set to be used in subsequent modeling
        steps.

        INPUT:  variable descriptions:

                run_name:		    project name
                expdata:        input NCNR data file (*.sub)
                io:             I(0)
                ioe:            Error in I(0)
                dq:             new delta q
                maxpoints:      number of new points

        OUTPUT:

                file is stored in "run_name"/data_interpolation directory

                ofile:                  output filename

        '''

        log = self.log
        pgui = self.run_utils.print_gui
        splint = self.splint
        spline = self.spline
        log.debug('in interpolate')

        mvars = self.module_variables
        divars = self.data_interpolation_variables

        divars.yp1 = 1.0
        divars.ypn = 1.0

        log.debug('calculating splines')

        try:
            spline(divars.y)
            divars.y2 = divars.array2
            spline(divars.z)
            divars.z2 = divars.array2
        except:
            message = 'Failed to interpolate data: is data file corrupt?'
            pgui(message)
            message += '\nstopping here\n\n'
            pgui(message)
            return

        log.debug('back from splines')

        outfile2 = open(divars.interpath + mvars.ofile, 'w')
        outfile3 = open(divars.interpath + 'stn_' + mvars.ofile, 'w')

        divars.io_tally = []
        outfile2.write('%f\t%f\t%f\n' % (0.0, mvars.io, mvars.ioe))
        outfile3.write('%f\t%f\t%f\n' % (0.0, mvars.io, mvars.ioe))
        divars.io_tally.append([0.0, mvars.io, mvars.ioe])
        divars.ux = 0.00
        pgui("\nSignal to noise cutoff value: %s\n" % (str(divars.cutval)))

        for i in range(mvars.maxpoints - 1):
            divars.ux = divars.ux + mvars.dq
            splint(divars.y, divars.y2)
            divars.ny = divars.nyval
            splint(divars.z, divars.z2)
            divars.nz = divars.nyval
            divars.io_tally.append([divars.ux, divars.ny, divars.nz])
            outfile2.write('%f\t%f\t%f\n' % (divars.ux, divars.ny, divars.nz))
            if(divars.ux <= divars.cutval):
                outfile3.write(
                    '%f\t%f\t%f\n' % (divars.ux, divars.ny, divars.nz))

        outfile2.close()
        outfile3.close()

        self.save_data_to_plot_as_json()

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
        divars = self.data_interpolation_variables
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        pgui("\nInterpolated data were written to %s\n" %
             ('./' + divars.interpath + mvars.ofile))
        pgui("\nInterpolated data with S/N > 2 were written to %s\n\n" %
             ('./' + divars.interpath + 'stn_' + mvars.ofile))
        pgui("\ndelta q = %f (1/A)\n\nnumber of q-points = %i\n\nq-range: 0 to %f (1/A)\n" %
             (mvars.dq, mvars.maxpoints, (mvars.maxpoints - 1) * mvars.dq))


        pgui('\nDATA INTERPOLATION IS DONE\n')

        st = ''.join(['=' for x in range(60)])
        pgui("\n%s \n" % (st))

        time.sleep(2.0)

        ''' display progress '''

        fraction_done = 1
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        time.sleep(2.0)

        self.run_utils.clean_up(log)

        return
