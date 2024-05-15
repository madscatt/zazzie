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
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig

'''
        INTERPOLATE is the module that calculates an approximate data set
        with a defined number of points and grid spacing.  The output of
        this file is used to compare to synthetic profiles generated in
        subsequent modules.

        This module is called from Data Interpolation from the main
        GUI through the graphical_interpolate.py script.


        REFERENCE:

        Numerical Recipes: The Art of Scientific Computing
        Third Edition (2007), 1256 pp.
        Cambridge University Press
        ISBN-10: 0521880688
'''


if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'data_interpolation'


class module_variables():

    def __init__(self, parent=None):
        self.app = app


class data_interpolation_variables():

    def __init__(self, parent=None):
        pass


class data_interpolation():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.module_variables = module_variables()

        self.data_interpolation_variables = data_interpolation_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.interpolate()

        self.epilogue()

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

    def readfile(self, data_file):
        '''
        READFILE is the function to read NCNR SANS data files
        '''
        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in readfile')

        mvars = self.module_variables
        divars = self.data_interpolation_variables

        fake_error = False
        error_magnitude = 0.1

        divars.x = []
        divars.y = []
        divars.z = []
        divars.nval = 1
        divars.x.append(0.0)
        divars.y.append(mvars.io)
        divars.z.append(mvars.ioe)

        for line in data_file:
            #this_line = string.split(line)
            this_line = line.split()
            try:
                qval = locale.atof(this_line[0])
                ival = locale.atof(this_line[1])
                # ival=abs(locale.atof(this_line[1]))

                try:
                    error_value = locale.atof(this_line[2])
                except:
                    error_value = error_magnitude * ival
                    fake_error = True
                # if(error_value == 0.0):
                #    error_value = error_value + 1E-5

                divars.x.append(qval)
                divars.y.append(ival)
                divars.z.append(error_value)
                divars.nval = divars.nval + 1

            except:
                pass

        if fake_error:
            message = 'Error values in I(q) were set to be ' + str(100 * error_magnitude) + '% of I(0) values since appropriate values were not found'
            pgui(message)

        return

    def spline(self, array):
        '''
        SPLINE is the function to calculate an approximate data set
        '''

        log = self.log
        log.debug('in spline')

        mvars = self.module_variables
        divars = self.data_interpolation_variables

        u = numpy.zeros(divars.nval)
        divars.array2 = numpy.zeros(divars.nval)
        maxval = 0.99E30
        nmax = 500
        if(divars.yp1 > maxval):
            divars.array2[0] = 0.0
            u[0] = 0.0
        else:
            divars.array2[0] = -0.5
            u[0] = (3.0 / (divars.x[1] - divars.x[0])) * ((array[
                    1] - array[0]) / (divars.x[1] - divars.x[0]) - divars.yp1)
        for i in range(1, divars.nval - 2):
            sig = (divars.x[i] - divars.x[i - 1]) / (divars.x[
                   i + 1] - divars.x[i - 1])
            p = sig * divars.array2[i - 1] + 2.0
            divars.array2[i] = (sig - 1.0) / p
            u[i] = (array[i + 1] - array[i]) / (divars.x[i + 1] - divars.x[i]) - (array[i] - array[i - 1]) / (divars.x[i] - divars.x[i - 1])
            u[i] = (6.0 * u[i] / (divars.x[i + 1] - divars.x[
                    i - 1]) - sig * u[i - 1]) / p

        if(divars.ypn > maxval):
            qn = 0.0
            un = 0.0
        else:
            qn = 0.5
            un = (3.0 / (divars.x[divars.nval - 1] - divars.x[divars.nval - 2])) * (divars.ypn -
                                                                                    (array[divars.nval - 1] - array[divars.nval - 2]) / (divars.x[divars.nval - 1] - divars.x[divars.nval - 2]))

        divars.array2[divars.nval - 1] = (un - qn * u[divars.nval - 2]) / (
            qn * divars.array2[divars.nval - 2] + 1.0)

        for k in range(divars.nval - 2, -1, -1):
            divars.array2[k] = divars.array2[k] * divars.array2[k + 1] + u[k]

        return

    def splint(self, array, array2):

        log = self.log
        log.debug('in splint')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        divars = self.data_interpolation_variables

        klo = 0
        khi = divars.nval - 1
        while(khi - klo > 1):
            if((khi - klo) > 1.0):
                k = int((khi + klo) / 2)
                if(divars.x[k] > divars.ux):
                    khi = k
                else:
                    klo = k
        h = divars.x[khi] - divars.x[klo]
        if(h == 0.0):
            pgui('ERROR: BAD INPUT TO ROUTINE SPLINT')
        a = (divars.x[khi] - divars.ux) / h
        b = (divars.ux - divars.x[klo]) / h
        divars.nyval = a * array[klo] + b * array[khi] + \
            ((a * a * a - a) * array2[klo] +
             (b * b * b - b) * array2[khi]) * (h * h) / 6.0

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

#        ''' display progress '''
#
#        fraction_done = 1
#        report_string = 'STATUS\t' + str(fraction_done)
#        pgui(report_string)


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
