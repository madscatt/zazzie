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
import Gnuplot
import Gnuplot.PlotItems
import Gnuplot.funcutils
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

class data_interpolation_input_variables():

    def __init__(self, parent=None):
        pass

class data_interpolation():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = data_interpolation_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.interpolate(input_variables, txtOutput)

        self.epilogue()

        return

    def unpack_variables(self,variables):
        '''
        method to extract variables into system wise class instance
        '''

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')
        
        mvars.runname = variables['runname'][0]
        mvars.expdata = variables['expdata'][0]
        mvars.ofile = variables['ofile'][0]
        mvars.io = variables['io'][0]
        mvars.ioe = variables['ioe'][0]
        mvars.dq = variables['dq'][0]
        mvars.maxpoints = variables['maxpoints'][0]
        mvars.plotflag = variables['plotflag'][0]

        log.debug(vars(mvars))

        return


    def wait(self,str=None, prompt='Plot will clear in 10 seconds ...\n'):
        '''
        WAIT is the function to prompt the user to clear a plot on a screen
        '''

        if str is not None:
            print str
        try:
            if(platform.system() == "Linux"):
                import curses
                stdscr = curses.initscr()
                stdscr.addstr('press a key to continue')
                c = stdscr.getch()
                curses.endwin()
        except:
            time.sleep(2)


#   pgui performs this function
#    def print_failure(message, txtOutput):
#
#        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#        txtOutput.put(message)
#
#        return


    def readfile(self,data_file):
        '''
        READFILE is the function to read NCNR SANS data files
        '''
        log = self.log
        pgui = self.run_utils.print_gui        
        log.debug('in readfile')

        mvars = self.mvars
        avars = self.avars

        fake_error = False
        error_magnitude = 0.1

        avars.x = []
        avars.y = []
        avars.z = []
        avars.nval = 1
        avars.x.append(0.0)
        avars.y.append(mvars.io)
        avars.z.append(mvars.ioe)

        for line in data_file:
            this_line = string.split(line)
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

                avars.x.append(qval)
                avars.y.append(ival)
                avars.z.append(error_value)
                avars.nval = avars.nval + 1

            except:
                pass

        if fake_error:
            message = 'Error values in I(q) were set to be ' + str(100 * error_magnitude) + '% of I(0) values since appropriate values were not found'
            pgui(message)
            

        return

    def spline(self,array):
        '''
        SPLINE is the function to calculate an approximate data set
        '''

        log = self.log
        log.debug('in spline')

        mvars = self.mvars
        avars = self.avars             

        u = numpy.zeros(avars.nval)
        avars.array2 = numpy.zeros(avars.nval)
        maxval = 0.99E30
        nmax = 500
        if(avars.yp1 > maxval):
            avars.array2[0] = 0.0
            u[0] = 0.0
        else:
            avars.array2[0] = -0.5
            u[0] = (3.0 / (avars.x[1] - avars.x[0])) * ((array[1] - array[0]) / (avars.x[1] - avars.x[0]) - avars.yp1)
        for i in range(1, avars.nval - 2):
            sig = (avars.x[i] - avars.x[i - 1]) / (avars.x[i + 1] - avars.x[i - 1])
            p = sig * avars.array2[i - 1] + 2.0
            avars.array2[i] = (sig - 1.0) / p
            u[i] = (array[i + 1] - array[i]) / (avars.x[i + 1] - avars.x[i]) - (array[i] - array[i - 1]) / (avars.x[i] - avars.x[i - 1])
            u[i] = (6.0 * u[i] / (avars.x[i + 1] - avars.x[i - 1]) - sig * u[i - 1]) / p

        if(avars.ypn > maxval):
            qn = 0.0
            un = 0.0
        else:
            qn = 0.5
            un = (3.0 / (avars.x[avars.nval - 1] - avars.x[avars.nval - 2])) * (avars.ypn -
                                                    (array[avars.nval - 1] - array[avars.nval - 2]) / (avars.x[avars.nval - 1] - avars.x[avars.nval - 2]))

        avars.array2[avars.nval - 1] = (un - qn * u[avars.nval - 2]) / (qn * avars.array2[avars.nval - 2] + 1.0)

        for k in range(avars.nval - 2, -1, -1):
            avars.array2[k] = avars.array2[k] * avars.array2[k + 1] + u[k]
    
        return 


    def splint(self,array,array2):

        log = self.log
        log.debug('in splint')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        klo = 0
        khi = avars.nval - 1
        while(khi - klo > 1):
            if((khi - klo) > 1.0):
                k = int((khi + klo) / 2)
                if(avars.x[k] > avars.ux):
                    khi = k
                else:
                    klo = k
        h = avars.x[khi] - avars.x[klo]
        if(h == 0.0):
            pgui('ERROR: BAD INPUT TO ROUTINE SPLINT')
        a = (avars.x[khi] - avars.ux) / h
        b = (avars.ux - avars.x[klo]) / h
        avars.nyval = a * array[klo] + b * array[khi] + \
            ((a * a * a - a) * array2[klo] +
            (b * b * b - b) * array2[khi]) * (h * h) / 6.0

        return 



    def initialization(self):
        '''
        method to prepare for data interpolation
        '''
        
#mvars:    runname, expdata, ofile, io, ioe, dq, maxpoints, plotflag 

        log = self.log
        log.debug('in initialization')
        readfile = self.readfile
        pgui = self.run_utils.print_gui
        mvars = self.mvars
        avars = self.avars

        avars.interpath = mvars.runname + '/data_interpolation/'
        direxist = os.path.exists(avars.interpath)
        if(direxist == 0):
            os.system('mkdir -p ' + avars.interpath)

        # ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))
        pgui("Input file used : %s\n" % (mvars.expdata))

        data_file = open(mvars.expdata, 'r').readlines()
#       data_file=open(mvars.expdata,'r').read().splitlines()

    #   following line replaces control M (\r) with \n

        data_file = [x.replace("\r", "\n") for x in data_file]

        readfile(data_file)

        avars.odata = []
        flag = 0
        avars.cut = []
        for i in range(avars.nval):
            avars.odata.append([avars.x[i], avars.y[i], avars.z[i]])
            if(avars.y[i] / avars.z[i] < 2.0 and flag == 0):
                avars.cut.append([avars.x[i - 1], avars.y[i - 1], mvars.io])
                avars.cutval = avars.x[i - 1]
                flag = 1
            elif((i == avars.nval - 1) and flag == 0):
                avars.cut.append([avars.x[i - 1], avars.y[i - 1], mvars.io])
                avars.cutval = avars.x[i - 1]
                flag = 1

        return


    def interpolate(self,variables, txtOutput):
        '''
        INTERPOLATE is the function to read in variables from GUI input and 
        calculate an approximate data set to be used in subsequent modeling
        steps. 

        INPUT:  variable descriptions:

                runname:		    project name
                expdata:        input NCNR data file (*.sub)
                io:             I(0) 
                ioe:            Error in I(0) 
                dq:             new delta q 
                maxpoints:      number of new points

        OUTPUT:

                file is stored in "runname"/data_interpolation directory

                ofile:                  output filename 

        '''

#mvars:   runname, expdata, ofile, io, ioe, dq, maxpoints, plotflag

        log = self.log
        pgui = self.run_utils.print_gui
        splint = self.splint
        spline = self.spline
        log.debug('in interpolate')

        mvars = self.mvars
        avars = self.avars

        avars.yp1 = 1.0
        avars.ypn = 1.0

        log.debug('calculating splines')

        try:
            spline(avars.y)
            avars.y2 = avars.array2
            spline(avars.z)
            avars.z2 = avars.array2
        except:
            message = 'Failed to interpolate data: is data file corrupt?'
            pgui(message)
            message += '\nstopping here\n\n'
            pgui(message)
            return

        log.debug('back from splines')

        outfile2 = open(avars.interpath + mvars.ofile, 'w')
        outfile3 = open(avars.interpath + 'stn_' + mvars.ofile, 'w')

        avars.io_tally = []
        outfile2.write('%f\t%f\t%f\n' % (0.0, mvars.io, mvars.ioe))
        outfile3.write('%f\t%f\t%f\n' % (0.0, mvars.io, mvars.ioe))
        avars.io_tally.append([0.0, mvars.io, mvars.ioe])
        avars.ux = 0.00
        pgui("\nSignal to noise cutoff value: %s\n" %(str(avars.cutval)))
        for i in range(mvars.maxpoints - 1):
            avars.ux = avars.ux + mvars.dq
            splint(avars.y, avars.y2)
            avars.ny = avars.nyval
            splint(avars.z, avars.z2)
            avars.nz = avars.nyval
            avars.io_tally.append([avars.ux, avars.ny, avars.nz])
            outfile2.write('%f\t%f\t%f\n' % (avars.ux, avars.ny, avars.nz))
            if(avars.ux <= avars.cutval):
                outfile3.write('%f\t%f\t%f\n' % (avars.ux, avars.ny, avars.nz))

        outfile2.close()
        outfile3.close()

        ''' display progress '''

        fraction_done = 1
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        wait = self.wait

        log.debug('in epilogue')

#        print'Interpolated data were written to %s\n' % ('./' + avars.interpath + mvars.ofile)
#        print 'Interpolated data with S/N > 2 were written to %s\n' % ('./' + avars.interpath + 'stn_' + mvars.ofile)
#        print '\ndelta q = %f\t : number of q-points = %i\t : q-range: q = 0 to %f\n' % (mvars.dq, mvars.maxpoints, (mvars.maxpoints - 1) * mvars.dq)
        pgui("\nInterpolated data were written to %s\n" %
                     ('./' + avars.interpath + mvars.ofile))
        pgui("\nInterpolated data with S/N > 2 were written to %s\n\n" %
                    ('./' + avars.interpath + 'stn_' + mvars.ofile))
        pgui("\ndelta q = %f (1/A)\n\nnumber of q-points = %i\n\nq-range: 0 to %f (1/A)\n" %
                    (mvars.dq, mvars.maxpoints, (mvars.maxpoints - 1) * mvars.dq))
    
        if(mvars.plotflag == 1):
            graph = Gnuplot.Gnuplot(debug=1)
            graph.clear()
            graph('set title "Interpolation Results"')
            graph.xlabel('Q (1/A)')
            graph.ylabel('I(Q)')
            graph('set logscale y')
        
            graph.plot(Gnuplot.Data(avars.odata, using='1:2 w p ps 4', title='Original Data'), Gnuplot.Data(avars.io_tally, using='1:2 w lp ps 2',
                title='Interpolated Data'), Gnuplot.Data(avars.cut, title='[I(Q)/(std.dev. I(Q))] < 2', using='1:2:3 w yerrorbars'))

        time.sleep(2)

        if(mvars.plotflag == 1):
            wait('\n')

        self.run_utils.clean_up(log)

        pgui('INTERPOLATE IS DONE')

        pgui("%s \n" % ('=' * 60))
        time.sleep(1.0)

        return

