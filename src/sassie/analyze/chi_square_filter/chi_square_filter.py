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
import sys
import time
import string
import locale
import random
import numpy
import time
import platform
import glob
import sassie.util.basis_to_python as basis_to_python
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
#import chi_helper as chi_helper
import sassie.analyze.chi_square_filter.chi_helper as chi_helper


'''
    CHI_SQUARE_FILTER is the main module that contains the functions
    that are used to compare interpolated experimental data
    to the synthetic data sets generated from structures files
    using the modules in CALCULATE.

    Some filtering options based on Rg are provided.

    This module is called from Chi-Square Filter from the main
    GUI through the graphical_chi_square_filter.py script.
'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'chi_square_filter'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class chi_square_filter_input_variables():

    def __init__(self, parent=None):
        pass

class chi_square_filter():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = chi_square_filter_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.find_best()

        self.epilogue()

        return


    def wait(str=None, prompt='Plot will clear in 10 seconds ...\n'):
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
            time.sleep(1)


#   pgui performs this function
#   def print_failure(message, txtOutput):
#
#       txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#       txtOutput.put(">>>> RUN FAILURE <<<<\n")
#       txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#       txtOutput.put(message)
#
#       return


    def unpack_variables(self,variables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.saspaths = variables['saspaths'][0]
        mvars.sasintfiles = variables['sasintfiles'][0]
        mvars.io = variables['io'][0]
        mvars.number_of_weight_files = variables['number_of_weight_files'][0]
        mvars.basis_string = variables['basis_string'][0]
        mvars.weight_file_names = variables['weight_file_names'][0]
        mvars.sastype = variables['sastype'][0]
        mvars.reduced_x2 = variables['reduced_x2'][0]
        mvars.plotflag = variables['plotflag'][0]

        log.debug(vars(mvars))
        
        return


    def initialization(self):
        '''
        method to prepare for find_best
        '''

#mvars:    runname, saspaths, sasintfiles, io, number_of_weight_files, basis_string, weight_file_names, sastype, reduced_x2, plotflag


        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui
        mvars = self.mvars
        avars = self.avars


        avars.epsilon = 0.00001

        avars.bestworstfile = 'bestworstfile.txt'
        avars.averagefile = 'averagefile.txt'
        avars.x2file = 'x2file.txt'

        avars.sas_paths = [x.strip() for x in mvars.saspaths.split(',')]
        avars.sasintfiles = [x.strip() for x in mvars.sasintfiles.split(',')]

#        print 'sas_paths: ', avars.sas_paths
#        print 'sasintfiles: ', avars.sasintfiles
#        print 'io: ', mvars.io

        base_paths = []
#        print 'sas_paths = ',avars.sas_paths
        for sas_path in avars.sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

#        print 'base_paths = ', base_paths
#        print 'sas_paths = ',avars.sas_paths

        avars.spectra_paths = []
        avars.filter_paths = []

        if(mvars.sastype == 0):
            for base_path in base_paths:
                avars.filter_paths.append(os.path.join(
                    mvars.runname, 'chi_square_filter', base_path))
                avars.spectra_paths.append(os.path.join(
                    mvars.runname, 'chi_square_filter', base_path, 'spectra'))

        else:
            avars.spectra_paths.append(os.path.join(
                mvars.runname, 'chi_square_filter', 'spectra'))
            avars.filter_paths.append(os.path.join(mvars.runname, 'chi_square_filter'))

#        print 'filter_paths: ', avars.filter_paths
#        print 'spectra_paths: ', avars.spectra_paths

        for filter_path in avars.filter_paths:
            direxist = os.path.exists(filter_path)
            if(direxist == 0):
                os.system('mkdir -p ' + filter_path)

        for spectra_path in avars.spectra_paths:
            direxist = os.path.exists(spectra_path)
            if(direxist == 0):
                os.system('mkdir -p ' + spectra_path)

# SK notes:
# Right now, the same weight files would be used for all sas_paths.  The user would choose weight files to cover all sas_paths, i.e., some would be used for one input path and not others.
# We should look at a way to work some kind of global x2 into the filtering,i.e., find cases where the best fit x2 overlaps the files in all of the sas_paths. Perhaps have an advanced option flag to indicate whether you want a global weight file to be calculated. Then you would have to choose which weight file to consider for each sas_path, since the user will likely want to consider different x2 and Rg conditions for different contrasts.  One could then look at the average value of the weight column for each weight file across all the sas_paths.  Then write a new global weight file and flag only the structures where the average was 1.0.
# Another way to do this is to allow the number of weight files, weight file names and basis string to be a list of lists so they can vary with
# sas_path so the user can tailor each weight file to cover the best fit structures at each contrast.  If so, then the "if" statement below needs
# to go inside the loop over sas_paths in the find_best method.

        if (mvars.number_of_weight_files > 0):
#            print 'basis string: ', mvars.basis_string
            avars.python_basis = basis_to_python.clean_up_weight_basis(mvars.basis_string)
            avars.weight_file_names = string.split(mvars.weight_file_names, ',')

        log.debug(vars(mvars))
        log.debug(vars(avars))


    def find_best(self):
        '''
        FIND_BEST is the function to read in variables from GUI input and compare
        interpolated experimental data to synthetic data sets generated from
        structure files using the modules in CALCULATE.

        INPUT:	variable descriptions:

            runname:                used for automatic path generation
            saspaths:               paths to synthetic SAS files
            sasintfiles:            name of interpolated data sets
            io:                     I(0) values
            sastype:                designate scattering calculator (0=sascalc, 1=xtal2sas, 2=cryson, 3=crysol)
            reduced_x2:		        set option for using reduced X2 (0=chi-squared, 1=reduced chi-squared, 2=Pearson's chi-squared, 3=R-factor)
            number_of_weight_files: number of weight files to be generated
            weight_file_names:      names of weight files
            basis_string:           basis strings for weight files, i.e., 'x2<1', '(x2<5) and (Rg<30)'

        OUTPUT:
            txtOutput:              TK handler for output to GUI textbox

            files stored in output file directory:

            bestworstfile.txt:      best and worst sas profiles from single structure
            averagefile.txt:        average profile for the entire ensemble
            x2file.txt:             list of X2 and Rg with each filename
            x2_vs_rg_plot.txt:      Rg and X2 without filename (for plotting)
            sas_spectra_plot.txt:   text file of sas profile plot from Gnuplot

            spectra/                spectra files scaled to I(0)
        '''


        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in find_best')

        mvars = self.mvars
        avars = self.avars

        ttxt = time.asctime(time.gmtime(time.time()))

        st = ''.join(['=' for x in xrange(60)])
        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

# big loop over sas_paths

        count = 0

        print 'sas_paths in find best = ',avars.sas_paths
        print 'filter_paths: ', avars.filter_paths
        print 'spectra_paths: ', avars.spectra_paths
        print 'sasintfile: ', avars.sasintfiles
        

        for sas_path in avars.sas_paths:

            print 'sas_path: ', sas_path

#           read the data file
            ginfile = open(avars.sasintfiles[count], 'r').readlines()
            ngpoints = len(ginfile)
            lowlin = string.split(ginfile[0])
            low = locale.atof(lowlin[0])
            lowplusonelin = string.split(ginfile[1])
            lowplusone = locale.atof(lowplusonelin[0])
            highlin = string.split(ginfile[-1])
            high = locale.atof(highlin[0])

            step = lowplusone - low

            runst = avars.filter_paths[count]
            print 'runst: ', runst

            x2cl = []
            x2ch = []
            rgcl = []
            rgch = []

            gridq = []
            igoal = []
            isigma = []
            for i in xrange(ngpoints):
                lin = string.split(ginfile[i])
                gval = locale.atof(lin[0])
                val = locale.atof(lin[1])
                val2 = locale.atof(lin[2])
                gridq.append(gval)
                igoal.append(val)
                isigma.append(val2)

            goal = numpy.array(igoal)
            sigma = numpy.array(isigma)
            sigma2 = sigma * sigma

            outfile = open(os.path.join(runst, avars.averagefile), 'w')
            x2outfile = open(os.path.join(runst, avars.x2file), 'w')
            bestworstoutfile = open(os.path.join(runst, avars.bestworstfile), 'w')

            nfilestring = 'find ' + sas_path + ' -name "*.log" | grep -c log'
#            print 'nfilestring: ', nfilestring
            nfout = os.popen(nfilestring, 'r').readlines()
            nf = locale.atoi(nfout[0])
            pgui("\n# READING Rg VALUES FROM LOG FILES")
            log = []

            for name in glob.glob(os.path.join(sas_path, '*.log')):
                log.append(name[:-4])

            log.sort()

#            print 'len(log) = ', len(log)

            ns = nf	        # number of trajectory steps
            rgarray = []
            keep = []
            allx2 = []
            for i in xrange(ns):
                inf = open(log[i] + '.log', 'r').readlines()
                ln = len(inf)
                for k in xrange(ln):
                    lin = string.split(inf[k])
                    if(len(lin) > 0):
                        if(mvars.sastype == 0 or mvars.sastype == 1):
                            if(lin[1] == 'Rg'):
                                trg = lin[3]
                                ftrg = locale.atof(trg)
                                keep.append(i)
                                rgarray.append(ftrg)

                        if(mvars.sastype == 2):
                            if(lin[0] == 'Atomic'):
                                trg = lin[3]
                                ftrg = locale.atof(trg)
                                keep.append(i)
                                rgarray.append(ftrg)
                        if(mvars.sastype == 3):
                            if(lin[0] == 'Electron'):
                                trg = lin[3]
                                ftrg = locale.atof(trg)
                                keep.append(i)
                                rgarray.append(ftrg)

            if(len(keep) == 0):
                message = 'did not find any profiles that met the input criterion'
                message += ' :  stopping here'
                pgui(message)
                return

            first = 0
            norm = 1.0  # nh=1  # number of header lines
            for i in xrange(ns):
                if(i in keep):
                    if(mvars.sastype == 0 or mvars.sastype == 1):
                        sinf = open(log[i] + '.iq', 'r').readlines()
                        nh = 0
                    if(mvars.sastype == 2 or mvars.sastype == 3):
                        sinf = open(log[i] + '.int', 'r').readlines()
                        nh = 1
                    sln = len(sinf)
                    jj = 0
                    if(first == 0):
                        allspec = numpy.zeros((len(keep), (sln - nh)), numpy.float)
                        spectra1 = numpy.zeros((sln - nh), numpy.float)
                        spectra2 = numpy.zeros((sln - nh), numpy.float)
                        spectra3 = numpy.zeros((sln - nh), numpy.float)
                        spectra4 = numpy.zeros((sln - nh), numpy.float)
                        qvar = numpy.zeros((sln - nh), numpy.float)
                        first = 1
                        pgui("# READING %i SAS INT FILES" % (ns))
                    for j in xrange(sln):
                        slin = string.split(sinf[j])
                        if(j >= nh):
                            qval = locale.atof(slin[0])
                            qvar[jj] = qval
                            val1 = locale.atof(slin[1])
                            val2 = locale.atof(slin[2])
                            if(mvars.sastype != 0 and mvars.sastype != 1):
                                val3 = locale.atof(slin[3])
                                val4 = locale.atof(slin[4])
                        #   if(mvars.sastype==3 and jj==0):
                        #	    norm=val4
                        #   elif(jj==0):
                            if(jj == 0):
                                norm = val1
                            nval1 = mvars.io[count] * val1 / norm
                            nval2 = mvars.io[count] * val2 / norm
                        #   if(mvars.sastype==2):
                            if(mvars.sastype != 0 and mvars.sastype != 1):
                                nval3 = mvars.io[count] * val3 / norm
                                nval4 = mvars.io[count] * val4 / norm
                                allspec[i][jj] = nval1
                        #   elif(mvars.sastype==3):
                        #	    nval3=mvars.io[count]*val3/norm
                        #	    nval4=mvars.io[count]*val4/norm
                        #	    allspec[i][jj]=nval4
                            else:
                                allspec[i][jj] = nval1

                            spectra1[jj] = spectra1[jj] + nval1
                            spectra2[jj] = spectra2[jj] + nval2
                            if(mvars.sastype != 0 and mvars.sastype != 1):
                                spectra3[jj] = spectra3[jj] + nval3
                                spectra4[jj] = spectra4[jj] + nval4
                            jj = jj + 1

                if(ns > 10):
                    if((i + 1) % (ns / 10) == 0):
                        print_string = '# ' + str(i + 1) + ' FILES READ :\t' + str(int(numpy.ceil((float(i + 1) / float(ns)) * 100.0))) + ' percent done'
                        pgui(print_string)
                        fraction_done = (float(i + 1) / float(ns))
                        progress_string = 'COMPLETED ' + \
                            str(i + 1) + ' of ' + str(ns) + ' : ' + \
                            str(fraction_done * 100.0) + ' % done'
                        pgui('%s\n' % progress_string)
                        report_string = 'STATUS\t' + str(fraction_done)
                        pgui(report_string)

            report_string = 'STATUS\t' + str(1)
            pgui(report_string)

            spectra1 = spectra1 / len(keep)
            spectra2 = spectra2 / len(keep)
            if(mvars.sastype != 0 and mvars.sastype != 1):
                spectra3 = spectra3 / len(keep)
                spectra4 = spectra4 / len(keep)

            outfile.write('# AVERAGE SPECTRA : file generated on ' + ttxt + '\n')
            for k in xrange(sln - nh):
                if(mvars.sastype != 0 and mvars.sastype != 1):
                    outfile.write('%f\t%f\t%f\t%f\t%f\n' % (qvar[k], spectra1[
                                k], spectra2[k], spectra3[k], spectra4[k]))
                else:
                    outfile.write('%f\t%f\t%f\n' %
                                (qvar[k], spectra1[k], spectra2[k]))
            outfile.close()

            x2outfile.write('# file generated on ' + ttxt +
                           '\n# structure, X2, Rg, filename\n')

            numspec = len(keep)
            bq = numpy.zeros(ngpoints, numpy.float)
            best = numpy.zeros(ngpoints, numpy.float)
            worst = numpy.zeros(ngpoints, numpy.float)
            work = numpy.zeros(ngpoints, numpy.float)
            diff = numpy.zeros(ngpoints, numpy.float)
            spec = numpy.zeros((numspec, ngpoints), numpy.float)
            j = 0
            for k in xrange(sln - nh):
                #   print 'k = ',k,' qvar[k] = ',qvar[k]

                if(qvar[k] >= low and qvar[k] <= high + avars.epsilon):
                    #    print 'low = ',low,'high = ',high
                    #    print 'j = ',j
                    bq[j] = qvar[k]
               #   if(mvars.sastype==2):
                    if(mvars.sastype != 0 and mvars.sastype != 1):
                        best[j] = spectra1[k]
               #    elif(mvars.sastype==3):
               #	    best[j]=spectra4[k]
                    else:
                        best[j] = spectra1[k]
                    for i in xrange(numspec):
                        spec[i][j] = allspec[i][k]
                    j += 1
                #   print 'qvar[k] = ',qvar[k]
                #   print 'j = ',j

            if(ngpoints != j):
                log.debug('bug')
                log.debig('ngpoints = %i, j = %i' % (ngpoints,j))
                log.debug('bug')
                log.debug('low = %f, high = %f' % (low,high))
                message = 'number of grid points ' + str(ngpoints) + ' BUG!!'
                message += 'number of j points ' + str(j) + ' BUG!!'
                message += ' :  stopping here'
                pgui(message)

            pgui('# ANALYZING SAS SPECTRA')
            search = 1
            nsteps = 0
            bestspectra = 0
            x2array = []
            best_tally = []
            worst_tally = []
            goal_tally = []
            x2rg_tally = []
            avg_tally = []
            avg_sum = []

            vrarray = []

            for i in xrange(0, nf):
                j = 0
                if(mvars.sastype != 0 and mvars.sastype != 1):
                    suf = '.xiq'
                else:
                    suf = '.ciq'
                thisfile = open(avars.spectra_paths[count] +
                                '/spec_' + str(i + 1) + suf, 'w')
                tspec = numpy.zeros(ngpoints, numpy.float)
                for k in xrange(sln - nh):
                    if(qvar[k] >= low and qvar[k] <= high + avars.epsilon):
                        tspec[j] = allspec[i][j]
                        thisfile.write('%f\t%f\t%f\n' % (qvar[j], tspec[j], 0.0))
                        j = j + 1
                thisfile.close()
                diff = tspec - goal

                vr = numpy.sum(abs(diff / ((tspec + goal) / 2.0)))
                vrarray.append(vr)

                diff2 = diff * diff
                if(mvars.reduced_x2 == 1):     # reduced X2 (assuming n = 0)
                    X2 = numpy.sum(diff2 / sigma2) / (ngpoints - 1)
                elif(mvars.reduced_x2 == 0):
                    X2 = numpy.sum(diff2 / sigma2)  # chi-square distribution
                elif(mvars.reduced_x2 == 2):
                    X2 = numpy.sum(diff2 / goal)  # Pearson's X2 test-statistic
                elif(mvars.reduced_x2 == 3):
                    X2 = numpy.sum(numpy.abs(diff)) / \
                        numpy.sum(numpy.abs(goal))  # R-factor

                x2rg_tally.append([rgarray[i], X2])
                x2outfile.write('%d\t%f\t%f\t%s\t%f\n' %
                            (i + 1, X2, rgarray[i], log[i], vr))
                x2outfile.flush()

                x2array.append(X2)
                if(i == 0):
                    bX2 = X2
                    wX2 = X2
                    pgui('#INITIAL CHI-SQUARED = %f' % (X2))
                    worstfilename = log[i]
                    best = tspec.copy()
                    worst = tspec.copy()
                    avg = tspec.copy()
                    bestspectra = i
                    bestspectra = i + 1
                    worstspectra = i
                    worstspectra = i + 1
                    bestfilename = log[i]
                    for j in xrange(ngpoints):
                        best_tally.append([gridq[j], tspec[j]])
                        worst_tally.append([gridq[j], tspec[j]])
                        goal_tally.append([gridq[j], goal[j]])
                   #    avg_sum.append([gridq[j],tspec[j]])
                        avg_sum.append(tspec[j])
                else:
                    for j in xrange(ngpoints):
                        avg_sum[j] = avg_sum[j] + tspec[j]
                if(X2 < bX2):
                    best = tspec.copy()
                    bestspectra = i + 1
                    pgui('new best:  this X2 = %f' % (X2))
                    best_tally = []
                    for j in xrange(ngpoints):
                        best_tally.append([gridq[j], tspec[j]])
                    bestfilename = log[i]
                    bX2 = X2
                if(X2 > wX2):
                    wX2 = X2
                    pgui('new worst:  this X2 = %f' % (X2))
                    worst = tspec.copy()
                    worst_tally = []
                    worstspectra = i + 1
                    for j in xrange(ngpoints):
                        worst_tally.append([gridq[j], tspec[j]])
                    worstfilename = log[i]

            for j in xrange(ngpoints):
                avg_tally.append([gridq[j], avg_sum[j] / nf])

#            sasplot=open('./'+runst+'sas_spectra_plot.txt','w')
            sasplot = open(os.path.join(runst, 'sas_spectra_plot.txt'), 'w')
            sasplot.write('# q-value, goal, best, worst, average\n')
            for j in xrange(ngpoints):
                sasplot.write('%f\t%f\t%f\t%f\t%f\n' % (gridq[j], goal_tally[j][
                 1], best_tally[j][1], worst_tally[j][1], avg_tally[j][1]))
            sasplot.close()

#            x2plot=open('./'+runst+'x2_vs_rg_plot.txt','w')
            x2plot = open(os.path.join(runst, 'x2_vs_rg_plot.txt'), 'w')
            nnpoints = len(x2rg_tally)
            for j in xrange(nnpoints):
                x2plot.write('%f\t%f\n' % (x2rg_tally[j][0], x2rg_tally[j][1]))
                x2plot.flush()
            x2plot.close()
            x2outfile.close()

            if (mvars.number_of_weight_files > 0):
#                print 'weight file names: ', avars.weight_file_names
#                print 'python_basis: ', avars.python_basis
                chi_helper.write_weight_files(
                    avars.python_basis, x2array, rgarray, avars.weight_file_names, runst)

            pgui('\n# BEST SPECTRA = %i' % (bestspectra))
            pgui('# BEST CHI-SQUARED = %f' % (bX2))
            pgui('# WORST SPECTRA = %i' % (worstspectra))
            pgui('# WORST CHI-SQUARED = %f' % (wX2))
            bestworstoutfile.write('# FILE GENERATED ON ' + ttxt +
                               '\n# gridq,best,goal,diff(best-goal),worst,wdiff(worst-goal)\n')
            bestworstoutfile.write('#bestfilename = %s : bX2 = %f, worstfilename = %s : wX2 = %f\n' % (
                bestfilename, bX2, worstfilename, wX2))
            for k in xrange(ngpoints):
                bestworstoutfile.write('%f\t%f\t%f\t%f\t%f\t%f\n' % (gridq[k], best[k], goal[
                k], best[k] - goal[k], worst[k], worst[k] - goal[k]))
            bestworstoutfile.close()

            pgui("Data stored in directory: %s\n\n" % ('./' + runst))

            pgui("PROCESSED %i SAS FILES:\n\n" % (nf))
            pgui(
                ">> The BEST and WORST SAS spectra are in the file named : %s \n" % (avars.bestworstfile))
            pgui(
                ">> The AVERAGE SAS spectra is in the file named : %s \n" % (avars.averagefile))
            pgui(
                ">> Chi-square, Rg, and filename are in the file named : %s \n" % (avars.x2file))


            path, bfilename = os.path.split(bestfilename)
            path, wfilename = os.path.split(worstfilename)

            pgui("\nBEST SINGLE STRUCTURE IS NUMBER %i WITH X2 = %f : \t spectra: %s\n" % (
                bestspectra, bX2, bfilename))
            pgui("\nWORST SINGLE STRUCTURE IS NUMBER %i WITH X2 = %f: \t spectra: %s\n" % (
                worstspectra, wX2, wfilename))
            pgui("\n%s \n" % (st))

            if(mvars.plotflag == 1):
                wait('\n')
            else:
                time.sleep(2)

            count += 1

# end of big loop over sas_paths

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

        self.run_utils.clean_up(log)

        pgui('CHI_SQUARE_FILTER IS DONE')

        time.sleep(1.0)

        return
