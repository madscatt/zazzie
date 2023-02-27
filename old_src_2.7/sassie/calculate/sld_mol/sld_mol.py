'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D., Hirsh Nanda, Ph.D.

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

import sys
import string
import math
import os
import time
import multiprocessing
import time
import numpy
import scipy.optimize
import sassie.calculate.sld_mol.modified_scipy_interpolate as interpolate
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management
import Gnuplot
import Gnuplot.PlotItems
import Gnuplot.funcutils
#import matplotlib.pyplot as plt

#       SLD_MOL
#
#    11/10/2004    --    initial coding                :    jc
#    12/10/2004    --    prototype for alpha sassie        :    jc
#    12/06/2011    --    re-written with new volumetric code    :    hn
#    12/07/2011    --    added SASMOL support            :    jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    SLD_MOL is the module that calculates the scattering length density profile
    from a dcd/pdb file

    This method compares an experimentally derived SLD profile with heavy atom
    distribution from a pdb or dcd file containing protein structure(s).
    It performs a fit allowing a normalization factor, z-pos & constant shift.
    The SLD profile is convolved with a Gaussian of user defined width to mimic instrument
    resolution and roughness and outputs a text file with frame number and fit_error.
    
'''
if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'sld_mol'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class sld_mol_input_variables():

    def __init__(self, parent=None):
        pass


class sld_mol():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = sld_mol_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization() 

        self.sld_calc()

        self.epilogue()

        return


#    pgui performs this function
#    def print_failure(message,txtOutput):
#
#            txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#            txtOutput.put(">>>> RUN FAILURE <<<<\n")
#            txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#            txtOutput.put(message)
#
#            return



    def unpack_variables(self,variables):
        '''
        method to extract variables into system wise class instance
        '''
        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')    

        mvars.runname = variables['runname'][0]
        mvars.path = variables['path'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.dcdfile = variables['dcdfile'][0]
        mvars.expdatafile = variables['expdatafile'][0]
        mvars.outputfile = variables['outputfile'][0]

        mvars.runtype = variables['runtype'][0]
        mvars.bulk_sld = variables['bulk_sld'][0]
        mvars.xon = variables['xon'][0]

        mvars.num_deut_regions = variables['num_deut_regions'][0]
        mvars.deut_low_res = variables['deut_low_res'][0]
        mvars.deut_high_res = variables['deut_high_res'][0]

        mvars.dbin = variables['dbin'][0]
        mvars.width = variables['width'][0]

        mvars.sldfit = variables['sldfit'][0]
        mvars.sldoffset = variables['sldoffset'][0]

        mvars.zfit0 = variables['zfit0'][0]
        mvars.zfitmin = variables['zfitmin'][0]
        mvars.zfitmax = variables['zfitmax'][0]
        mvars.zevalmin = variables['zevalmin'][0]
        mvars.zevalmax = variables['zevalmax'][0]
        mvars.A0 = variables['A0'][0]
        mvars.Amin = variables['Amin'][0]
        mvars.Amax = variables['Amax'][0]
        mvars.plotflag = variables['plotflag'][0]

        
#mvars: runname, path, pdbfile, dcdfile, expdatafile, outputfile, plotflag, runtype, xon, num_deut_regions, deut_low_res, deut_high_res, dbin, width, sldfit, sldoffset, zfit0, zfitmin, zfitmax, zevalmin, zevalmax, A0, Amin, Amax

#avars:  zexp, interpolated_experimental_sld, interpolated_variance_data, m1, number_of_frames, intype, sldpath, sld_output_files, dcdinputfile, minmaxz, new_dbin, deuterated_residues, resultsfile, radius, volumes, residue_sld_hydrogen, residue_sld_deuterium, residue_sld_electron


    def plot_me(self, plot_type, plot_variables):
        '''
        method to plot the experimental SLD
        '''
        if(plot_type == 'gnuplot'):
            graph = plot_variables[0]
            sld_data = plot_variables[1]
            graph.plot(Gnuplot.Data(sld_data, using='1:2 w p ps 1', title='experimental SLD'))

        elif(plot_type == 'matplotlib'):
            pass

        return


    def read_experimental_sld_file(self, graph):
        '''
        method to read the and interpolate the experimental SLD file
        '''
        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in read_experimental_sld_file')    

        avars.zexp = []
        sld_value = []
        variance_value = []
        inputfile = open(mvars.expdatafile, 'r').readlines()

        pgui('READING EXPERIMENTAL SLD')

        exp_dat = []

        for i in range(len(inputfile)):
            lin = string.split(inputfile[i])
            if(lin[0][0] != "#" and len(lin) >= 2):
#       this is a z-offset and NOT a sld offset as it is named
                avars.zexp.append(float(lin[0]) + mvars.sldoffset)
                sld_value.append(float(lin[1])) 
                exp_dat.append([float(lin[0]) + mvars.sldoffset, float(lin[1])])
                if(len(lin) == 3):
                    variance_value.append(float(lin[2]))

        #interp1d gives us an object (function) to be evaluated at the desired z values later

        avars.interpolated_experimental_sld = interpolate.interp1d(avars.zexp, sld_value, 'cubic')

        if(len(variance_value) == len(sld_value)):
            avars.interpolated_variance_data = interpolate.interp1d(
                avars.zexp, variance_value, 'cubic')
        else:
            avars.interpolated_variance_data = None

        fnew = avars.interpolated_experimental_sld(avars.zexp)  #used to plot experimental data with matplotlib

        if(mvars.plotflag == 1):
            plt.plot(avars.zexp, sld_value, '-', avars.zexp, fnew, 'o')
            raw_input('Please press return to continue...\n')
        elif(mvars.plotflag == 2):
            plot_type = 'gnuplot'
            plot_variables = [graph, exp_dat]
            self.plot_me(plot_type, plot_variables)

        return


    def update_coor(self, residues, resids_mask, frame):

        mvars = self.mvars
        avars = self.avars

        error, new_coor = avars.m1.get_coor_using_mask(frame, resids_mask)

        residues.setCoor(new_coor)

        return


    def print_status(self, step):

        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui     

        fraction_done = (float(step + 1) / float(avars.number_of_frames))
        progress_string = '\nCOMPLETED ' + str(step + 1) + ' of ' + str(
            avars.number_of_frames) + ' : ' + str(fraction_done * 100.0) + ' % done'
        pgui(progress_string)
#        print('%s\n' % progress_string)
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)
#        print report_string

        return


    def residue_com(self, residues, frame):

        mvars = self.mvars
        avars = self.avars    

        all_residue_com = []

        for i in xrange(avars.m1.number_of_resids()):
            resids_mask = avars.m1.resids_mask()[i]
            self.update_coor(residues[i], resids_mask, frame)
            com = residues[i].calccom(0)
            all_residue_com.append(com)

        return numpy.array(all_residue_com)


    def setup_sld_calculation(self):
        '''
            xon: xray (xon=1) or neutron (xon=0)
            deuterated_residues: list of resids for deuterated residues
        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in setup_sld_calculation')    


        amino_acid_properties = avars.m1.amino_acid_sld()

        resid = avars.m1.resid()

        avars.m1.initialize_children()

        residues = avars.m1.init_child("resids")

        number_of_resids = avars.m1.number_of_resids()

        isotopic_residues = numpy.zeros(number_of_resids, numpy.short)

        avars.residue_sld_hydrogen = numpy.zeros(number_of_resids)
        avars.residue_sld_deuterium = numpy.zeros(number_of_resids)
        avars.residue_sld_electron = numpy.zeros(number_of_resids)
        avars.volumes = numpy.zeros(number_of_resids)
        avars.radius = numpy.zeros(number_of_resids)

        hydrogen_scattering_length = -3.7409
        deuterium_scattering_length = 6.671

        first_residue = resid[0]

        for i in avars.deuterated_residues:
            isotopic_residues[i - first_residue] = 1

        for i in range(number_of_resids):
            residue_name = residues[i].resname()[0]
            residue_value = amino_acid_properties[residue_name]
            avars.volumes[i] = residue_value[0]
            if(mvars.xon == 1):
                avars.residue_sld_electron[i] = residue_value[1]
            elif(mvars.xon == 0):
                avars.residue_sld_hydrogen[i] = (residue_value[
                                       2 + isotopic_residues[i]] + residue_value[4] * hydrogen_scattering_length)
                avars.residue_sld_deuterium[i] = (residue_value[
                                        2 + isotopic_residues[i]] + residue_value[4] * deuterium_scattering_length)

        avars.radius = numpy.power(avars.volumes * 3 / (4 * math.pi), 1.0 / 3.0)

        return residues


    def calc_sl(self, frame, residues):
        '''
            This method calculates and returns the scattering lengths and volumes.
            The SLD is calculated in the main routine.
            frames: list of frame numbers to calculate SLD
            xon: xray (xon=1) or neutron (xon=0)
            radius, volumes, residue_sld_hydrogen, residue_sld_deuterium,
            residue_sld_electron, residues all generated in setup_sld_calculation
    '''
        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in calc_sl')    

        log.debug('intype in calc_sl: %s \n' % (avars.intype))
        number_of_resids = avars.m1.number_of_resids()
        log.debug('number of resids: %i\n' % (number_of_resids))
        com_zpos = numpy.array(())  

        log.debug('minmaxz: %s \n' % (avars.minmaxz))
        minz = avars.minmaxz[0]
        maxz = avars.minmaxz[1]

        binnum = int((maxz - minz) / mvars.dbin) + 1

        log.debug('binnum, dbin: %i\t%f\n' % (binnum, mvars.dbin))

        slprofileE = numpy.zeros((binnum))
        slprofileH = numpy.zeros((binnum))
        slprofileD = numpy.zeros((binnum))
        volprofile = numpy.zeros((binnum))

        if(mvars.xon == 0):
            fact = 1.0E-5
        elif(mvars.xon == 1):
            fact = 2.8E-5

        if (avars.intype == 'dcd'):
            avars.m1.read_dcd_step(avars.dcdinputfile, 0)
            all_residue_com = self.residue_com(residues, 0)
        elif (avars.intype == 'pdb'):
            all_residue_com = self.residue_com(residues, frame)

        com_zpos_tmp = numpy.zeros((number_of_resids))

        for i in range(number_of_resids):
            com_zpos_tmp[i] = all_residue_com[i][2]
        
        #getting com_zpos for each resid
        com_zpos = numpy.append(com_zpos, com_zpos_tmp)

        com_zpos = numpy.reshape(com_zpos, (-1, number_of_resids))

        for i in numpy.arange(binnum):
            local_a = mvars.dbin * math.pi * \
                (avars.radius**2 - ((minz + i * mvars.dbin) - com_zpos[0])**2)
            local_a = local_a.clip(min=0.0)
            volprofile[i] = numpy.sum(local_a)
            if(mvars.xon == 1):
                slprofileE[i] = numpy.sum(avars.residue_sld_electron * (local_a / avars.volumes))
            elif(mvars.xon == 0):
                slprofileH[i] = numpy.sum(avars.residue_sld_hydrogen * (local_a / avars.volumes))
                slprofileD[i] = numpy.sum(avars.residue_sld_deuterium * (local_a / avars.volumes))

#        print 'slprofileH: ', slprofileH
#        print 'volprofile: ', volprofile
#        print 'slprofileH*fact: ', slprofileH*fact

        if(mvars.xon == 0):
            return volprofile, slprofileH * fact, slprofileD * fact
        elif(mvars.xon == 1):
            return volprofile, slprofileE * fact, slprofileD * fact


    def smooth(self, x, window_len=11, window='hanning'):
        """smooth the data using a window with requested size.

        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.

        input:
            x: the input signal
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

        output:
            the smoothed signal

        example:

        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)

        see also:

        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter

        TODO: the window parameter could be the window itself if an array instead of a string
        """

        if x.ndim != 1:
            raise ValueError, "smooth only accepts 1 dimension arrays."

        if x.size < window_len:
            raise ValueError, "Input vector needs to be bigger than window size."

        if window_len < 3:
            return x

        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

        s = numpy.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]

        if window == 'flat':  # moving average
            w = numpy.ones(window_len, 'd')
        else:
            w = eval('numpy.' + window + '(window_len)')

        y = numpy.convolve(w / w.sum(), s, mode='valid')
        return y


    def smooth_spectra(self, hvydist):

        mvars = self.mvars
        avars = self.avars

#       the value of 10 has been hardwired here, negating the input width parameter
#        newdist = self.smooth(hvydist, 10, 'blackman')

#       the code below makes use of the input width parameter but doesn't allow values < 10 
#       since widths less then 10 don't produce good smoothing    
        if avars.iwidth < 10:
            avars.iwidth = 10
        newdist = self.smooth(hvydist, avars.iwidth, 'blackman')
        
        return numpy.r_[newdist]


    def calc_error(self, lFitParam, zmaxprot, fnvfprot, fnsldprot):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in calc_error')    

        zfit, Afit = lFitParam    
        err = []     
        zmaxpr = zfit + zmaxprot - mvars.dbin   
        npoints = 0
        for z in numpy.arange(mvars.zevalmin, mvars.zevalmax, mvars.dbin):
            if zfit > z:
#               there is nothing to minimize in this case if sldfit = 1
                if(interpolated_variance_data == None):
                    err.append(abs(avars.interpolated_experimental_sld(z) - mvars.bulk_sld) * 1e6)
                else:
                    err.append(((avars.interpolated_experimental_sld(z) - mvars.bulk_sld)
                           ** 2) / (avars.interpolated_variance_data(z)**2))
                    npoints += 1
            elif zfit <= z and zmaxpr > z:
#               zfit and Afit are minimized if sldfit = 1
                vfprot = fnvfprot(z - zfit); sldprot = fnsldprot(z - zfit)
                SLDz = Afit * (vfprot * sldprot) + Afit * (1 - vfprot) * \
                           mvars.bulk_sld + (1 - Afit) * mvars.bulk_sld
                if(avars.interpolated_variance_data == None):
                    err.append(abs(avars.interpolated_experimental_sld(z) - SLDz) * 1e6)
                else:
                    err.append(((avars.interpolated_experimental_sld(z)-SLDz)**2)/(avars.interpolated_variance_data(z)**2))
                    npoints += 1
            elif zmaxpr <= z:
#               there is nothing to minimize in this case if sldfit = 1
                if(avars.interpolated_variance_data == None):
                    err.append(abs(avars.interpolated_experimental_sld(z) - mvars.bulk_sld) * 1e6)
                else:
                    err.append(((avars.interpolated_experimental_sld(z) - mvars.bulk_sld)
                           ** 2) / (avars.interpolated_variance_data(z)**2))
                    npoints += 1
            else:
                message = 'failed in error calculation\n' + str(zfit)
                pgui(message)

#        print 'npoints, err: ', npoints, err

        if(npoints > 1):
            sum_error = sum(err) / (npoints - 1)
        else:
            sum_error = sum(err)

#        print zfit, Afit, sum_error

        return sum_error


    def plot_and_save_sld_fit(self, param, vfprot, sldprot, graph2, sum_sld_fit, this_frame, fit_error, results):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in plot_and_save_sld')    

        z, A= param
        SLDfit=A*(vfprot*sldprot)+A*(1-vfprot)*mvars.bulk_sld+(1-A)*mvars.bulk_sld
        SLDexp=avars.interpolated_experimental_sld(avars.zexp)

        if(mvars.plotflag == 1):

            plt.figure(3)
            plt.plot(avars.zexp,SLDexp,'-',numpy.arange(len(SLDfit))*avars.new_dbin+z,SLDfit)
            raw_input('Please press return to continue...\n')

#        save experimental data
        exp_data = []    
        for i in xrange(len(avars.zexp)):
            exp_data.append([avars.zexp[i],SLDexp[i]])    
    
        fit_data = []
        avg_sld_fit = []
        for i in xrange(len(SLDfit)):
            this_z = i*avars.new_dbin+z
            fit_data.append([this_z,SLDfit[i]])
            sum_sld_fit[i] += SLDfit[i]
            avg_sld_fit.append([this_z,(sum_sld_fit[i]/this_frame)])
    
#       if runtype==1, get fit error for average SLD
#       NOT YET IMPLEMENTED
        if (mvars.runtype == 1):
            avg_fit_error = 0.0

#       make plots and write files
        if(mvars.plotflag == 2):
            if(mvars.runtype == 0):
                if mvars.sldfit == 0:
                    title0 = 'average SLD'
                elif mvars.sldfit == 1:
                    title0 = 'optimized average SLD'
                graph2.plot(Gnuplot.Data(exp_data, using='1:2 w p ps 1', title='experimental SLD'),
                        Gnuplot.Data(avg_sld_fit, using='1:2 w p ps 1', title=title0))
#           for runtype = 1, each SLD and a running average are plotted as they are calculated
#           takes some time if there are a large number of frames
#           final best, worst and average structures are plotted instead
#           elif(mvars.runtype == 1):
#               if mvars.sldfit == 0:
#                   title0 = 'SLD for structure '+str(this_frame)
#                   title1 = 'average SLD for structures 1 - '+str(this_frame)
#               elif mvars.sldfit == 1:
#                   title0 = 'optimized SLD for structure '+str(this_frame)
#                   title1 = 'average of optimized SLDs (1 - '+str(this_frame)+')'
#               graph2.plot(Gnuplot.Data(exp_data, using='1:2 w p ps 1', title='experimental SLD'),Gnuplot.Data(fit_data, using='1:2 w p ps 1', title=title0), Gnuplot.Data(avg_sld_fit, using='1:2 w p ps 1', title=title1))

        if(mvars.runtype == 1):
            outputfile = open(avars.sld_output_files[this_frame - 1], 'w')
        average_outputfile = open(avars.sldpath + '/average_sld.txt', 'w')


#       write header for average_sld.txt file
        if(mvars.runtype == 0):
            if mvars.sldfit == 0:
                st = "#  Z POSITION = " + str(results[0][0]) + " : SURFACE COVERAGE = " + str(results[0][1]) + " : CHI_SQUARED = %0.4f" % (fit_error) + "\n"
            elif mvars.sldfit == 1:
                st = "#  OPTIMIZED Z POSITION = " + str(results[0][0]) + " : OPTIMIZED SURFACE COVERAGE = " + str(results[0][1]) + " : CHI_SQUARED = %0.4f" % (fit_error) + "\n"      

#NOTE:  for runtype == 1 and sldfit == 1, the optimized z position and surface coverage is the value for the last structure since the values will differ for each structure in this case.  
        elif(mvars.runtype == 1):
            if mvars.sldfit == 0:
                st = "#  STRUCTURES = 1 to " + str(this_frame) + ": Z POSITION = " + str(results[0][0]) + " : SURFACE COVERAGE = " + str(results[0][1]) + " : CHI_SQUARED = %0.4f"  % (avg_fit_error) + "\n"
            elif mvars.sldfit == 1:
                st = "#  STRUCTURES = 1 to " + str(this_frame) + ": OPTIMIZED Z POSITION (last structure) = " + str(results[0][0]) + " : OPTIMIZED SURFACE COVERAGE (last structure) = " + str(results[0][1]) + " : CHI_SQUARED = %0.4f"  % (avg_fit_error) + "\n"              


        average_outputfile.write(st)

#       write header for calculated individual SLD files (runtype == 1 only)
        if(mvars.runtype == 1):
            if mvars.sldfit == 0:
                st1 = "#  STRUCTURE = " + str(this_frame) + ": Z POSITION = " + str(results[0][
                             0]) + " : SURFACE COVERAGE = " + str(results[0][1]) + " : CHI_SQUARED = %0.4f"  % (fit_error) + "\n"
            elif mvars.sldfit == 1:
                st1 = "#  STRUCTURE = " + str(this_frame) + ": OPTIMIZED Z POSITION = " + str(results[0][
                             0]) + " : OPTIMIZED SURFACE COVERAGE = " + str(results[0][1]) + " : CHI_SQUARED = %0.4f"  % (fit_error) + "\n"                         

            outputfile.write(st1)

#       write column headings for calculated and experimental SLD files
        st2 = "# z (angstroms) : SLD\n"
        if(mvars.runtype == 1):
            outputfile.write(st2)
        average_outputfile.write(st2)

        for i in xrange(len(SLDfit)):
            if(mvars.runtype == 1):
                outputfile.write("%f\t%e\n" % (fit_data[i][0],fit_data[i][1]))
            average_outputfile.write("%f\t%e\n" % (avg_sld_fit[i][0],avg_sld_fit[i][1]))

        if(mvars.runtype == 1):
            outputfile.close()
        average_outputfile.close()

        return

    def plot_and_save_bestworstavg(self, graph1):
        '''
        finds SLD with the best chi-squared and the worst chi-squared
        only called if runtype = 1
        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in plot_and_save_bestworstavg')    
    
        infile=open(avars.resultsfile,'r').readlines()
        npts = 0
        bestnum = 0
        worstnum = 0
        for i in xrange(len(infile)):
            lin = string.split(infile[i])
#            print 'line: ', lin
#            print lin[0]
#            print lin[3]
            if(lin[0][0] != "#" and len(lin) >= 2):
                number=int(lin[0])
                errval=float(lin[3])
                if(npts == 0):
                    best = errval
                    bestnum = number
                    worst = errval
                    worstnum = number
#                    print 'initial best: ', number, errval
#                    print 'initial worst: ', number, errval
                npts += 1
#                print 'number, npts, errval: ', number,errval
                if(errval < best):
                    best = errval
                    bestnum = npts
#                    print 'new best: ', number, errval
                elif(errval > worst):
                    worst = errval
                    worstnum = npts
#                    print 'new worst: ', number, errval    
            
#        print 'final best: ', bestnum, best
#        print 'final worst: ', worstnum, worst           

        bestfile = avars.sld_output_files[bestnum-1]
#        print 'bestfile: ', bestfile
        worstfile = avars.sld_output_files[worstnum-1]
#        print 'worstfile: ', worstfile
        avgfile = avars.sldpath + '/average_sld.txt'
        expfile = mvars.expdatafile

        bfile=open(bestfile, 'r').readlines()
        zbest = []
        sldbest = []
        bestdata = []
        bnpts = 0
        for i in xrange(len(bfile)):
            blin = string.split(bfile[i])
#            print 'blin: ', blin    
            if(blin[0][0] != "#" and len(blin) >= 2):
                bnpts += 1
                zbest.append(float(blin[0]))
                sldbest.append(float(blin[1]))
#        print 'bnpts: ', bnpts
#        print 'zbest: ', zbest
#        print 'sldbest: ', sldbest
        for i in xrange(bnpts):
            bestdata.append([zbest[i], sldbest[i]])
        
        wfile=open(worstfile, 'r').readlines()
        zworst = []
        sldworst = []
        worstdata = []
        wnpts = 0
        for i in xrange(len(wfile)):
            wlin = string.split(wfile[i])
#            print 'wlin: ', wlin    
            if(wlin[0][0] != "#" and len(wlin) >= 2):
                wnpts += 1
                zworst.append(float(wlin[0]))
                sldworst.append(float(wlin[1]))
#        print 'wnpts: ', wnpts
#        print 'zworst: ', zworst
#        print 'sldworst: ', sldworst
        for i in xrange(wnpts):
            worstdata.append([zworst[i], sldworst[i]])


        afile=open(avgfile, 'r').readlines()
        zavg = []
        sldavg = []
        avgdata = []
        anpts = 0
        for i in xrange(len(afile)):
            alin = string.split(afile[i])
#           print 'alin: ', alin    
            if(alin[0][0] != "#" and len(alin) >= 2):
                anpts += 1
                zavg.append(float(alin[0]))
                sldavg.append(float(alin[1]))
#        print 'anpts: ', anpts
#        print 'zavg: ', zavg
#        print 'sldavg: ', sldavg
        for i in xrange(anpts):
            avgdata.append([zavg[i], sldavg[i]])
                
        efile=open(expfile, 'r').readlines()
        z_exp = []
        sldexp = []
        expdata = []
        enpts = 0
        for i in xrange(len(efile)):
            elin = string.split(efile[i])
#            print 'elin: ', elin    
            if(elin[0][0] != "#" and len(elin) >= 2):
                enpts += 1
                z_exp.append(float(elin[0])+mvars.sldoffset)
                sldexp.append(float(elin[1]))
#        print 'enpts: ', enpts
#        print 'z_exp: ', z_exp
#        print 'sldexp: ', sldexp
        for i in xrange(enpts):
            expdata.append([z_exp[i], sldexp[i]])

        if(mvars.plotflag == 2):
            graph1.plot(Gnuplot.Data(expdata, using='1:2 w p ps 1', title='experimental'), Gnuplot.Data(
                bestdata, using='1:2 w p ps 1', title='best'), Gnuplot.Data(worstdata, using='1:2 w p ps 1', title='worst'), Gnuplot.Data(avgdata, using='1:2 w p ps 1', title='average'))


        bestworst_file = open(avars.sldpath + '/bestworstfile.txt', 'w')
        st = '\n#BEST FILE = '+ bestfile +' : CHI-SQUARED = '+str(best)+'\n'
        bestworst_file.write(st)
        st1 = '\n#WORST FILE = '+ worstfile +' : CHI-SQUARED = '+str(worst)+'\n'
        bestworst_file.write(st1)

        pgui(st)
        pgui(st1)

#        print ('%s\n' % st)
#        print ('%s\n' % st1)

        return 

    def initialization(self):
        '''
        method to prepare for sld_calc
        '''    

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in initialization')    
    

        avars.m1 = sasmol.SasMol(0)
        avars.m1.read_pdb(mvars.path+mvars.pdbfile)

#        print 'calculating minmax for all frames'

        try:
            if(mvars.dcdfile[-3:] == 'dcd'):
                this_dcdinputfile = avars.m1.open_dcd_read(mvars.path+mvars.dcdfile)
                avars.number_of_frames = this_dcdinputfile[2] 
                avars.intype = 'dcd'
            elif(mvars.dcdfile[-3:] == 'pdb'):
                avars.m1.read_pdb(mvars.path+mvars.dcdfile)
                avars.number_of_frames = avars.m1.number_of_frames()
#                print 'number of pdb frames: ', avars.number_of_frames
                avars.intype = 'pdb'
        except:
            message = 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
            message += ' :  stopping here'
#           print message
            pgui(message)
            sys.exit()

#        print 'intype: ', avars.intype

        for i in xrange(avars.number_of_frames):
            if(avars.intype == 'dcd'):
                avars.m1.read_dcd_step(this_dcdinputfile, i)
                minmax = avars.m1.calcminmax_frame(0)
            elif(avars.intype == 'pdb'):
                minmax = avars.m1.calcminmax_frame(i)

        avars.minmaxz = [minmax[0][2],minmax[1][2]]
        log.debug('minmaxz in main: %s \n' % (str(avars.minmaxz)))

        binnum = int((minmax[1][2]-minmax[0][2])/mvars.dbin)+1

#        print 'binnum in main: ', binnum
#        print 'number_of_frames = ',avars.number_of_frames

        if(avars.intype == 'dcd'):
            avars.m1.close_dcd_read(this_dcdinputfile[0])   

        direxist = os.path.exists(mvars.runname)
        if(direxist == 0):
            os.system('mkdir -p ' + mvars.runname + '/')

        avars.sldpath = os.path.join(mvars.runname, 'sld_mol')
        direxist = os.path.exists(avars.sldpath)
        if(direxist == 0):
            os.system('mkdir -p ' + avars.sldpath)

        avars.sld_output_files = []

        for i in range(avars.number_of_frames):
            nst = str(i + 1).zfill(5)
            avars.sld_output_files.append(avars.sldpath + '/sldfile_' + nst + '.txt')

#        print 'sldpath: ', avars.sldpath
#        print 'sld_output_files: ', avars.sld_output_files



        return 

    def get_deuterated_atoms(self):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in get_deuterated_atoms')    

        avars.deuterated_residues = []
        for i in xrange((mvars.num_deut_regions)):
            for k in xrange(mvars.deut_low_res[i],mvars.deut_high_res[i]+1):
                avars.deuterated_residues.append(k)    
        return


    def sld_calc(self):
        '''
        SLD_MOL is the module that calculates the scattering length density profile
        from a dcd/pdb file

        This method compares an experimentally derived SLD profile with heavy atom
        distribution from a pdb or dcd file containing protein structure(s).
        It performs a fit allowing a normalization factor, z-pos & constant shift.
        The SLD profile is convolved with a Gaussian of user defined width to mimic instrument
        resolution and roughness and outputs a text file with frame number and fit_error.
    
        INPUT:  variable descriptions:

            runname:            project name                                          
            pdbfile:            reference PDB file
            dcdfile:            input filename (DCD or PDB)
            expdatafile:        experimental SLD data file name
            outputfile:         output file name 
            runtype:            0 = average SLD over all structures, 1 = best fit SLD for each individual structure
            bulk_sld:           SLD for bulk solvent
            xon:                scattering type: 0 = neutron, 1 = x-ray
            num_deut_regions:   number of fully deuterated regions in molecule
            deut_low_res:       low residue number(s) for deuterated region(s)
            deut_high_res:      high residue number(s) for deuterated region(s)
            sldfit:             0 = no optimization of z0 and A0, 1 = optimize z0 and A0
            sldoffset:          offset correction to experimental SLD
            dbin:               bin width for z values
            width:              SLD width for Gaussian smoothing
            zfit0:              z0 value (anchoring position) for calculated SLDs initial estimate
            zfitmin:            minimum z value used during optimization
            zfitmax:            maximum z value used during optimization
            zevalmin:           minimum z to evaluate experimental SLD during error calculation
            zevalmax:           maximum z to evaluate experimental SLD during error calculation
            A0:                 fraction surface coverage for calculated SLDs initial estimate
            Amin:               minimum A0 value used during optimization
            Amax:               maximum A0 value used during optimization


        OUTPUT:

            files stored in ./"runname"/sld_mol/ directory:

            outputfile:         output file containing z0, A0 and fit-error for each calculated SLD
            average_sld.txt:    file containing average SLD
            sldfile*.txt:       files containing individual SLDs for each frame (runtype = 1 only)
            bestworstfile.txt:  file containing the filenames and error values for the best- and worst-fitting structures

        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui        
        log.debug('in main sld_calc')    

        st =''.join(['=' for x in xrange(60)])
        ttxt=time.asctime( time.gmtime( time.time() ) ) 
        pgui("\n%s \n" %(st))
        pgui("DATA FROM RUN: %s \n\n" %(ttxt))

        graph = 0
        graph1 = 0
        graph2 = 0

        if(mvars.plotflag == 2):
            plot_type = 'gnuplot'
            graph = Gnuplot.Gnuplot(debug=1)
            graph.clear()
            graph('set title "Experimental SLD"')
            graph.xlabel('Z position (Angstroms)')
            graph.ylabel('Scattering Length Density')

            graph1 = Gnuplot.Gnuplot(debug=1)
            graph1.clear()
            graph1('set title "Experimental and Calculated SLDs"')
            graph1.xlabel('Z position (Angstroms)')
            graph1.ylabel('Scattering Length Density')
    
            graph2 = Gnuplot.Gnuplot(debug=1)
            graph2.clear()
            graph2('set title "Experimental and Calculated SLDs"')
            graph2.xlabel('Z position (Angstroms)')
            graph2.ylabel('Scattering Length Density')

        avars.iwidth = int(mvars.width/mvars.dbin)

        log.debug('iwidth: %i \n' %(avars.iwidth))

        SLDh2o = -0.567E-6
        SLDd2o = 6.33E-6

        deut_ratio = (mvars.bulk_sld - SLDh2o) / (SLDd2o - SLDh2o)

        log.debug('deut_ratio: %f \n' % (deut_ratio))

        self.get_deuterated_atoms()

        pgui('deuterated_residues: %s\n' % (str(avars.deuterated_residues)))

        fitparam0 = [mvars.zfit0,mvars.A0]
        parambounds = [[mvars.zfitmin, mvars.zfitmax], [mvars.Amin, mvars.Amax]]

        log.debug('fitparam0,parambounds: %s%s\n' % (str(fitparam0), str(parambounds)))
    
        self.read_experimental_sld_file(graph)

        outfile=open(avars.sldpath+'/'+mvars.outputfile,'w')
        if(mvars.sldfit == 0):
            st = "#  STRUCTURE : Z POSITION : SURFACE COVERAGE : CHI-SQUARED\n"
        elif(mvars.sldfit == 1):
            st = "#  STRUCTURE : OPTIMIZED Z POSITION : OPTIMIZED SURFACE COVERAGE : CHI-SQUARED\n"
        outfile.write(st)

        residues = self.setup_sld_calculation()

        if(mvars.runtype == 0):

            pgui("CALCULATING AVERAGE SLD FOR ENSEMBLE\n")
#            print 'CALCULATING AVERAGE SLD FOR ENSEMBLE'

            sum_flag = 0

#           Calculate initial individual SLDs and average SLD 
            if (avars.intype == 'dcd'):
                avars.dcdinputfile = avars.m1.open_dcd_read(mvars.path+mvars.dcdfile)
            elif (avars.intype == 'pdb'):
                avars.dcdinputfile = 'None'

            for i in xrange(avars.number_of_frames):
            
                if(mvars.xon==0):
                    volprof,slprofH,slprofD = self.calc_sl(i,residues)
                    slprof = (1.0-deut_ratio)*slprofH + deut_ratio*slprofD
                else:
                    volprof,slprof,slprofD = self.calc_sl(i,residues)
        
                sldprof = numpy.zeros(len(volprof),numpy.float)
                for ii in xrange(len(volprof)):
                    if(volprof[ii] > 0):
                        sldprof[ii] = slprof[ii]/volprof[ii]
        
                if(sum_flag == 0):
                    sumsldprof = numpy.zeros(len(volprof),numpy.float)
                    sum_flag = 1

                sumsldprof += sldprof 

            avgsldprof = numpy.zeros(len(sumsldprof),numpy.float)
            avgsldprof = sumsldprof/avars.number_of_frames        

#           smoothing average SLD, z and volume profiles
            avgsldprofsmooth = self.smooth_spectra(avgsldprof)
            zprofsmooth = numpy.arange(len(avgsldprofsmooth))*mvars.dbin
            zmaxprofsmooth = len(avgsldprofsmooth)*mvars.dbin
            volprofsmooth = self.smooth_spectra(volprof)
            vfprofsmooth = volprofsmooth/numpy.max(volprofsmooth)
#           print 'zmaxprofsmooth: ', zmaxprofsmooth

#           interpolation function for average SLD and volume
            fnsldprofsmooth = interpolate.interp1d(zprofsmooth,avgsldprofsmooth,'cubic')
            fnvfprofsmooth = interpolate.interp1d(zprofsmooth,vfprofsmooth,'cubic')

            if(mvars.sldfit == 1):    
                pgui("OPTIMIZING AVERAGE SLD\n")
#                print 'OPTIMIZING AVERAGE SLD'
                results = scipy.optimize.fmin_tnc(self.calc_error, fitparam0, fprime=None, args=(zmaxprofsmooth, fnvfprofsmooth, fnsldprofsmooth), approx_grad=1, bounds=parambounds)
        
            else:
                results=[]
                results.append(fitparam0)

#            print 'results: ', results[0]   

            fit_error = self.calc_error(results[0],zmaxprofsmooth, fnvfprofsmooth, fnsldprofsmooth)
        
            this_frame = 1

#            print results[0][0], results[0][1], fit_error
            outfile.write('%s\t%f\t%f\t%0.4f\n' %(this_frame, results[0][0], results[0][1], fit_error))

            avars.new_dbin = float(len(volprof))*mvars.dbin/float(len(volprofsmooth))

#            print 'dbin, new_dbin: ', mvars.dbin, avars.new_dbin
    
            sum_sld_fit = numpy.zeros(len(vfprofsmooth),numpy.float)
            self.plot_and_save_sld_fit(results[0], vfprofsmooth, avgsldprofsmooth, graph2, sum_sld_fit, this_frame, fit_error, results)

            i = avars.number_of_frames - 1
            self.print_status(i)

        elif(mvars.runtype == 1):
    
            sum_flag = 0

            pgui("CALCULATING SLD FOR EACH STRUCTURE\n")
#            print 'CALCULATING SLD FOR EACH STRUCTURE'

            if (avars.intype == 'dcd'):
                avars.dcdinputfile = avars.m1.open_dcd_read(mvars.path+mvars.dcdfile)
            elif (avars.intype == 'pdb'):
                avars.dcdinputfile = 'None'

            for i in xrange(avars.number_of_frames):
            
                if(mvars.xon==0):
                    volprof,slprofH,slprofD = self.calc_sl(i,residues)
                    slprof = (1.0-deut_ratio)*slprofH + deut_ratio*slprofD
                else:
                    volprof,slprof,slprofD = self.calc_sl(i,residues)
        
                sldprof = numpy.zeros(len(volprof),numpy.float)
                for ii in xrange(len(volprof)):
                    if(volprof[ii] > 0):
                        sldprof[ii] = slprof[ii]/volprof[ii]

                volprofsmooth = self.smooth_spectra(volprof)
                sldprofsmooth = self.smooth_spectra(sldprof)

                vfprofsmooth = volprofsmooth/numpy.max(volprofsmooth)
                zprofsmooth = numpy.arange(len(sldprofsmooth))*mvars.dbin
                zmaxprofsmooth = len(sldprofsmooth)*mvars.dbin
#                print 'zmaxprofsmooth: ', zmaxprofsmooth            
        
                fnsldprofsmooth = interpolate.interp1d(zprofsmooth,sldprofsmooth,'cubic')
                fnvfprofsmooth = interpolate.interp1d(zprofsmooth,vfprofsmooth,'cubic')
        
                if(mvars.sldfit == 1):    
                    results = scipy.optimize.fmin_tnc(self.calc_error, fitparam0, fprime=None, args=(zmaxprofsmooth, fnvfprofsmooth, fnsldprofsmooth), approx_grad=1, bounds=parambounds)
                else:
                    results=[]
                    results.append(fitparam0)
        
                fit_error = self.calc_error(results[0],zmaxprofsmooth, fnvfprofsmooth, fnsldprofsmooth)
            
                this_frame = i+1
#                print results[0][0], results[0][1], fit_error
                outfile.write('%s\t%f\t%f\t%0.4f\n' %(this_frame, results[0][0], results[0][1], fit_error))
#                outfile.write(`this_frame`+'\t'+`results[0][0]`+'\t'+`results[0][1]`+'\t'+`fit_error`+'\n')
    
                if(sum_flag == 0):
                    sum_sld_fit = numpy.zeros(len(vfprofsmooth),numpy.float)
                    sum_flag = 1
        
                avars.new_dbin = float(len(volprof))*mvars.dbin/float(len(volprofsmooth))
        
                self.plot_and_save_sld_fit(results[0], vfprofsmooth, sldprofsmooth, graph2, sum_sld_fit, this_frame, fit_error, results)

                if(((i+1)%(float(avars.number_of_frames)/100.0)==0 or (i<10))):   
                    self.print_status(i)
    
        outfile.close()   

        if(mvars.runtype == 1):
#           plot best, worst, average slds and write bestworstfile.txt
            avars.resultsfile = avars.sldpath+'/'+mvars.outputfile
            self.plot_and_save_bestworstavg(graph1)

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
        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        pgui("\nscattering length density profile(s) and statistics are saved in %s directory\n\n" % ('./'+avars.sldpath))

        st=''.join(['=' for x in xrange(60)])
        pgui("\n%s \n" %(st))
        pgui('SLD IS DONE')
        time.sleep(2)



        return        


