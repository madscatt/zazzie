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
import string
import locale
import bisect
import time
import numpy
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.analyze.cube as cube
import sassie.analyze.renorm as renorm
import sassie.interface.input_filter as input_filter

#       DENSITY_PLOT
#
#       11/20/2005       --      initial coding         :	jc
#	    03/15/2012	 --	 added segment based cube files :	jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    DENSITY_PLOT is the module that contains the functions
    that are used to compare interpolated experimental data
    to the synthetic data sets generated from structures files
    that are compared using the modules in CHI_SQUARE_FILTER. 
	This is done by populating voxels in three-dimensions based
	on occupancy and/or weighting functions representing quality
	of fit (X2) or Rg ranges.

    This module is called from Density Plot from the main
    GUI through the graphical_density.py script.

	This module calls to C / Python extension modules to speed up
	the binning (see: cube.c and renorm.c).

'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'density_plot'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class density_plot_input_variables():

    def __init__(self, parent=None):
        pass

class density_plot():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, segvariables, txtOutput):

        self.mvars = module_variables()

        self.avars = density_plot_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.process_input_variables(segvariables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.density(input_variables, segvariables, txtOutput)

        self.epilogue()

        return

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
        mvars.dcdfile = variables['dcdfile'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.ofile = variables['ofile'][0]
        mvars.xlength = variables['xlength'][0]
        mvars.gridsp = variables['gridsp'][0]
        mvars.ylength = variables['ylength'][0]
        mvars.save_occupancy = variables['save_occupancy'][0]
        mvars.zlength = variables['zlength'][0]
        mvars.nsegments = variables['nsegments'][0]
        mvars.equalweights = variables['equalweights'][0]
        mvars.weightsfile = variables['weightsfile'][0]

        log.debug(vars(mvars))

        return


    def process_input_variables(self,segvariables):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in process_input_variables')

        log.debug('segvariables: %s' %(segvariables))

        mvars.allsnumregions = []
        mvars.allslow = []
        mvars.allshigh = []
        mvars.allsbasis = []
        mvars.allsname = []

        log.debug('len(segvariables): %i' % len(segvariables))
        log.debug('nsegments: %i' % (mvars.nsegments))

        for i in range(len(segvariables)):
            mvars.allsnumregions.append(segvariables[i][0])
            mvars.allslow.append(segvariables[i][1])
            mvars.allshigh.append(segvariables[i][2])
            mvars.allsbasis.append(segvariables[i][3])
            mvars.allsname.append(segvariables[i][4])

        avars.anregions = []

        for i in range(len(mvars.allsnumregions)):
            nr = locale.atoi(mvars.allsnumregions[i])
            avars.anregions.append(nr)

        avars.alow = []
        avars.ahigh = []

        for i in range(len(mvars.allslow)):
            linlow = string.split(mvars.allslow[i], ',')
            linhigh = string.split(mvars.allshigh[i], ',')
            rlow = []
            rhigh = []
            for k in range(len(linlow)):
                tlow = locale.atoi(linlow[k])
                thigh = locale.atoi(linhigh[k])
                rlow.append(tlow)
                rhigh.append(thigh)
            avars.alow.append(rlow)
            avars.ahigh.append(rhigh)

        log.debug(vars(mvars))
        log.debug(vars(avars))

        return



    def write_cube_header(self,number_atoms,outfile, remark=''):
        '''
        method to write cube file header
        '''
        log = self.log
        log.debug('in write_cube_header')

        mvars = self.mvars
        avars = self.avars

        st1 = remark + 'SASSIE/GAUSSIAN CUBE FILE\n'
        st2 = remark + 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n'
        outfile.write(remark + 'calpha_' + st1 + st2)
        outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                    (number_atoms, avars.xmin * avars.ang2au, avars.ymin * avars.ang2au, avars.zmin * avars.ang2au))
        outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                    (avars.nxgp, mvars.gridsp * avars.ang2au, 0.00, 0.00))
        outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                    (avars.nygp, 0.00, mvars.gridsp * avars.ang2au, 0.00))
        outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                    (avars.nzgp, 0.00, 0.00, mvars.gridsp * avars.ang2au))

        return            


    def write_cube(self, fltr, abin, number_atoms, outfile):
        '''
        method to write cube data
        '''
        log = self.log
        log.debug('in write_cube')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        coor = avars.m1.coor()
        name = avars.m1.name()
        resid = avars.m1.resid()
        segname = avars.m1.segname()

        self.write_cube_header(number_atoms, outfile)

        for natom in range(len(name)):
            if eval(fltr):
                outfile.write('%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                            0, natom, 0] * avars.ang2au, coor[0, natom, 1] * avars.ang2au, coor[0, natom, 2] * avars.ang2au))

        for i in range(avars.nxgp):

            if((((i + 1) % avars.nxgp) / 10.0) == 0 or (avars.nxgp < 10)):
                fraction_done = (float(i + 1) / float(avars.nxgp))
                progress_string = '\nWROTE ' + str(i + 1) + ' of ' + str(
                    avars.nxgp) + ' grid points : ' + str(fraction_done * 100.0) + ' % done'
                pgui('%s\n' % progress_string)
                report_string = 'STATUS\t' + str(fraction_done)
                pgui(report_string)

            for j in range(avars.nygp):
                for k in range(avars.nzgp):
                    outfile.write('%lf ' % (abin[i][j][k] * avars.normcbin))
                    if((k % 6) == 5):
                        outfile.write('\n')
                outfile.write('\n')

        filename = str(outfile).split("'")[1]
        pgui("Wrote the density data to : %s\n\n" % filename)

        return


#    def write_occupancy(self,m, nxgp, nygp, nzgp, nsegments, anregions, cbin, sbin, arbin, number_atoms_full, number_atoms_segments, #number_atoms_regions,   allsbasis, allsname, alow, ahigh, xmin, ymin, zmin, ang2au, gridsp, outfile, txtOutput):
    def write_occupancy(self, outfile):
        '''
        method to write occupancy data
        '''
        log = self.log
        log.debug('in write_occupancy')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        coor = avars.m1.coor()
        name = avars.m1.name()
        resid = avars.m1.resid()
        segname = avars.m1.segname()

        fltr = ''
        for i in xrange(mvars.nsegments):
            fltr += 'segname[natom] == "' + mvars.allsname[i] + \
                '" and name[natom] == "' + mvars.allsbasis[i] + '"'
            if i < (mvars.nsegments - 1):
                fltr += ' or '
        outfile.write('\n############################\n#Complete\n')
#        write_cube_header(number_atoms_full, nxgp, nygp, nzgp,
#                        xmin, ymin, zmin, ang2au, gridsp, outfile, remark='#')
        log.debug('fltr 1 in write_occupancy: %s' % (fltr))
        self.write_cube_header(avars.number_atoms_full, outfile, remark='#')
        for natom in range(len(name)):
            if eval(fltr):
                outfile.write('#%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                            0, natom, 0] * avars.ang2au, coor[0, natom, 1] * avars.ang2au, coor[0, natom, 2] * avars.ang2au))

        for ii in xrange(mvars.nsegments):
            outfile.write('\n############################\n#Seg_%d  \n' % ii)
            fltr = "name[natom]=='" + mvars.allsbasis[ii] + \
                "' and segname[natom]=='" + mvars.allsname[ii] + "'"
#            write_cube_header(number_atoms_segments[
#                            ii], nxgp, nygp, nzgp, xmin, ymin, zmin, ang2au, gridsp, outfile, remark='#')
            log.debug('fltr 2 in write_occupancy: %s' % (fltr))
            self.write_cube_header(avars.number_atoms_segments[ii], outfile, remark='#')                            
            for natom in range(len(name)):
                if eval(fltr):
                    outfile.write('#%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                                0, natom, 0] * avars.ang2au, coor[0, natom, 1] * avars.ang2au, coor[0, natom, 2] * avars.ang2au))
            for jj in range(avars.anregions[ii]):
                outfile.write(
                    '\n############################\n#Seg_%d_Reg_%d  \n' % (ii, jj))
#                write_cube_header(number_atoms_regions[ii][
#                                jj], nxgp, nygp, nzgp, xmin, ymin, zmin, ang2au, gridsp, outfile, remark='#')
                self.write_cube_header(avars.anumber_atoms_regions[ii][jj], outfile, remark='#')
                fltr = "name[natom]=='" + mvars.allsbasis[ii] + "' and segname[natom]=='" + mvars.allsname[ii] + \
                    "' and (resid[natom]>=" + str(avars.alow[ii][jj]) + \
                    " and resid[natom]<=" + str(avars.ahigh[ii][jj]) + ")"
                log.debug('fltr 3 in write_occupancy: %s' % (fltr))
                for natom in range(len(name)):
                    if eval(fltr):
                        outfile.write('#%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                                    0, natom, 0] * avars.ang2au, coor[0, natom, 1] * avars.ang2au, coor[0, natom, 2] * avars.ang2au))

        outfile.write('\n############################\n')
        outfile.write('# X       Y        Z        Complete  ')
        for ii in xrange(mvars.nsegments):
            outfile.write('Seg_%d  ' % ii)
            for jj in range(avars.anregions[ii]):
                outfile.write('Seg_%d_Reg_%d ' % (ii, jj))
        outfile.write('\n')

        for i in range(avars.nxgp):

            if((((i + 1) % avars.nxgp) / 10.0) == 0 or (avars.nxgp < 10)):
                fraction_done = (float(i + 1) / float(avars.nxgp))
                progress_string = 'WROTE ' + \
                    str(i + 1) + ' of ' + str(avars.nxgp) + ' : ' + \
                    str(fraction_done * 100.0) + ' % done'
                pgui('%s\n' % progress_string)
                report_string = 'STATUS\t' + str(fraction_done)
                pgui(report_string)

            for j in range(avars.nygp):
                for k in range(avars.nzgp):
                    outfile.write('%8.3f %8.3f %8.3f ' % (
                        avars.xmin + i * mvars.gridsp, avars.ymin + j * mvars.gridsp, avars.zmin + k * mvars.gridsp))
                    outfile.write('%.3e ' % (avars.cbin[i][j][k]))
                    for ii in xrange(mvars.nsegments):
                        outfile.write('%.3e ' % avars.sbin[ii][i][j][k])
                        for jj in range(avars.anregions[ii]):
                            outfile.write('%.3e ' % avars.arbin[ii][jj][i][j][k])
                    outfile.write('\n')

        filename = str(outfile).split("'")[1]
        pgui("Wrote occupancy data to : %s\n\n" % filename)

        return


    def initialization(self):
        '''
        method to prepare for density
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        if (mvars.runname[-1] == '/'):
            log.debug('runname(1) = %s' % (mvars.runname))
            avars.densitypath = mvars.runname + 'density_plot/'
            log.debug('densitypath = %s' % (avars.densitypath))
        else:
            log.debug('runname(2) = %s' % (mvars.runname))
            avars.densitypath = mvars.runname + '/density_plot/'
            log.debug('densitypath = %s' % (avars.densitypath))

        direxist = os.path.exists(avars.densitypath)
        if(direxist == 0):
            os.system('mkdir -p ' + avars.densitypath)

        if(mvars.equalweights == 0):
            avars.wst = '_unequalweights'
        else:
            avars.wst = '_equalweights'

        avars.ang2au = 1.0 / 0.5291772108

        avars.m1 = sasmol.SasMol(0)
        avars.m1.read_pdb(mvars.pdbfile)

        pgui('calculating min and max over all frames')
        if(mvars.dcdfile[-3:] == 'dcd'):
            minmax_array = avars.m1.calc_minmax_all_steps(mvars.dcdfile)
            avars.ldcdfile = avars.m1.open_dcd_read(mvars.dcdfile)
            avars.nf = avars.ldcdfile[2]
            avars.intype = 'dcd'
        elif(mvars.dcdfile[-3:] == 'pdb'):
            avars.m1.read_pdb(mvars.dcdfile)
            avars.nf = avars.m1.number_of_frames()
            minmax_array = avars.m1.calc_minmax_all_steps(mvars.dcdfile, pdb='pdb')
            avars.intype = 'pdb'

        log.debug('number of frames = %i' % (avars.nf))
        log.debug('intype = %s' % (avars.intype))

        avars.outfile = (open(avars.densitypath + mvars.ofile + '_' +
                        str(mvars.gridsp) + avars.wst + '_complete.cube', 'w'))

        for i in xrange(mvars.nsegments):
            log.debug('\nsegment = %s :\n\tnum regions = %i' % (str(i+1),avars.anregions[i]))
            log.debug('\tlow = %s' % (str(avars.alow[i])))
            log.debug('\thigh = %s\n' % (str(avars.ahigh[i])))
    
        avars.soutfile = []
#        sstring = []       not used
        for i in xrange(mvars.nsegments):
            filest = avars.densitypath + mvars.ofile + '_' + \
                str(mvars.gridsp) + avars.wst + '_segment_' + str(i + 1) + '_complete.cube'
#            sstring.append(filest)
            avars.soutfile.append(open(filest, 'w'))

        avars.aroutfile = []
#        arstring = []      not used
        for i in xrange(mvars.nsegments):
            routfile = []
            rstring = []
            for j in range(avars.anregions[i]):
                filest = avars.densitypath + mvars.ofile + '_' + \
                    str(mvars.gridsp) + avars.wst + '_segment_' + str(i + 1) + \
                    '_region_' + str(j + 1) + '.cube'
                rstring.append(filest)
                routfile.append(open(filest, 'w'))
            avars.aroutfile.append(routfile)
#            arstring.append(rstring)

        if(mvars.equalweights == 1):
            avars.weights = numpy.ones(avars.nf, numpy.float32)
#            avars.sumweights = float(avars.nf)     not used
        else:
            winfile = open(mvars.weightsfile, 'r').readlines()
            pweights = []
            for i in range(len(winfile)):
                # if(i>1):
                try:
                    lin = string.split(winfile[i])
                    lw = locale.atof(lin[2])
                    pweights.append(lw)
                except:
                    message = 'could not read line = ' + str(i + 1)
                    message += ': stopping here'
                    pgui(message)

            avars.weights = numpy.array(pweights)
#            avars.sumweights = numpy.sum(weights)
            log.debug('weights:  %s' % (str(avars.weights)))


        total_min_array = minmax_array[0]
        total_max_array = minmax_array[1]

        pgui('total_min_array = %s' % (str(total_min_array)))
        pgui('total_max_array = %s' % (str(total_max_array)))

        actual_xlength = total_max_array[0] - total_min_array[0]
        actual_ylength = total_max_array[1] - total_min_array[1]
        actual_zlength = total_max_array[2] - total_min_array[2]

        pgui('actual xlength = %f' % (actual_xlength))
        pgui('actual ylength = %f' % (actual_ylength))
        pgui('actual zlength = %f' % (actual_zlength))

        avars.xmin = numpy.floor(-mvars.xlength / 2.0)
        xmax = numpy.ceil(mvars.xlength / 2.0)
        avars.ymin = numpy.floor(-mvars.ylength / 2.0)
        ymax = numpy.ceil(mvars.ylength / 2.0)
        avars.zmin = numpy.floor(-mvars.zlength / 2.0)
        zmax = numpy.ceil(mvars.zlength / 2.0)

        if(avars.xmin > total_min_array[0] - mvars.gridsp):
            avars.xmin = total_min_array[0] - mvars.gridsp  # make it slightly smaller
            pgui('reset xmin to fit molecule: %s' % (str(avars.xmin)))
        if(avars.ymin > total_min_array[1] - mvars.gridsp):
            avars.ymin = total_min_array[1] - mvars.gridsp  # make it slightly smaller
            pgui('reset ymin to fit molecule: %s' % (str(avars.ymin)))
        if(avars.zmin > total_min_array[2] - mvars.gridsp):
            avars.zmin = total_min_array[2] - mvars.gridsp  # make it slightly smaller
            pgui('reset zmin to fit molecule: %s'  % (str(avars.zmin)))
        if(xmax < total_max_array[0] + mvars.gridsp):
            xmax = total_max_array[0] + mvars.gridsp  # make it slightly larger
            pgui('reset xmax to fit molecule: %s' % (str(xmax)))
        if(ymax < total_max_array[1] + mvars.gridsp):
            ymax = total_max_array[1] + mvars.gridsp  # make it slightly larger
            pgui('reset ymax to fit molecule: %s' % (str(ymax)))
        if(zmax < total_max_array[2] + mvars.gridsp):
            zmax = total_max_array[2] + mvars.gridsp  # make it slightly larger
            pgui('reset zmax to fit molecule: %s' % (str(zmax)))

        avars.nxgp = 1 * int((xmax - avars.xmin) / mvars.gridsp)
        avars.nygp = 1 * int((ymax - avars.ymin) / mvars.gridsp)
        avars.nzgp = 1 * int((zmax - avars.zmin) / mvars.gridsp)

        basis_filter_1 = ''
        for i in xrange(mvars.nsegments):
            basis_filter_1 += 'segname[i] == "' + mvars.allsname[i] + \
                '" and name[i] == "' + mvars.allsbasis[i] + '"'
            if i < (mvars.nsegments - 1):
                basis_filter_1 += ' or '

        error, avars.full_mask = avars.m1.get_subset_mask(basis_filter_1)
        log.debug('full_mask = %s %s %s' % (str(avars.full_mask[4]), str(avars.full_mask[21]), str(avars.full_mask[22])))
        avars.number_atoms_full = numpy.sum(avars.full_mask)
        log.debug('number of atoms full = %i' % (avars.number_atoms_full))

        avars.segments_mask = []
        avars.number_atoms_segments = []

        for i in xrange(mvars.nsegments):
            basis_filter_1s = 'segname[i] == "' + mvars.allsname[i] + \
                '" and name[i] == "' + mvars.allsbasis[i] + '"'
            error, segment_mask = avars.m1.get_subset_mask(basis_filter_1s)
            avars.segments_mask.append(segment_mask)
            avars.number_atoms_segments.append(numpy.sum(segment_mask))
            log.debug('number of atoms segment %s = %s' % (str(i), str(numpy.sum(segment_mask))))

        avars.aregions_mask = []
        avars.anumber_atoms_regions = []

        for i in xrange(mvars.nsegments):
            regions_mask = []
            number_atoms_regions = []
            for j in range(avars.anregions[i]):
                basis_filter_2 = 'name[i] == "' + mvars.allsbasis[i] + '" and segname[i] == "' + mvars.allsname[
                    i] + '" and (resid[i] >= ' + str(avars.alow[i][j]) + ' and resid[i] <= ' + str(avars.ahigh[i][j]) + ')'
                error, region_mask = avars.m1.get_subset_mask(basis_filter_2)
                regions_mask.append(region_mask)
                number_atoms_regions.append(numpy.sum(region_mask))
                log.debug('number of atoms segment %s region %s = %s' % (str(i), str(j), str(numpy.sum(region_mask))))

            avars.aregions_mask.append(regions_mask)
            avars.anumber_atoms_regions.append(number_atoms_regions)

        log.debug(vars(mvars))
        log.debug(vars(avars))

        avars.cbin = numpy.zeros((avars.nxgp, avars.nygp, avars.nzgp), numpy.float32)

        avars.sbin = []
        for i in xrange(mvars.nsegments):
            avars.sbin.append(numpy.zeros((avars.nxgp, avars.nygp, avars.nzgp), numpy.float32))

        avars.arbin = []
        for i in xrange(mvars.nsegments):
            rbin = []
            for j in range(avars.anregions[i]):
                rbin.append(numpy.zeros((avars.nxgp, avars.nygp, avars.nzgp), numpy.float32))
            avars.arbin.append(rbin)




    def density(self,variables, segvariables, txtOutput):
        '''
        DENSITY_PLOT is the function to read in variables from GUI input and compare
        generate three-dimensional volumetric data files using the GAUSSIAN file
        format. 

        INPUT:  variable descriptions:

                runname:        project name
                dcdfile:        input filename
                pdbfile:        reference pdb name
                xlength:        x boxlength
                gridsp:         grid spacing (angstroms)
                ylength:        y boxlength
                save_occupancy: whether save the unweighted raw cube data or not
                zlength:        z boxlength
                nsegments:	number of segments
                nregions:       number of regions per segment
                lowregions:     low region array per segment
                highregions:    high region array per segment
                equalweights:   use equalweights (1=yes) or weights from file (0=no)
                weightsfile:    filename containing weights per structure (see BEST)

        OUTPUT:
                ofile:              	output filename prefix

                files stored in ~/runname/filter directory:

                *_complete.cube:	Gaussian volumetric cube file of all basis atoms
                *_region_X.cube:	Gaussian volumetric cube file of basis atoms in region X


        The output density will be normalized as follows against the maximum density value from the composite map including all atoms:
        rho(norm)[i][j][k] = rho[i][j][k]*100.0/max(rho)
        where i=1,2,...,Nsegments, j=1,2,...,Nregions, and k=1,2,...,Ngridpoints
        '''
#mvars:     runname, dcdfile, pdbfile, ofile, xlength, gridsp, ylength, save_occupancy, zlength, nsegments, equalweights, weightsfile
#           anregions, alow, ahigh, allsbasis, allsname

        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in density')

        mvars = self.mvars
        avars = self.avars



    # ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        for i in range(avars.nf):

            if(avars.intype == 'dcd'):
                avars.m1.read_dcd_step(avars.ldcdfile, i)
                error, cz = avars.m1.get_coor_using_mask(0, avars.full_mask)
            elif(avars.intype == 'pdb'):
                error, cz = avars.m1.get_coor_using_mask(i, avars.full_mask)

            weight = avars.weights[i]

            try:
 
                cube.cube(cz[0], avars.cbin, weight, avars.xmin, avars.ymin,
                        avars.zmin, mvars.gridsp, avars.nxgp, avars.nygp, avars.nzgp)

                for j in range(mvars.nsegments):
                    if(avars.intype == 'dcd'):
                        error, lcz = avars.m1.get_coor_using_mask(0, avars.segments_mask[j])
                    elif(avars.intype == 'pdb'):
                        error, lcz = avars.m1.get_coor_using_mask(i, avars.segments_mask[j])
                    cube.cube(lcz[0], avars.sbin[j], weight, avars.xmin, avars.ymin,
                            avars.zmin, mvars.gridsp, avars.nxgp, avars.nygp, avars.nzgp)

                for j in range(mvars.nsegments):
                    for k in xrange(avars.anregions[j]):
                        if(avars.intype == 'dcd'):
                            error, lcz = avars.m1.get_coor_using_mask(
                                0, avars.aregions_mask[j][k])
                        elif(avars.intype == 'pdb'):
                            error, lcz = avars.m1.get_coor_using_mask(
                                i, avars.aregions_mask[j][k])
                        cube.cube(lcz[0], avars.arbin[j][k], weight, avars.xmin,
                                avars.ymin, avars.zmin, mvars.gridsp, avars.nxgp, avars.nygp, avars.nzgp)
            except:
                message = 'try increasing boxsize via xlength, ylength, and zlength'
                message += ' : stopping here'
                pgui(message)
                return

            if(((i + 1) % (float(avars.nf) / 10.0) == 0 or (avars.nf < 10))):
#                pstring = 'i = ', str(i + 1), 'of ', str(avars.nf), ' steps: ', str((float(i + 1) / float(avars.nf)) * 100.0), ' percent done'
#                pgui(pstring)
                fraction_done = (float(i + 1) / float(avars.nf))                
                progress_string = 'COMPLETED ' + \
                    str(i + 1) + ' of ' + str(avars.nf) + ' : ' + \
                    str(fraction_done * 100.0) + ' % done'
                pgui('%s\n' % progress_string)
                report_string = 'STATUS\t' + str(fraction_done)
                pgui(report_string)

        if(mvars.save_occupancy == "Y"):
            outfile_ro = (open(avars.densitypath + mvars.ofile + '_' +
                            str(mvars.gridsp) + avars.wst + '_occupancy.txt', 'w'))
 #           write_occupancy(m1, nxgp, nygp, nzgp, nsegments, anregions, cbin, sbin, arbin, number_atoms_full, number_atoms_segments,
  #                          anumber_atoms_regions, allsbasis, allsname, alow, ahigh, xmin, ymin, zmin, ang2au, gridsp, outfile_ro, txtOutput)
            self.write_occupancy(outfile_ro)
            outfile_ro.close()

        pgui('\nrenormalizing bins')

        # now renormalize the bins so that density 0 to 100%

        maxcbin = renorm.renorm(avars.cbin, avars.nxgp, avars.nygp, avars.nzgp)
        maxsbin = []
        for i in xrange(mvars.nsegments):
            maxsbin.append(renorm.renorm(avars.sbin[i], avars.nxgp, avars.nygp, avars.nzgp))
        amaxreg = []
        for i in xrange(mvars.nsegments):
            maxreg = []
            for j in range(avars.anregions[i]):
                maxreg.append(renorm.renorm(avars.arbin[i][j], avars.nxgp, avars.nygp, avars.nzgp))
            amaxreg.append(maxreg)

        if(maxcbin == 0.0):
            pgui('maxcbin = %f' % (maxcbin))
            message = 'divide by zero error in renormalization step (1)'
            message += ' maxbin = ' + str(maxcbin)
            message += ' :  stopping here'
            pgui(message)
            return

        else:
            avars.normcbin = 100.0 / maxcbin


#NOTE:  normsbin and anormreg are not used in original calls to write_cube below (now commented out).  Only normcbin is used.  Is this correct?
#If so, then we don't need normsbin or anormreg

        normsbin = []
        for i in xrange(mvars.nsegments):
            if(maxsbin[i] == 0.0):
                pgui("maxsbin[%i] = %f\n" % (i, maxsbin[i]))
                message = 'divide by zero error in renormalization step (2)'
                message += 'i = ' + str(i) + ' maxsbin = ' + str(maxsbin[i])
                message += ' :  stopping here'
                pgui(message)
                return
            else:
                normsbin.append(100.0 / maxsbin[i])

        anormreg = []
        for i in xrange(mvars.nsegments):
            normreg = []
            for j in range(avars.anregions[i]):
                if(amaxreg[i][j] == 0.0):
                    pgui("amaxreg[%i][%i] = %f\n" % (i, j, amaxreg[i][j]))
                    message = 'divide by zero error in renormalization step (2)'
                    message += 'i = ' + str(i) + ' j = ' + \
                        str(j) + ' amaxreg = ' + str(amaxreg[i][j])
                    message += ' :  stopping here'
                    pgui(message)
                    return

                else:
                    normreg.append(100.0 / amaxreg[i][j])

            anormreg.append(normreg)

#        z = 0.0  not used

        fraction_done = 0.0
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        pgui('writing cube files to disk')

        fltr = ''
        for i in xrange(mvars.nsegments):
            fltr += 'segname[natom] == "' + mvars.allsname[i] + \
                '" and name[natom] == "' + mvars.allsbasis[i] + '"'
            if i < (mvars.nsegments - 1):
                fltr += ' or '
#        write_cube(m1, fltr, nxgp, nygp, nzgp, cbin, normcbin, number_atoms_full,
#                xmin, ymin, zmin, ang2au, gridsp, outfile, txtOutput)
        log.debug('fltr in first call to write_cube: %s' % (fltr))
        self.write_cube(fltr, avars.cbin, avars.number_atoms_full, avars.outfile)
        avars.outfile.close()

        for i in xrange(mvars.nsegments):

            fraction_done = 0.0
            report_string = 'STATUS\t' + str(fraction_done)
            pgui(report_string)

            fltr = 'name[natom]=="' + mvars.allsbasis[i] + \
                '" and segname[natom]=="' + mvars.allsname[i] + '"'
#            write_cube(m1, fltr, nxgp, nygp, nzgp, sbin[i], normcbin, number_atoms_segments[
#                    i], xmin, ymin, zmin, ang2au, gridsp, soutfile[i], txtOutput)
            log.debug('fltr in 2nd call to write_cube:  %s' % (fltr))
            self.write_cube(fltr, avars.sbin[i], avars.number_atoms_segments[i], avars.soutfile[i])
            avars.soutfile[i].close()

        for i in xrange(mvars.nsegments):

            for j in range(avars.anregions[i]):

                fraction_done = 0.0
                report_string = 'STATUS\t' + str(fraction_done)
                pgui(report_string)

                fltr = 'name[natom]=="' + mvars.allsbasis[i] + '" and segname[natom]=="' + mvars.allsname[i] + \
                    '" and (resid[natom]>=' + str(avars.alow[i][j]) + \
                    ' and resid[natom]<=' + str(avars.ahigh[i][j]) + ')'
#                write_cube(m1, fltr, nxgp, nygp, nzgp, arbin[i][j], normcbin, anumber_atoms_regions[
#                        i][j], xmin, ymin, zmin, ang2au, gridsp, aroutfile[i][j], txtOutput)
                log.debug('flter in 3rd call to write_cube: %s' % (fltr))
                self.write_cube(fltr, avars.arbin[i][j], avars.anumber_atoms_regions[i][j], avars.aroutfile[i][j])
                avars.aroutfile[i][j].close()

        if(avars.intype == 'dcd'):
            avars.m1.close_dcd_read(avars.ldcdfile[0])

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        st = ''.join(['=' for x in xrange(60)])
        pgui("\n%s \n\n" % (st))        
        pgui('DENSITY_PLOT IS DONE')

        time.sleep(1.0)

        return

