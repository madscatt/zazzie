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
import sassie.sasmol.sasmol as sasmol
import sassie.analyze.cube as cube
import sassie.analyze.renorm as renorm
import sassie.interface.input_filter as input_filter

#       DENSITY_PLOT
#
#       11/20/2005       --      initial coding                 :	jc
#	03/15/2012	 --	 added segment based cube files :	jc
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
	the binning (see: code.c and renorm.c).

'''


def print_failure(message, txtOutput):

    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
    txtOutput.put(message)

    return


def unpack_variables(variables):

    runname = variables['runname'][0]
    dcdfile = variables['dcdfile'][0]
    pdbfile = variables['pdbfile'][0]
    ofile = variables['ofile'][0]
    xlength = variables['xlength'][0]
    gridsp = variables['gridsp'][0]
    ylength = variables['ylength'][0]
    save_occupancy = variables['save_occupancy'][0]
    zlength = variables['zlength'][0]
    nsegments = variables['nsegments'][0]
    equalweights = variables['equalweights'][0]
    weightsfile = variables['weightsfile'][0]

    return runname, dcdfile, pdbfile, ofile, xlength, gridsp, ylength, save_occupancy, zlength, nsegments, equalweights, weightsfile


def process_input_variables(segvariables, nsegments):

    #nregions        = variables['nregions'][0]
    #lowregions      = variables['lowregions'][0]
    #highregions     = variables['highregions'][0]

    allsnumregions = []
    allslow = []
    allshigh = []
    allsbasis = []
    allsname = []

    for i in range(len(segvariables)):
        allsnumregions.append(segvariables[i][0])
        allslow.append(segvariables[i][1])
        allshigh.append(segvariables[i][2])
        allsbasis.append(segvariables[i][3])
        allsname.append(segvariables[i][4])

    anregions = []

    for i in range(len(allsnumregions)):
        nr = locale.atoi(allsnumregions[i])
        anregions.append(nr)

    alow = []
    ahigh = []

    for i in range(len(allslow)):
        linlow = string.split(allslow[i], ',')
        linhigh = string.split(allshigh[i], ',')
        rlow = []
        rhigh = []
        for k in range(len(linlow)):
            tlow = locale.atoi(linlow[k])
            thigh = locale.atoi(linhigh[k])
            rlow.append(tlow)
            rhigh.append(thigh)
        alow.append(rlow)
        ahigh.append(rhigh)

    return anregions, alow, ahigh, allsbasis, allsname


def stopme():

    print '\n\nDEBUGGING CODE: STOPPING HERE'
    print 'DEBUGGING CODE: STOPPING HERE'
    print 'DEBUGGING CODE: STOPPING HERE\n\n'
    sys.exit()


def write_cube_header(number_atoms, nxgp, nygp, nzgp, xmin, ymin, zmin, ang2au, gridsp, outfile, remark=''):
    st1 = remark + 'SASSIE/GAUSSIAN CUBE FILE\n'
    st2 = remark + 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n'
    outfile.write(remark + 'calpha_' + st1 + st2)
    outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                  (number_atoms, xmin * ang2au, ymin * ang2au, zmin * ang2au))
    outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                  (nxgp, gridsp * ang2au, 0.00, 0.00))
    outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                  (nygp, 0.00, gridsp * ang2au, 0.00))
    outfile.write(remark + '%d\t%f\t%f\t%f\n' %
                  (nzgp, 0.00, 0.00, gridsp * ang2au))


def write_cube(m, fltr, nxgp, nygp, nzgp, abin, normbin, number_atoms, xmin, ymin, zmin, ang2au, gridsp, outfile, txtOutput):

    coor = m.coor()
    name = m.name()
    resid = m.resid()
    segname = m.segname()

    write_cube_header(number_atoms, nxgp, nygp, nzgp, xmin,
                      ymin, zmin, ang2au, gridsp, outfile)

    for natom in range(len(name)):
        if eval(fltr):
            outfile.write('%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                          0, natom, 0] * ang2au, coor[0, natom, 1] * ang2au, coor[0, natom, 2] * ang2au))

    for i in range(nxgp):

        if((((i + 1) % nxgp) / 10.0) == 0 or (nxgp < 10)):
            fraction_done = (float(i + 1) / float(nxgp))
            progress_string = '\nWROTE ' + str(i + 1) + ' of ' + str(
                nxgp) + ' grid points : ' + str(fraction_done * 100.0) + ' % done'
            print('%s\n' % progress_string)
            report_string = 'STATUS\t' + str(fraction_done)
            txtOutput.put(report_string)

        for j in range(nygp):
            for k in range(nzgp):
                outfile.write('%lf ' % (abin[i][j][k] * normbin))
                if((k % 6) == 5):
                    outfile.write('\n')
            outfile.write('\n')

    filename = str(outfile).split("'")[1]
    txtOutput.put("Wrote the density data to : %s\n\n" % filename)


def write_occupancy(m, nxgp, nygp, nzgp, nsegments, anregions, cbin, sbin, arbin, number_atoms_full, number_atoms_segments, number_atoms_regions, allsbasis, allsname, alow, ahigh, xmin, ymin, zmin, ang2au, gridsp, outfile, txtOutput):

    coor = m.coor()
    name = m.name()
    resid = m.resid()
    segname = m.segname()

    fltr = ''
    for i in xrange(nsegments):
        fltr += 'segname[natom] == "' + allsname[i] + \
            '" and name[natom] == "' + allsbasis[i] + '"'
        if i < (nsegments - 1):
            fltr += ' or '
    outfile.write('\n############################\n#Complete\n')
    write_cube_header(number_atoms_full, nxgp, nygp, nzgp,
                      xmin, ymin, zmin, ang2au, gridsp, outfile, remark='#')
    for natom in range(len(name)):
        if eval(fltr):
            outfile.write('#%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                          0, natom, 0] * ang2au, coor[0, natom, 1] * ang2au, coor[0, natom, 2] * ang2au))

    for ii in xrange(nsegments):
        outfile.write('\n############################\n#Seg_%d  \n' % ii)
        fltr = "name[natom]=='" + allsbasis[ii] + \
            "' and segname[natom]=='" + allsname[ii] + "'"
        write_cube_header(number_atoms_segments[
                          ii], nxgp, nygp, nzgp, xmin, ymin, zmin, ang2au, gridsp, outfile, remark='#')
        for natom in range(len(name)):
            if eval(fltr):
                outfile.write('#%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                              0, natom, 0] * ang2au, coor[0, natom, 1] * ang2au, coor[0, natom, 2] * ang2au))
        for jj in range(anregions[ii]):
            outfile.write(
                '\n############################\n#Seg_%d_Reg_%d  \n' % (ii, jj))
            write_cube_header(number_atoms_regions[ii][
                              jj], nxgp, nygp, nzgp, xmin, ymin, zmin, ang2au, gridsp, outfile, remark='#')
            fltr = "name[natom]=='" + allsbasis[ii] + "' and segname[natom]=='" + allsname[ii] + \
                "' and (resid[natom]>=" + str(alow[ii][jj]) + \
                " and resid[natom]<=" + str(ahigh[ii][jj]) + ")"
            for natom in range(len(name)):
                if eval(fltr):
                    outfile.write('#%-7s %f\t%f\t%f\t%f\n' % ('12', 0.00, coor[
                                  0, natom, 0] * ang2au, coor[0, natom, 1] * ang2au, coor[0, natom, 2] * ang2au))

    outfile.write('\n############################\n')
    outfile.write('# X       Y        Z        Complete  ')
    for ii in xrange(nsegments):
        outfile.write('Seg_%d  ' % ii)
        for jj in range(anregions[ii]):
            outfile.write('Seg_%d_Reg_%d ' % (ii, jj))
    outfile.write('\n')

    for i in range(nxgp):

        if((((i + 1) % nxgp) / 10.0) == 0 or (nxgp < 10)):
            fraction_done = (float(i + 1) / float(nxgp))
            progress_string = 'WROTE ' + \
                str(i + 1) + ' of ' + str(nxgp) + ' : ' + \
                str(fraction_done * 100.0) + ' % done'
            print('%s\n' % progress_string)
            report_string = 'STATUS\t' + str(fraction_done)
            txtOutput.put(report_string)

        for j in range(nygp):
            for k in range(nzgp):
                outfile.write('%8.3f %8.3f %8.3f ' % (
                    xmin + i * gridsp, ymin + j * gridsp, zmin + k * gridsp))
                outfile.write('%.3e ' % (cbin[i][j][k]))
                for ii in xrange(nsegments):
                    outfile.write('%.3e ' % sbin[ii][i][j][k])
                    for jj in range(anregions[ii]):
                        outfile.write('%.3e ' % arbin[ii][jj][i][j][k])
                outfile.write('\n')

    filename = str(outfile).split("'")[1]
    txtOutput.put("Wrote occupancy data to : %s\n\n" % filename)


def density(variables, segvariables, txtOutput):
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
    runname, dcdfile, pdbfile, ofile, xlength, gridsp, ylength, save_occupancy, zlength, nsegments, equalweights, weightsfile = unpack_variables(
        variables)

    anregions, alow, ahigh, allsbasis, allsname = process_input_variables(
        segvariables, nsegments)

    if (runname[-1] == '/'):
        print 'runname(1) = ', runname
        densitypath = runname + 'density_plot/'
        print 'densitypath = ', densitypath
    else:
        print 'runname(2) = ', runname
        densitypath = runname + '/density_plot/'
        print 'densitypath = ', densitypath

    direxist = os.path.exists(densitypath)
    if(direxist == 0):
        os.system('mkdir -p ' + densitypath)

    if(equalweights == 0):
        wst = '_unequalweights'
    else:
        wst = '_equalweights'

    ang2au = 1.0 / 0.5291772108

    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdbfile)


    print 'calculating min and max over all frames'
    if(dcdfile[-3:] == 'dcd'):
        minmax_array = m1.calc_minmax_all_steps(dcdfile)
        ldcdfile = m1.open_dcd_read(dcdfile)
        nf = ldcdfile[2]
        intype = 'dcd'
    elif(dcdfile[-3:] == 'pdb'):
        m1.read_pdb(dcdfile)
        nf = m1.number_of_frames()
        minmax_array = m1.calc_minmax_all_steps(dcdfile, pdb='pdb')
        intype = 'pdb'
#
    print 'number of frames = ', nf
    print 'intype = ', intype

    outfile = (open(densitypath + ofile + '_' +
                    str(gridsp) + wst + '_complete.cube', 'w'))

    for i in xrange(nsegments):
        print '\nsegment = ', str(i + 1), ':\n\tnum regions = ', anregions[i]
        print '\tlow = ', alow[i]
        print '\thigh = ', ahigh[i], '\n'

    soutfile = []
    sstring = []
    for i in xrange(nsegments):
        filest = densitypath + ofile + '_' + \
            str(gridsp) + wst + '_segment_' + str(i + 1) + '_complete.cube'
        sstring.append(filest)
        soutfile.append(open(filest, 'w'))

    aroutfile = []
    arstring = []
    for i in xrange(nsegments):
        routfile = []
        rstring = []
        for j in range(anregions[i]):
            filest = densitypath + ofile + '_' + \
                str(gridsp) + wst + '_segment_' + str(i + 1) + \
                '_region_' + str(j + 1) + '.cube'
            rstring.append(filest)
            routfile.append(open(filest, 'w'))
        aroutfile.append(routfile)
        arstring.append(rstring)

    if(equalweights == 1):
        weights = numpy.ones(nf, numpy.float32)
        sumweights = float(nf)
    else:
        winfile = open(weightsfile, 'r').readlines()
        pweights = []
        # if((len(winfile)-2)!=nf):
        #	print 'nf = ',nf
        #	print 'len(winfile) = ',len(winfile)
        #	message = 'length of weightsfile, '+weightsfile+' indicates there are not enought data to match the number of frames in input file'+dcdfile+' nf = '+str(nf)
        #	message+=' :  stopping here'
        #	print_failure(message,txtOutput)
        #	return

        for i in range(len(winfile)):
            # if(i>1):
            try:
                lin = string.split(winfile[i])
                lw = locale.atof(lin[2])
                pweights.append(lw)
            except:
                message = 'could not read line = ' + str(i + 1)
                message += ': stopping here'
                print_failure(message, txtOutput)

#This is now checked in density_plot_filter
#        if(len(pweights) != nf):
#            print 'nf = ', nf
#            print 'len(pweights) = ', len(pweights)
#            message = 'number of frames ' + str(len(pweights)) + ' in weightsfile, ' + weightsfile + \
#                ' does not match the number of frames in input file' + \
#                dcdfile + ' nf = ' + str(nf)
#            message += ' :  stopping here'
#            print_failure(message, txtOutput)
#            return

        weights = numpy.array(pweights)
        sumweights = numpy.sum(weights)
        print 'weights:  ', weights

    # print 'calculating min and max over all frames'

    total_min_array = minmax_array[0]
    total_max_array = minmax_array[1]

    print 'total_min_array = ', total_min_array
    print 'total_max_array = ', total_max_array

    actual_xlength = total_max_array[0] - total_min_array[0]
    actual_ylength = total_max_array[1] - total_min_array[1]
    actual_zlength = total_max_array[2] - total_min_array[2]

    print 'actual xlength = ', actual_xlength
    print 'actual ylength = ', actual_ylength
    print 'actual zlength = ', actual_zlength

    xmin = numpy.floor(-xlength / 2.0)
    xmax = numpy.ceil(xlength / 2.0)
    ymin = numpy.floor(-ylength / 2.0)
    ymax = numpy.ceil(ylength / 2.0)
    zmin = numpy.floor(-zlength / 2.0)
    zmax = numpy.ceil(zlength / 2.0)

    if(xmin > total_min_array[0] - gridsp):
        xmin = total_min_array[0] - gridsp  # make it slightly smaller
        print 'reset xmin to fit molecule: ' + str(xmin)
    if(ymin > total_min_array[1] - gridsp):
        ymin = total_min_array[1] - gridsp  # make it slightly smaller
        print 'reset ymin to fit molecule: ' + str(ymin)
    if(zmin > total_min_array[2] - gridsp):
        zmin = total_min_array[2] - gridsp  # make it slightly smaller
        print 'reset zmin to fit molecule: ' + str(zmin)
    if(xmax < total_max_array[0] + gridsp):
        xmax = total_max_array[0] + gridsp  # make it slightly larger
        print 'reset xmax to fit molecule: ' + str(xmax)
    if(ymax < total_max_array[1] + gridsp):
        ymax = total_max_array[1] + gridsp  # make it slightly larger
        print 'reset ymax to fit molecule: ' + str(ymax)
    if(zmax < total_max_array[2] + gridsp):
        zmax = total_max_array[2] + gridsp  # make it slightly larger
        print 'reset zmax to fit molecule: ' + str(zmax)

    nxgp = 1 * int((xmax - xmin) / gridsp)
    nygp = 1 * int((ymax - ymin) / gridsp)
    nzgp = 1 * int((zmax - zmin) / gridsp)

    name = m1.name()
    resid = m1.resid()
    segname = m1.segname()

    basis_filter_1 = ''
    for i in xrange(nsegments):
        basis_filter_1 += 'segname[i] == "' + allsname[i] + \
            '" and name[i] == "' + allsbasis[i] + '"'
        if i < (nsegments - 1):
            basis_filter_1 += ' or '

    error, full_mask = m1.get_subset_mask(basis_filter_1)
    print 'full_mask = ', full_mask[4], full_mask[21], full_mask[22]
    number_atoms_full = numpy.sum(full_mask)
    print 'number of atoms full = ', number_atoms_full

    segments_mask = []
    number_atoms_segments = []

    for i in xrange(nsegments):
        basis_filter_1s = 'segname[i] == "' + allsname[i] + \
            '" and name[i] == "' + allsbasis[i] + '"'
        error, segment_mask = m1.get_subset_mask(basis_filter_1s)
        segments_mask.append(segment_mask)
        number_atoms_segments.append(numpy.sum(segment_mask))
        print 'number of atoms segment ' + str(i) + ' = ', numpy.sum(segment_mask)

    aregions_mask = []
    anumber_atoms_regions = []

    for i in xrange(nsegments):
        regions_mask = []
        number_atoms_regions = []
        for j in range(anregions[i]):
            basis_filter_2 = 'name[i] == "' + allsbasis[i] + '" and segname[i] == "' + allsname[
                i] + '" and (resid[i] >= ' + str(alow[i][j]) + ' and resid[i] <= ' + str(ahigh[i][j]) + ')'
            error, region_mask = m1.get_subset_mask(basis_filter_2)
            regions_mask.append(region_mask)
            number_atoms_regions.append(numpy.sum(region_mask))
            print 'number of atoms segment ' + str(i) + ' region ' + str(j) + ' = ', numpy.sum(region_mask)

        aregions_mask.append(regions_mask)
        anumber_atoms_regions.append(number_atoms_regions)

    cbin = numpy.zeros((nxgp, nygp, nzgp), numpy.float32)

    sbin = []
    for i in xrange(nsegments):
        sbin.append(numpy.zeros((nxgp, nygp, nzgp), numpy.float32))

    arbin = []
    for i in xrange(nsegments):
        rbin = []
        for j in range(anregions[i]):
            rbin.append(numpy.zeros((nxgp, nygp, nzgp), numpy.float32))
        arbin.append(rbin)

    # ttxt=time.ctime()
    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" % (st))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    for i in range(nf):

        if(intype == 'dcd'):
            m1.read_dcd_step(ldcdfile, i)
            error, cz = m1.get_coor_using_mask(0, full_mask)
        elif(intype == 'pdb'):
            error, cz = m1.get_coor_using_mask(i, full_mask)

        weight = weights[i]
        try:
            cube.cube(cz[0], cbin, weight, xmin, ymin,
                      zmin, gridsp, nxgp, nygp, nzgp)

            for j in range(nsegments):
                if(intype == 'dcd'):
                    error, lcz = m1.get_coor_using_mask(0, segments_mask[j])
                elif(intype == 'pdb'):
                    error, lcz = m1.get_coor_using_mask(i, segments_mask[j])
                cube.cube(lcz[0], sbin[j], weight, xmin, ymin,
                          zmin, gridsp, nxgp, nygp, nzgp)

            for j in range(nsegments):
                for k in xrange(anregions[j]):
                    if(intype == 'dcd'):
                        error, lcz = m1.get_coor_using_mask(
                            0, aregions_mask[j][k])
                    elif(intype == 'pdb'):
                        error, lcz = m1.get_coor_using_mask(
                            i, aregions_mask[j][k])
                    cube.cube(lcz[0], arbin[j][k], weight, xmin,
                              ymin, zmin, gridsp, nxgp, nygp, nzgp)
        except:
            message = 'try increasing boxsize via xlength, ylength, and zlength'
            message += ' : stopping here'
            print_failure(message, txtOutput)
            return

        if(((i + 1) % (float(nf) / 10.0) == 0 or (nf < 10))):
            print 'i = ', i + 1, 'of ', nf, ' steps: ', (float(i + 1) / float(nf)) * 100.0, ' percent done'
            fraction_done = (float(i + 1) / float(nf))
            progress_string = 'COMPLETED ' + \
                str(i + 1) + ' of ' + str(nf) + ' : ' + \
                str(fraction_done * 100.0) + ' % done'
            print('%s\n' % progress_string)
            report_string = 'STATUS\t' + str(fraction_done)
            txtOutput.put(report_string)

    if(save_occupancy == "Y"):
        outfile_ro = (open(densitypath + ofile + '_' +
                           str(gridsp) + wst + '_occupancy.txt', 'w'))
        write_occupancy(m1, nxgp, nygp, nzgp, nsegments, anregions, cbin, sbin, arbin, number_atoms_full, number_atoms_segments,
                        anumber_atoms_regions, allsbasis, allsname, alow, ahigh, xmin, ymin, zmin, ang2au, gridsp, outfile_ro, txtOutput)
        outfile_ro.close()

    '''
    try:
        cbin=cbin/float(number_atoms_full)

        for i in xrange(nsegments):
            sbin[i]=sbin[i]/float(number_atoms_segments[i])

        for i in xrange(nsegments):
            for j in range(anregions[i]):
                arbin[i][j]=arbin[i][j]/float(anumber_atoms_regions[i][j])
    except:
        message='divide by zero error in cbin/rbin renormalization step (0)'
        message+=' number atoms full = '+str(number_atoms_full)
        message+=' number atoms regions '+str(number_atoms_regions[i])
        message+=' this means you found no atoms in your basis'
        message+=' :  stopping here'
        print_failure(message,txtOutput)
        return
    '''

    print '\nrenormalizing bins'

    # now renormalize the bins so that density 0 to 100%

    maxcbin = renorm.renorm(cbin, nxgp, nygp, nzgp)
    maxsbin = []
    for i in xrange(nsegments):
        maxsbin.append(renorm.renorm(sbin[i], nxgp, nygp, nzgp))
    amaxreg = []
    for i in xrange(nsegments):
        maxreg = []
        for j in range(anregions[i]):
            maxreg.append(renorm.renorm(arbin[i][j], nxgp, nygp, nzgp))
        amaxreg.append(maxreg)

    if(maxcbin == 0.0):
        print 'maxcbin = ', maxcbin
        message = 'divide by zero error in renormalization step (1)'
        message += ' maxbin = ' + str(maxcbin)
        message += ' :  stopping here'
        print_failure(message, txtOutput)
        return

    else:
        normcbin = 100.0 / maxcbin

    normsbin = []
    for i in xrange(nsegments):
        if(maxsbin[i] == 0.0):
            print ("maxsbin[%i] = %f\n" % (i, maxsbin[i]))
            message = 'divide by zero error in renormalization step (2)'
            message += 'i = ' + str(i) + ' maxsbin = ' + str(maxsbin[i])
            message += ' :  stopping here'
            print_failure(message, txtOutput)
            return
        else:
            normsbin.append(100.0 / maxsbin[i])

    anormreg = []
    for i in xrange(nsegments):
        normreg = []
        for j in range(anregions[i]):
            if(amaxreg[i][j] == 0.0):
                print ("amaxreg[%i][%i] = %f\n" % (i, j, amaxreg[i][j]))
                message = 'divide by zero error in renormalization step (2)'
                message += 'i = ' + str(i) + ' j = ' + \
                    str(j) + ' amaxreg = ' + str(amaxreg[i][j])
                message += ' :  stopping here'
                print_failure(message, txtOutput)
                return

            else:
                normreg.append(100.0 / amaxreg[i][j])

        anormreg.append(normreg)
    z = 0.0

    fraction_done = 0.0
    report_string = 'STATUS\t' + str(fraction_done)
    txtOutput.put(report_string)

    print 'writing cube files to disk'

    fltr = ''
    for i in xrange(nsegments):
        fltr += 'segname[natom] == "' + allsname[i] + \
            '" and name[natom] == "' + allsbasis[i] + '"'
        if i < (nsegments - 1):
            fltr += ' or '
    write_cube(m1, fltr, nxgp, nygp, nzgp, cbin, normcbin, number_atoms_full,
               xmin, ymin, zmin, ang2au, gridsp, outfile, txtOutput)
    outfile.close()

    for i in xrange(nsegments):

        fraction_done = 0.0
        report_string = 'STATUS\t' + str(fraction_done)
        txtOutput.put(report_string)

        fltr = 'name[natom]=="' + allsbasis[i] + \
            '" and segname[natom]=="' + allsname[i] + '"'
        write_cube(m1, fltr, nxgp, nygp, nzgp, sbin[i], normcbin, number_atoms_segments[
                   i], xmin, ymin, zmin, ang2au, gridsp, soutfile[i], txtOutput)
        # write_cube(m1,fltr,nxgp,nygp,nzgp,sbin[i],normsbin[i],number_atoms_segments[i],xmin,ymin,zmin,ang2au,gridsp,soutfile[i],txtOutput)
        soutfile[i].close()

    for i in xrange(nsegments):

        for j in range(anregions[i]):

            fraction_done = 0.0
            report_string = 'STATUS\t' + str(fraction_done)
            txtOutput.put(report_string)

            fltr = 'name[natom]=="' + allsbasis[i] + '" and segname[natom]=="' + allsname[i] + \
                '" and (resid[natom]>=' + str(alow[i][j]) + \
                ' and resid[natom]<=' + str(ahigh[i][j]) + ')'
            write_cube(m1, fltr, nxgp, nygp, nzgp, arbin[i][j], normcbin, anumber_atoms_regions[
                       i][j], xmin, ymin, zmin, ang2au, gridsp, aroutfile[i][j], txtOutput)
            # write_cube(m1,fltr,nxgp,nygp,nzgp,arbin[i][j],anormreg[i][j],anumber_atoms_regions[i][j],xmin,ymin,zmin,ang2au,gridsp,aroutfile[i][j],txtOutput)
            aroutfile[i][j].close()

    if(intype == 'dcd'):
        m1.close_dcd_read(ldcdfile[0])

    txtOutput.put("\n%s \n\n" % (st))
    time.sleep(2.0)
    print 'DENSITY_PLOT IS DONE'

    return()

if __name__ == "__main__":

    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT
    equalweights = 1
    gridsp = 5.0
    nsegments = 1
    ofile = 'test'
    path = ''
    pdbfile = 'shorter_min3.pdb'
    dcdfile = 'run_0.dcd'
    runname = 'run_0'
    save_occupancy = 'N'
    weightsfile = ' '
    xlength = 100.0
    ylength = 100.0
    zlength = 100.0
    segvariables = [[u'1', '6', '123', u'CA', u'GAG']]
    # END USER EDIT
    # END USER EDIT
    # END USER EDIT

    svariables = {}
    svariables['dcdfile'] = (str(dcdfile), 'string')
    svariables['equalweights'] = (str(equalweights), 'int')
    svariables['gridsp'] = (str(gridsp), 'float')
    svariables['nsegments'] = (str(nsegments), 'int')
    svariables['ofile'] = (str(ofile), 'string')
    svariables['path'] = (str(path), 'string')
    svariables['pdbfile'] = (str(pdbfile), 'string')
    svariables['runname'] = (str(runname), 'string')
    svariables['save_occupancy'] = (str(save_occupancy), 'string')
    svariables['weightsfile'] = (str(weightsfile), 'string')
    svariables['xlength'] = (str(xlength), 'float')
    svariables['ylength'] = (str(ylength), 'float')
    svariables['zlength'] = (str(zlength), 'float')

    error, variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print 'error = ', error
        sys.exit()

    import sassie.interface.density_plot_filter as density_plot_filter

    error = density_plot_filter.check_density_plot(
        variables, segvariables, no_file_check="true")
    if len(error) > 0:
        print 'error = ', error
        sys.exit()

    import multiprocessing
    txtQueue = multiprocessing.JoinableQueue()
    density(variables, segvariables, txtQueue)
