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
import locale
import string
import sassie.interface.input_filter as input_filter
import sasmol.sasmol as sasmol


#used for SASSIE_gui?  Not tested.
def check_psegvariables_exist(psegvariables, glbls):
    error = []
    if psegvariables not in glbls:
        error.append(
            'Please enter the segment variables by clicking on the button "click to enter segment information"!')
    return error


#check segment variables
def check_segvariables(pdbfile, segvariables, nsegments):
    error = []

    otmp = sasmol.SasMol(0)
    otmp.read_pdb(pdbfile, fastread=True)
    totnsegments = len(set(otmp.segname()))

#    print 'nsegments, len(segvariables), totnsegments: ', nsegments, len(segvariables), totnsegments

    if len(segvariables) != nsegments:
        error.append(
            'The number of segments in segvariables is not equal to the value of nsegments!')
        return error

    if len(segvariables) != totnsegments:
        error.append(
            'The number of segments in segvariables is not equal to that from the pdb file!')
        return error

    allsnumranges = []
    allsrlow = []
    allsrhigh = []
    allsegbasis = []
    allsegname = []
    for i in range(len(segvariables)):
        allsnumranges.append(segvariables[i][0])
        allsrlow.append(segvariables[i][1])
        allsrhigh.append(segvariables[i][2])
        allsegbasis.append(segvariables[i][3])
        allsegname.append(segvariables[i][4])

    anumranges = []
    arlow = []
    arhigh = []
    for i in range(len(allsnumranges)):
        try:
            nr = locale.atoi(allsnumranges[i])
        except:
            error.append('The number of ranges "' + allsnumranges[i] + '" for segment number ' + str(
                i) + ' in the segment input fields should be an integer type!')
            return error
        if nr < 1:
            error.append('The number of ranges "' + allsnumranges[i] + '" for segment number ' + str(
                i) + ' in the segment input fields should be equal/greater than 1!')
            return error
        anumranges.append(nr)

    for i in range(len(allsrlow)):
        linrlow = string.split(allsrlow[i], ',')
        linrhigh = string.split(allsrhigh[i], ',')
        rlow = []
        rhigh = []
        for k in range(len(linrlow)):
            try:
                trlow = locale.atoi(linrlow[k])
            except:
                error.append('The low resid "' + allsrlow[i] + '" for segment number ' + str(
                    i) + ' in the segment input fields should be an integer array!')
                return error
            rlow.append(trlow)
        for k in range(len(linrhigh)):
            try:
                trhigh = locale.atoi(linrhigh[k])
            except:
                error.append('The high resid "' + allsrhigh[i] + '" for segment number ' + str(
                    i) + ' in the segment input fields should be an integer array!')
                return error
            rhigh.append(trhigh)
        arlow.append(rlow)
        arhigh.append(rhigh)

    for i in range(len(segvariables)):
        numranges = anumranges[i]
        segbasis = allsegbasis[i]
        segname = allsegname[i]

        if segname.strip() not in otmp.segname():
            error.append('The segment name "' + segname +
                         '" is not found in the pdb file!')
            return error
        elif segbasis.strip() not in otmp.name():
            error.append('The segment basis "' + segbasis +
                         '" is not found in the pdb file!')
            return error

#This test seems to be covered when the number of ranges is tested for each segment above.  
#Not tested but kept just in case.
        if numranges < 1:
            error.append('number of ranges ' + str(numranges) +
                         ' should be equal or greater than 1!')
            return error

        if numranges != len(arlow[i]) or numranges != len(arhigh[i]):
            error.append(
                'the number of low/high residue input does not match the number of ranges!')
            error.append('arlow = ' + str(arlow) +
                         ' : arhigh = ' + str(arhigh))
            error.append('len(arlow) = ' + str(len(arlow)) +
                         ' : len(arhigh) = ' + str(len(arhigh)))
            error.append('numranges = ' + str(numranges))
            return error

        locvariables = ['resid']
        value, result = input_filter.get_pdb_complex_stats(
            pdbfile, segname, locvariables)
        if value == 1:
            resid = map(int, result[0])
#value==0 doesn't seem to occur if pdb file passes tests for valid pdb file above
#Not tested, but kept just in case
        elif value == 0:
            error.append(
                'cannot get the pdb statistics for segment ' + segname)
            return error
        number_aa = resid[-1] - resid[0] + 1

        for j in range(numranges):
            lowres1 = arlow[i][j]
            highres1 = arhigh[i][j]
            if(lowres1 not in resid):
                error.append('input pdb file, ' + str(os.path.basename(pdbfile)) + ' does not have low residue amino acid, ' +
                             str(lowres1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
                return error
            elif(highres1 not in resid):
                error.append('input pdb file, ' + str(os.path.basename(pdbfile)) + ' does not have high residue amino acid, ' +
                             str(highres1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
                return error
            elif(lowres1 >= highres1):
                error.append('the low-resid: ' + str(lowres1) +
                             ', is greater than/equal to the high-resid: ' + str(highres1) + ' in segname: ' + segname)
                return error
    return error


def check_density_plot(variables, psegvariables, **kwargs):

    runname = variables['runname'][0]
    path = variables['path'][0]
    dcdfile = variables['dcdfile'][0]
    pdbfile = variables['pdbfile'][0]
    ofile = variables['ofile'][0]
    nsegments = variables['nsegments'][0]
    xlength = variables['xlength'][0]
    gridsp = variables['gridsp'][0]
    ylength = variables['ylength'][0]
    save_occupancy = variables['save_occupancy'][0]
    zlength = variables['zlength'][0]
    equalweights = variables['equalweights'][0]
    weightsfile = variables['weightsfile'][0]

    error = []
    error = input_filter.check_name(runname)
    if(error != []):
        return error

    if 'no_file_check' not in kwargs:
        ev, rv, wv = input_filter.check_permissions(path)
        if(not ev or not rv or not wv):
            error.append('permission error in input file path ' +
                         path + '  [code = ' + str(ev) + str(rv) + str(wv) + ']')
            if(ev == False):
                error.append('path does not exist')
            elif(rv == False):
                error.append('read permission not allowed')
            elif(wv == False):
                error.append('write permission not allowed')
            return error

#check input pdb file
    ev, value = input_filter.check_pdb_dcd(pdbfile, 'pdb')

    if(ev == 0):
        error.append('input pdb file, ' +
                     os.path.basename(pdbfile) + ', does not exist')
        return error
    if(value == 0):
        error.append('input pdb file, ' +
                     os.path.basename(pdbfile) + ', is not a valid pdb file')
        return error

#check input trajectory (pdb or dcd) file
    newfile = os.path.basename(dcdfile)
    try:
        if (newfile[0] == '.' or newfile[-4:] not in ['.pdb', '.dcd']):
            error.append(
                'input trajectory filename cannot start with "." and must end with .pdb or .dcd : ' + newfile)
            return error
#could not find an example that was not already handled by input_filter
#Not tested, but kept just in case
    except:
        error.append(
            'error checking file name ', + dcdfile)
        return error

#if file exists, need to find number of frames for comparison to weight file below
    error = input_filter.check_file_exists(dcdfile)
    if(len(error) != 0):
        return error

    ev, value = input_filter.check_pdb_dcd(dcdfile, 'dcd')
    if(value == 0):
        ev, value = input_filter.check_pdb_dcd(dcdfile, 'pdb')
        if(value == 0):
            error.append('input file, ' + os.path.basename(dcdfile) +
                     ', is not a valid dcd or pdb file')
            return error
        infile_type = 'pdb'
    else:
        infile_type = 'dcd'

    if(infile_type == 'dcd'):        
        value = input_filter.certify_pdb_dcd(pdbfile, dcdfile)
#        print '1: value = ',value
        if(value == 0):
            error.append('input pdb ' + os.path.basename(pdbfile) + ' and dcd file ' +
                         os.path.basename(dcdfile) + ' are not compatible')
            return error
    elif(infile_type == 'pdb'):
        ev, value = input_filter.certify_pdb_pdb(pdbfile, dcdfile)
#        print '2: value = ',value
        if(value == 0):
            error.append('input pdb ' + os.path.basename(pdbfile) + ' and dcd file ' +
                        os.path.basename(dcdfile) + ' are not compatible')
            return error


# exception below not tested since other tests failed before reaching this point; if this point is reached, file is read successfully
# kept exception just in case
# note that number_of_frames is used later to compare to number of weights in weight file
    try:
        m1 = sasmol.SasMol(0)
        if(infile_type == 'dcd'):
            dcdinputfile = m1.open_dcd_read(dcdfile)
            number_of_frames = dcdinputfile[2]
        elif(infile_type == 'pdb'):
            m1.read_pdb(dcdfile)
            number_of_frames = m1.number_of_frames()
#        print 'infile_type, number_of_frames: ', infile_type, number_of_frames    
    except:
        error.append(
            'could not read frames in dcd file ', + dcdfile)
        return error


#   check xlength, ylength, zlength, gridsp, equalweights, save_occupancy
    if(xlength <= 0.0 or ylength <= 0.0 or zlength <= 0.0):
        error.append('boxlengths need to be > 0, lengths = ' +
                     str(xlength) + ',' + str(ylength) + ',' + str(zlength))
        return error
    elif(gridsp <= 0):
        error.append(
            'grid spacing needs to be > 0, grid spacing = ' + str(gridsp))
        return error
    elif(equalweights != 0 and equalweights != 1):
        error.append(
            'equalweights == 0 for "no" and 1 for "yes", equalweights = ' + str(equalweights))
        return error
    elif(save_occupancy != "Y" and save_occupancy != "N"):
        error.append(
            'save occupancy data needs to be either "Y" or "N": you entered: ' + save_occupancy)
        return error


#check weight file
    if(equalweights == 0):
        error = input_filter.check_file_exists(weightsfile)
        if(len(error) > 0):
            error.append('weight file not readable or does not exist')
            return error

        # check if the weights file has at least one non-zero value

## 	structure, X2, weight
# 1       11.462422       0.000000

        try:
            infile = open(weightsfile, 'r').readlines()
            if(len(infile) < 1):
                error.append("no lines in your weight file : " +
                             os.path.basename(weightsfile))
                return error             

        # check if the length of the weight is equal to the number of frames in the trajectory file
#            print 'len(infile): ', len(infile)
#            print 'number_of_frames: ', number_of_frames
            if(len(infile) != number_of_frames):
                error.append('number of lines ' + str(len(infile)) + ' in weightsfile, ' + os.path.basename(weightsfile) + 
                ', does not match the number of frames in input file, ' + os.path.basename(dcdfile) + '; number of frames = ' + str(number_of_frames))

            frames = []
            weights = []
            for lin in infile:
                words = lin.split()
                if not len(words):
                    continue
                if words[0][0] == '#':
                    continue
                else:
                    if locale.atof(words[2]) not in [0.0, 1.0]:
                        error.append(
                            'weight file column 3 can only have 0 or 1 : ' + words[2] + ' was found')
                        return error
                    weights.append(locale.atof(words[2]))
                    frames.append(locale.atoi(words[0]))
            frames.sort()

            local_count = 0

            for this_weight in weights:
                if this_weight == 0.0:
                    local_count += 1

            if local_count == len(weights):
                error.append(
                    "all weights in weights file are zero : " + os.path.basename(weightsfile))
                return error

            for i in range(len(frames)):
                if(frames[i] < 1):
                    error.append('weight file column 1 can only have positive integers : "' + str(
                        frames[i]) + '" was found in the weight file, ' + os.path.basename(weightsfile))
                    return error
#NOTE:  this test is different than the one above, i.e., it tests if the number in column 1 is greater than the number of frames in the trajectory file
                if(frames[i] > number_of_frames):
                    error.append('there are ' + str(number_of_frames) + ' frames in input file ' + os.path.basename(dcdfile) +
                ': "'+ str(frames[i]) + '" was found in weight file, ' + os.path.basename(weightsfile))
                    return error
                if(i > 0 and frames[i] == frames[i - 1]):
                    error.append('redundant number "' +
                                 str(frames[i]) + '" found in the weight file, ' + os.path.basename(weightsfile))
                    return error
        except:
            error.append(
                "unable to open and read your weight file : " + os.path.basename(weightsfile))
            return error

    lerror = check_segvariables(pdbfile, psegvariables, nsegments)
    if len(lerror):
#        print 'segvariables error: ', lerror
        error.append(lerror)
        return error

    return error


