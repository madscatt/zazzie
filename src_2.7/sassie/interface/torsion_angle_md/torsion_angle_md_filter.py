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
import sasmol.sasmol as sasmol
import sassie.interface.input_filter as input_filter


def process_input_variables(psegvariables, segname_molecules):

    all_flexible_segnames = []
    all_snumranges = []
    all_srlow = []
    all_srnum = []
    all_moltype = []
    all_segnames = []
    all_resid = []

    for i in range(len(psegvariables)):
        all_flexible_segnames.append(psegvariables[i][0])
        all_snumranges.append(psegvariables[i][1])
        all_srlow.append(psegvariables[i][2])
        all_srnum.append(psegvariables[i][3])
        all_moltype.append(psegvariables[i][4])

    all_numranges = []

    for i in range(len(all_snumranges)):
        nr = locale.atoi(all_snumranges[i])
        all_numranges.append(nr)

    all_rlow = []
    all_rnum = []

    for i in range(len(all_srlow)):
        linrlow = string.split(all_srlow[i], ',')
        linrnum = string.split(all_srnum[i], ',')
        rlow = []
        rnum = []
        for k in range(len(linrlow)):
            trlow = locale.atoi(linrlow[k])
            trnum = locale.atoi(linrnum[k])
            rlow.append(trlow)
            rnum.append(trnum)
        all_rlow.append(rlow)
        all_rnum.append(rnum)

    print 'all_flexible_segnames = ', all_flexible_segnames
    print 'all_numranges = ', all_numranges
    print 'all_rlow = ', all_rlow
    print 'all_rnum = ', all_rnum
    print 'all_moltype = ', all_moltype

    flexible_segment_variables = {}

    for i in xrange(len(all_flexible_segnames)):

        flexible_segment_variables[all_flexible_segnames[i]] = [
            all_numranges[i], all_rlow[i], all_rnum[i], all_moltype[i]]

    print 'flexible_segment_variables: ', flexible_segment_variables

    return flexible_segment_variables


def check_torsion_angle_md(variables, psegvariables, **kwargs):

    runname = variables['runname'][0]
    infile = variables['infile'][0]
    pdbfile = variables['pdbfile'][0]
    outfile = variables['outfile'][0]
    nsteps = variables['nsteps'][0]
    temperature = variables['temperature'][0]
    rgforce = variables['rgforce'][0]
    rgvalue = variables['rgvalue'][0]

    dna_segnames = variables['dna_segnames'][0]
    #reslow          = variables['reslow'][0]
    #numcont         = variables['numcont'][0]
    topfile = variables['topfile'][0]
    parmfile = variables['parmfile'][0]
    keepout = variables['keepout'][0]
    dcdfreq = variables['dcdfreq'][0]
    charmmexe = variables['charmmexe'][0]
    number_flexible_segments = variables['number_flexible_segments'][0]
    pretamd_min_steps = variables['pretamd_min_steps'][0]
    poll_frequency = variables['poll_frequency'][0]
    path = variables['path'][0]

    error = []

    error = input_filter.check_name(runname)
    if(error != []):
        return error

    pdbfile = os.path.join(path, pdbfile)
    ev, value = input_filter.check_pdb_dcd(pdbfile, 'pdb')
    if(ev == 0):
        error.append('input pdb file, ' + pdbfile + ', does not exist')
        return error
    elif(value == 0):
        error.append('input pdb file, ' + pdbfile +
                     ', is not a valid pdb file')
        return error

    infile = os.path.join(path, infile)
    error = input_filter.check_file_exists(infile)
    if(error != []):
        return error

    ev, value = input_filter.check_pdb_dcd(infile, 'dcd')
# ev == 0 not tested since non-existent file will trigger
# check_file_exists test above
    if(ev == 0):
        error.append('check input trajectory filename : ' + infile)
        return error
    elif(value == 0):
        ev, value = input_filter.check_pdb_dcd(infile, 'pdb')
        if(value == 0):
            error.append('input trajectory file, ' +
                         infile + ', is not a valid pdb file')
            return error
        infile_type = 'pdb'
    else:
        infile_type = 'dcd'

    if(infile_type == 'dcd'):
        value = input_filter.certify_pdb_dcd(pdbfile, infile)
        if(value == 0):
            error.append('input pdb file ' + pdbfile +
                         ' and dcd file ' + infile + ', are not compatible')
            return error
    elif(infile_type == 'pdb'):
        ev, value = input_filter.certify_pdb_pdb(pdbfile, infile)
        if(value == 0):
            error.append('input pdb file ' + pdbfile +
                         ' and pdb file ' + infile + ', are not compatible')
            return error

    # try:
    if True:
        m1 = sasmol.SasMol(0)
        if(infile_type == 'dcd'):
            dcdinputfile = m1.open_dcd_read(infile)
            number_of_frames = dcdinputfile[2]
        elif(infile_type == 'pdb'):
            m1.read_pdb(infile)
            number_of_frames = m1.number_of_frames()

    # except:
# not tested; file that can't be opened or one with no frames would not be valid or compatible with reference pdb and would trigger errors above
# kept just in case
    else:
        error.append('could not open trajectory file ' +
                     infile + ' to check number of frames')
        return error
    if(number_of_frames < 1):
        error.append('trajectory file has no frames : ' + infile)
        return error

    if(nsteps < 1):
        error.append('nsteps = ' + str(nsteps) + '?')
        return error
    elif(temperature < 0):
        error.append(
            'use a positive temperature, temperature = ' + str(temperature))
        return error
    elif(rgvalue < 0):
        error.append('use a positive Rg value, rg value = ' + str(rgvalue))
        return error
    elif(rgforce < 0):
        error.append(
            'use a positive Rg constraint force, rg force = ' + str(rgforce))
        return error
    elif(rgforce > 0 and rgvalue == 0):
        error.append('Rg value must be > 0 when Rg constraint force is > 0, rg value = ' +
                     str(rgvalue) + ', rg force = ' + str(rgforce))
        return error
    elif(keepout not in [0, 1]):
        error.append(
            'option to keep output files needs to be 0 or 1 = ' + str(keepout))
        return error
# not tested since poll frequency isn't input by the user
    elif(poll_frequency < 0):
        error.append(
            'poll frequency needs to be greater than 0 = ' + str(poll_frequency))
        return error
    elif(number_flexible_segments < 1):
        error.append('number of flexible segments needs to be > 0 : ' +
                     str(number_flexible_segments))
        return error
    if(len(dna_segnames) > 0):
        #dna_segnames = dna_segnames.split(',')
        dna_segnames = [segname.strip() for segname in dna_segnames.split(',')]
        print 'dna_segnames = ', dna_segnames
        m = sasmol.SasMol(0)
        m.read_pdb(pdbfile)
        segnames = m.segnames()
        for this_dna in dna_segnames:
            if this_dna not in segnames:
                error.append('dna segname : ' + str(this_dna) +
                             ' is not in your PDB file')
                return error
    elif(not os.path.exists(topfile)):
        error.append('topology file does not exist : ' + str(topfile))
        return error
    elif(not os.path.exists(parmfile)):
        error.append('parameter file does not exist : ' + str(parmfile))
        return error
    elif(not os.path.exists(charmmexe)):
        error.append(
            'CHARMM executable file does not exist : ' + str(charmmexe))
        return error
    elif(dcdfreq > nsteps):
        error.append(
            'dcdfreq must less than or equal to the number of TAMD steps: ' + str(dcdfreq))
        return error
    elif((nsteps % dcdfreq) != 0):
        error.append(
            'dcdfreq must be a divisor of number of TAMD steps: ' + str(dcdfreq))
        return error
    try:
        locale.atoi(pretamd_min_steps)
        if(locale.atoi(pretamd_min_steps) < 0):
            error.append(
                'preliminary minimization steps needs to be >= 0 : ' + str(pretamd_min_steps))
            return error
    except:
        error.append(
            'preliminary energy minimization steps must be an integer')
        return error

    # if(((int(pretamd_min_steps) % dcdfreq) != 0) or (int(pretamd_min_steps)< 0)):
    #    error.append( 'number of preliminary energy minimization steps must either be ZERO or a mutiple of number of the desired number of TAMD steps: '+str(pretamd_min_steps))
    #    return error

#    flexible_segment_variables = torsion_angle_md.process_input_variables(psegvariables)

    m = sasmol.SasMol(0)
    m.read_pdb(pdbfile)
    resid = m.resid()
    segname = m.segname()
    moltype = m.moltype()

# NOTE:  numerical psegvariables must be integers at this point;
# otherwise, process_input_variables will fail
    flexible_segment_variables = process_input_variables(
        psegvariables, segname)
    print 'fsv = ', flexible_segment_variables

    #d = {u'PEP': [1, [3], [5], 'protein']}
    for key in flexible_segment_variables:
        print flexible_segment_variables[key]
        print 'key = ', key
# Not tested.  Revision to genapp gui on 7/6/17 does not allow separate input of dna_segnames; it will only be populated if moltype is 'dna'.  So this error will not be encountered.
#        if(flexible_segment_variables[key][3] != 'dna'):
#            if key in dna_segnames:
#                error.append('you entered a non DNA segname as DNA : '+str(dna_segnames))
#                return error
        number_of_flexible_regions = flexible_segment_variables[key][0]
#        print 'number_of_flexible_regions: ', number_of_flexible_regions

        if(number_of_flexible_regions < 1):
            error.append(
                'must have at least one flexible region defined for segment ' + str(key))
            return error

        if(len(flexible_segment_variables[key][1]) != number_of_flexible_regions):
            error.append(
                'number of residue ranges does not match number of flexible regions')
            return error

        this_moltype = flexible_segment_variables[key][3]
        print 'this_moltype: ', this_moltype
        if this_moltype == 'rna':
            if 'dna' in moltype:
                print 'there are some dna residues in this pdb file'
                error.append('molecule type ' + flexible_segment_variables[key][
                             3] + ' provided for segment ' + str(key) + ' does not match the pdb file')
                return error
        if this_moltype not in moltype:
            error.append('molecule type ' + flexible_segment_variables[key][
                         3] + ' provided for segment ' + str(key) + ' does not match the pdb file')
            return error

        seg_mol = sasmol.SasMol(0)
        basis = 'segname[i] == "' + str(key) + '"'
        print 'basis = ', basis
        error, mask = m.get_subset_mask(basis)
        if(len(error) > 0):
            error.append('segment name not found in pdb file = ' + str(key))
#            print 'mask = ',mask
#            print 'basis = ',basis
            return error

        error = m.copy_molecule_using_mask(seg_mol, mask, 0)

        # retrieve residue ID for this segment
        resid = seg_mol.resid()
        print 'resid = ', resid
        sys.stdout.flush()
        number_aa = resid[-1] - resid[0] + 1
        print 'number_aa: ', number_aa
        low_residues = []
        num_cont = []
        high_residues = []
        for i in xrange(number_of_flexible_regions):
            this_low = flexible_segment_variables[key][1][i]
            this_numcont = flexible_segment_variables[key][2][i]
#            print 'this_low = ',this_low
#            print 'this_numcont = ',this_numcont
            low_residues.append(this_low)
            high_residues.append(this_low + this_numcont)
            if(this_low < 2):
                error.append(
                    'low residue can not include the n-terminus, reslow = ' + str(this_low))
                return error
            # this check won't work for dsDNA in which 2nd chain numbering
            # doesn't start at resid=1; disable for 'dna'
            if(this_low + this_numcont > number_aa - 1 and this_moltype != 'dna'):
                error.append('your low residue plus number contiguous exceeds the number of amino acids-1 (' +
                             str(number_aa - 1) + '), reslow = ' + str(this_low) + ' numcont = ' + str(this_numcont))
                return error
            if this_low not in resid:
                error.append('Input pdb file, ' + str(pdbfile) + ' does not have residue, ' +
                             str(this_low) + ' for segment name ' + str(key))
                return error
            if this_low + this_numcont not in resid:
                error.append('Input pdb file, ' + str(pdbfile) + ' does not have residue, ' +
                             str(this_low + this_numcont) + ' for segment name ' + str(key))
                return error
            for j in xrange(this_low, this_low + this_numcont):
                if j not in resid:
                    error.append('Input pdb file, ' + str(pdbfile) + ' does not have residue, "' + str(
                        j) + '" for segment name ' + str(key) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
                    return error

        for i in xrange(number_of_flexible_regions):
            try:
                if high_residues[i] > low_residues[i + 1]:
                    error.append('residue ranges overlap, low residue (i+1) = ' + str(
                        low_residues[i + 1]) + ' high residue (i) = ' + str(high_residues[i]))
                    return error
            except:
                pass

    return error
