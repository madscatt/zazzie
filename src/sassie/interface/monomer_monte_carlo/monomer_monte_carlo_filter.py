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
import sasmol.system as system
import sassie.simulate.constraints.constraints as constraints
import sassie.interface.input_filter as input_filter


def check_protein(variables, eflag, monflag, **kwargs):

    runname = variables['runname'][0]
    dcdfile = variables['dcdfile'][0]
    path = variables['path'][0]
    pdbfile = variables['pdbfile'][0]
    trials = variables['trials'][0]
    goback = variables['goback'][0]
    temp = variables['temp'][0]
    moltype = variables['moltype'][0]
    numranges = variables['numranges'][0]
    dtheta = variables['dtheta'][0]
    reslow = variables['reslow'][0]
    numcont = variables['numcont'][0]
    lowres1 = variables['lowres1'][0]
    highres1 = variables['highres1'][0]
    basis = variables['basis'][0]
    cutoff = variables['cutoff'][0]
    lowrg = variables['lowrg'][0]
    highrg = variables['highrg'][0]
    zflag = variables['zflag'][0]
    zcutoff = variables['zcutoff'][0]
    cflag = variables['cflag'][0]
    confile = variables['confile'][0]
    nonbondflag = variables['nonbondflag'][0]
    nonbondscale = variables['nonbondscale'][0]
    psffilepath = variables['psffilepath'][0]
    psffilename = variables['psffilename'][0]
    parmfilepath = variables['parmfilepath'][0]
    parmfilename = variables['parmfilename'][0]
    plotflag = variables['plotflag'][0]
    directedmc = variables['directedmc'][0]

    error = []
    error = input_filter.check_name(runname)
    if(error != []):
        return error

    if "no_file_check" not in kwargs:
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

    if(path != ""):
        if path[-1] != "/":
            pdbfile = path + '/' + pdbfile
        else:
            pdbfile = path + pdbfile
    else:
        pdbfile = path + pdbfile

    ev, value = input_filter.check_pdb_dcd(pdbfile, 'pdb')

    if(ev == 0):
        error.append('input pdb file, ' + pdbfile[3:] + ', does not exist')
        return error
    if(value == 0):
        error.append('input pdb file, ' +
                     #pdbfile[3:] + ', is not a valid pdb file')
                     pdbfile + ', is not a valid pdb file')
        return error

    if(trials < 1):
        error.append('trials = ' + str(trials) + '?')
        return error
    elif(temp < 0):
        error.append('use a positive temperature, temperature = ' + str(temp))
        return error
    elif(moltype != 'protein' and moltype != 'rna'):
        error.append(
            'only protein and rna backbone dihedral move sets are defined, you entered : ' + str(moltype))
        return error
    elif(cutoff < 0.001):
        error.append('use a larger cutoff value, cutoff = ' + str(cutoff))
        return error
    elif(zflag != 0 and zflag != 1):
        error.append(
            'ERROR in Z coordinate filter selection: zflag == 0 for "no" and 1 for "yes", zflag = ' + str(zflag))
        return error
    elif(cflag != 0 and cflag != 1):
        error.append(
            'ERROR in atomic constraints selection: cflag == 0 for "no" and 1 for "yes", cflag = ' + str(cflag))
        return error
    elif(cflag == 1):
        err = input_filter.check_file_exists(confile)
        if(err != []):
            lerr = ['ERROR in constraint filename selection: ']
            lerr.append(err)
            error.append(lerr[0] + err[0])
            return error
        filter_flag = 1
        m1 = system.Molecule(0)
        m1.read_pdb(pdbfile)
        err = constraints.read_constraints(m1, confile, filter_flag)
        if(err != []):
            error.append(err[0])
            return error

    m1 = system.Molecule(0)
    m1.read_pdb(pdbfile, fastread=True)
    name = m1.name()
    if(m1.number_of_moltypes() != 1):
        error.append(
            'ERROR: your PDB structure has more than one molecular type')
        return error
    elif(moltype != m1.moltype()[0]):
        error.append('ERROR: your PDB structure has been identified as ' +
                     m1.moltype()[0] + ' but you entered ' + moltype)
        return error
    elif(nonbondflag != 0 and nonbondflag != 1):
        error.append(
            'nonbondflag == 0 for "no" and 1 for "yes", nonbondflag = ' + str(nonbondflag))
        return error
    elif(lowrg > highrg):
        error.append('low Rg cutoff is larger than high Rg cutoff, lowrg = ' +
                     str(lowrg) + ' highrg = ' + str(highrg))
        return error
    elif(lowrg < 0 or highrg < 0):
        error.append('Rg cutoffs need to be >= zero, lowrg = ' +
                     str(lowrg) + ' highrg = ' + str(highrg))
        return error
    elif(len(dtheta) != numranges):
        error.append('the number of dtheta values does not match the number of ranges, dtheta = ' +
                     str(dtheta) + ' numranges = ' + str(numranges))
        return error
    elif(len(reslow) != numranges):
        error.append('the number of low residue values does not match the number of ranges, lowres = ' +
                     str(reslow) + ' numranges = ' + str(numranges))
        return error
    elif(len(numcont) != numranges):
        error.append('the number of contiguous residue values does not match the number of ranges, contiguous residues = ' +
                     str(numcont) + ' numranges = ' + str(numranges))
        return error

    for th in dtheta:
        if(th > 180.0):
            dtheta[th] = 180.0
        elif(th < 0.0):
            dtheta[th] = 0.0

    locvariables = ['resid']
    value, result = input_filter.get_pdb_stats(pdbfile, locvariables)
    resid = list(map(int, result[0]))

# if(resid[0] != 1):
#error.append('amino acid residues in starting pdbfile '+pdbfile+' must start at resid = 1 : '+str(resid[0]))
# return error

    number_aa = resid[-1] - resid[0] + 1

    for i in range(resid[0], number_aa):
        #	ti=i+1
        ti = i
        if ti not in resid:
            error.append('amino acid ' + str(ti) +
                         ' is missing from pdbfile' + pdbfile)
            print('amino acid ' + str(ti) + ' is missing from pdbfile' + pdbfile)
            return error

    for j in range(numranges):
        if (reslow[j] not in resid):
            error.append('Input pdb file, ' + str(pdbfile) + ' does not have low residue amino acid, "' + str(
                reslow[j]) + '" for segment number ' + str(j) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
            return error
        elif (reslow[j] + numcont[j] not in resid):
            error.append('Input pdb file, ' + str(pdbfile) + ' does not have residue amino acid, "' + str(reslow[
                         j] + numcont[j]) + '" for segment number ' + str(j) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
            return error
        elif (numcont[j] <= 0):
            error.append('The number of contiguous residues "' +
                         str(numcont[j]) + '" should be greater than 0!')
            return error
        elif (not ((lowres1 > (reslow[j] + numcont[j])) or (highres1 < reslow[j]))):
            error.append('alignment and flexible ranges should not overlap!')
            return error
        elif(reslow[j] < resid[0]):             #NOT TESTED.  First error in this "for" loop will be triggered if low residue number is lower than n-terminal aa.
            error.append(
                'low residue is lower than the n-terminal amino acid number, reslow = ' + str(reslow[j]))
            return error
        elif(reslow[j] + numcont[j] > reslow[j] + number_aa + 1):   #NOT TESTED.  First error in this "for" loop will be triggered before getting to this point. 
            error.append('your residue range exceeds the number of amino acids (' + str(
                number_aa) + '), reslow = ' + str(reslow[j]) + ' numcont = ' + str(numcont[j]))
            return error

    for i in range(numranges - 1):
        if(reslow[i] > reslow[i + 1]):
            error.append(
                'low residue values must increase from low to high, reslow = ' + str(reslow))
            return error
        elif(reslow[i] + numcont[i] > reslow[i + 1]):
            error.append('residue ranges overlap, reslow = ' +
                         str(reslow) + ' numcont = ' + str(numcont))
            return error

    if(reslow[-1] + numcont[-1] > reslow[j] + number_aa + 1):
        error.append('your residue range exceeds the number of amino acids (' +
                     str(number_aa) + '), reslow = ' + str(reslow) + ' numcont = ' + str(numcont))
        return error

    if(lowres1 not in resid):
        error.append('input pdb file, ' + str(pdbfile) + ' does not have low residue amino acid, ' +
                     str(lowres1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
        return error
    elif(highres1 not in resid):
        error.append('input pdb file, ' + str(pdbfile) + ' does not have high residue amino acid, ' +
                     str(highres1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
        return error
    elif(highres1 - lowres1 < 3):
        error.append(
            'alignment basis is too small (less than 3 points) or low residue > high residue')
        return error

    if(nonbondflag == 1 and "no_file_check" not in kwargs):

        ev, rv, wv = input_filter.check_permissions(psffilepath)
        if(ev == 0 or rv == 0 or wv == 0):
            error.append('permission error in psf file path ' + psffilepath +
                         '  [code = ' + str(ev) + str(rv) + str(wv) + ']')
            if(ev == 0):
                error.append('path does not exist')
            elif(rv == 0):
                error.append('read permission not allowed')
            elif(wv == 0):
                error.append('write permission not allowed')
            return error

        ev, rv, wv = input_filter.check_permissions(parmfilepath)
        if(ev == 0 or rv == 0 or wv == 0):
            error.append('permission error in parameter file path ' +
                         parmfilepath + '  [code = ' + str(ev) + str(rv) + str(wv) + ']')
            if(ev == 0):
                error.append('path does not exist')
            elif(rv == 0):
                error.append('read permission not allowed')
            elif(wv == 0):
                error.append('write permission not allowed')
            return error

        psffilename = psffilepath + '/' + psffilename
        ev, value = input_filter.certify_pdb_psf(pdbfile, psffilename)

        if(ev == 0):
            error.append('input psf file, ' +
                         psffilename[3:] + ', does not exist')
            return error
        if(value == 0):
            error.append('psffile ' + psffilename[3:] + ' and pdbfile ' + pdbfile[
                         3:] + ' are not compatible (different number of atoms)')
            return error

        parmfilename = parmfilepath + '/' + parmfilename
        error = input_filter.check_file_exists(parmfilename)

    if(plotflag != 0 and plotflag != 1):
        error.append(
            'plot flag needs to be 1 or 0 ... you entered: ' + str(plotflag))
        return error

    elif(directedmc < 0):
        error.append(
            'directed Monte Carlo needs to be 0 or a float > 0 (the "goal Rg") ... you entered: ' + str(directedmc))
        return error

        # check if overlap basis atom exits

    check_atoms = True
    # if True:

#    print 'basis: ', basis
    try:
        if(basis.lower() == 'backbone' or basis.lower() == 'heavy'):
            if(moltype == 'protein'):
                test_atoms = ["N", "CA", "C"]
            elif(moltype == 'rna'):
                test_atoms = ["P", "O5'", "C5'", "C3'", "O3'"]
        elif(basis.lower() == 'all'):
            m1 = system.Molecule(0)
            m1.read_pdb(pdbfile)
            m1.set_average_vdw()
            atom_vdw = m1.atom_vdw()
            name = m1.name()
            element = m1.element()
            natoms = m1.natoms()
#            for i in xrange(natoms):
#                print i,name[i],element[i],atom_vdw[i][0]
            if any(None in sub_list for sub_list in atom_vdw):          #NOT TESTED. There are no elements in the list that have "None" as a vdw parameter.
                atom_list = ''
                name = m1.name()
                for i in range(natoms):
                    if(atom_vdw[i][0] == None):
                        atom_list += ' ' + name[i]
                error.append(
                    "the following atoms do not have vdW parameters; use different overlap basis : " + atom_list)
                return error
            else:
                check_atoms = False
        else:
            string_basis = string.split(basis, ",")
            test_atoms = []
            for i in range(len(string_basis)):
                test_atoms.append(string_basis[i])

#        print 'check_atoms: ', check_atoms

        if check_atoms:
            for atom in test_atoms:
                if(atom not in name):
                    error.append("overlap basis atom " + basis +
                                 " is not in your PDB file")
                    return error

#   other errors are triggered before this statement is reached; kept exception just in case
    except:
        error.append("cannot parse overlap basis")
        return error

    # check residue topology and atom ordering

    m1 = system.Molecule(0)
    #error = m1.read_pdb(pdbfile,fastread=True,saspdbrx_topology=True)
    if(len(error) > 0):
        error.append(result)
        return error

    return error
