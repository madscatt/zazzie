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
import sasmol.system as system
import sassie.util.sasutil as sasutil


def check_and_convert_formula(formula_array):

    error = []
    formulas = []

    try:
        number_of_formulas = len(formula_array)
    except:
        error.append('unable to read formula')
        return error, formulas

    for i in xrange(number_of_formulas):

        error, formula_dictionary = sasutil.get_chemical_formula(formula_array[i])

        if(len(error) > 0):
            return error, formulas
        else:
            formulas.append(formula_dictionary)

    return error, formulas

def check_name(filename):
    bad_characters = ["<", ">", "|", "\\",
                      ":", "(", ")", "&", ";", "#", "?", "*"]
    error = []
    for i in xrange(len(filename)):
        character = filename[i]
        if character in bad_characters:
            error.append('file or path : ' + filename +
                         ' has incorrect character : ' + character)
#            print('file or path has incorrect character : ' + character)
            return error
    return error


def check_file_exists(infile):
    error = []

    try:
        error = check_name(infile)
        if(len(error) > 0):
            return error
    except:
        error.append("failed to check infile name")
        return error

    try:
        value = os.path.exists(infile)
    except:
        error.append('file/path : ' + infile + ' does not exist')
        return error

    if(value == 0):
        error.append('file : ' + infile + ' does not exist')
        return error

#	if(os.path.isdir(infile)):
#		error.append('file : '+infile+' is a directory!')
#		return error

    return error


def check_exe(exe):
    error = []
    if(os.path.isdir(exe)):
        error.append('Executable file : ' + exe + ' is a directory!')
        return error
    elif not os.path.isfile(exe):
        error.append('Executable file : ' + exe + ' is not a file')
        return error
    elif not os.access(exe, os.X_OK):
        error.append('Executable file : ' + exe + ' is not accessible')
        return error

    return error


def type_check_and_convert(svariables):

    error = []
    variables = {}

    for key in svariables:

        if (svariables[key][1] == 'string'):
            variables[key] = svariables[key]

        elif (svariables[key][1] == 'boolean'):
            variables[key] = svariables[key]

        elif (svariables[key][1] == 'float'):
            try:
                dum = locale.atof(svariables[key][0])
                variables[key] = (dum, 'float')
            except:
                error.append(
                    key + ' is not a ' + str(svariables[key][1]) + ' : ' + str(svariables[key][0]))
                return error, variables

        elif (svariables[key][1] == 'int'):
            try:
                dum = locale.atoi(svariables[key][0])
                variables[key] = (dum, 'int')
            except:
                error.append(key + ' is not a ' +
                             svariables[key][1] + ' : ' + svariables[key][0])
                return error, variables

        elif (svariables[key][1] == 'float_array'):
            value = svariables.get(key)

            try:
                lin = string.split(value[0], ',')

                duma = []
                for x in xrange(len(lin)):
                    try:
                        dum = locale.atof(lin[x])
                        duma.append(dum)
                    except:
                        error.append(
                            key + ' is not a ' + svariables[key][1] + ' : ' + svariables[key][0])
                        return error, variables
                variables[key] = (duma, 'float_array')
            except:
                error.append(key + ': could not read array of values')
                return error, variables

        elif (svariables[key][1] == 'int_array'):
            value = svariables.get(key)

            try:
                lin = string.split(value[0], ',')

                duma = []
                for x in xrange(len(lin)):
                    try:
                        dum = locale.atoi(lin[x])
                        duma.append(dum)
                    except:
                        error.append(
                            key + ' is not a ' + svariables[key][1] + ' : ' + svariables[key][0])
                variables[key] = (duma, 'int_array')
            except:
                error.append(key + ': could not read array of values')
                return error, variables

    return error, variables


def check_permissions(path):

    try:
        existvalue = os.access(path, os.F_OK)
    except:
        existvalue = False
    try:
        readvalue = os.access(path, os.R_OK)
    except:
        readvalue = False
    try:
        writevalue = os.access(path, os.W_OK)
    except:
        writevalue = False

    return existvalue, readvalue, writevalue


def check_binary(filename):

    textchars = ''.join(
        map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

    flag = is_binary_string(open(filename).read(1024))

    return flag


def check_pdb_dcd(infile, filetype):

    fileexist = 0
    value = 0
    try:
        fileexist = os.path.isfile(infile)
        if(fileexist):
            binary = check_binary(infile)
            print('binary = ', binary)
            test_mol = system.Molecule(0)
            fileexist = 1
            if(filetype == 'pdb' and not binary):
                test_mol.read_pdb(infile, fastread=True)
            elif(filetype == 'dcd' and binary):
                test_mol.read_single_dcd_step(infile, 0)
            else:
                return fileexist, value
            value = 1
        else:
            return fileexist, value
    except:
        value = 0

    return fileexist, value


def certify_pdb_pdb(pdbfile1, pdbfile2):

    fileexist = 0
    value = 0
    try:
        fileexist1 = os.path.isfile(pdbfile1)
        fileexist2 = os.path.isfile(pdbfile2)
        if(fileexist1 and fileexist2):
            fileexist = 1
            pdbmol1 = system.Molecule(0)
            pdbmol2 = system.Molecule(1)
            try:
                pdbmol1.read_pdb(pdbfile1, fastread=True)
                pdbmol2.read_pdb(pdbfile2, fastread=True)
                name1 = pdbmol1.name()
                name2 = pdbmol2.name()
                if(name1 == name2):
                    value = 1
            except:
                value = 0
        else:
            return fileexist, value
    except:
        value = 0

    return fileexist, value


def read_psf_file(psffile):

    segments = []
    names = []

    infile = open(psffile, 'r').readlines()
    nlines = len(infile)

#       1 FALA 1    ALA  CAY  CT3   -0.270000       12.0110           0

    st = string.split(infile[2])
    num_remarks = locale.atoi(st[0])
    print('remarks = ', num_remarks)
    offset1 = 2 + num_remarks + 2
    st = string.split(infile[offset1])
    natoms = locale.atoi(st[0])
    print('natoms = ', natoms)
    offset2 = offset1 + natoms + 2
    for i in xrange(offset1 + 1, offset1 + 1 + natoms):
        tal = string.split(infile[i])
        segments.append(tal[1])
        names.append(tal[4])

    return natoms, names


def certify_dcd_psf(dcdfile, psffile):
    '''
    This method checks that the number of atoms in the psf file
    is equal to the number of atoms in the dcd file.

    The method assumes that the psf and dcd files exist and are
    readable.
    '''

    fileexist = 0
    value = 0
    try:
        fileexist = os.path.isfile(psffile)
        if(fileexist):
            fileexist = 1
            try:
                natoms_psf, names_psf = read_psf_file(psffile)
                dcdmol = system.Molecule(1)

                dcdfile = dcdmol.open_dcd_read(dcdfile)
                natoms_dcd = dcdfile[1]

                if(natoms_psf == natoms_dcd):
                    value = 1
            except:
                value = 0
        else:
            return fileexist, value
    except:
        value = 0

    return fileexist, value


def certify_pdb_psf(pdbfile, psffile):

    fileexist = 0
    value = 0
    try:
        fileexist = os.path.isfile(psffile)
        if(fileexist):
            fileexist = 1
            try:
                natoms_psf, names_psf = read_psf_file(psffile)
                pdbmol = system.Molecule(1)
                pdbmol.read_pdb(pdbfile, fastread=True)
                natoms_pdb = pdbmol.natoms()
                names_pdb = pdbmol.name()
                # if((natoms_pdb == natoms_psf) and (names_pdb == names_psf)):
                if((natoms_pdb == natoms_psf)):
                    value = 1
            except:
                value = 0
        else:
            return fileexist, value
    except:
        value = 0

    return fileexist, value


def certify_pdb_dcd(pdbfile, dcdfile):
    '''
    This method checks that the number of atoms in the pdb file
    is equal to the number of atoms in the dcd file.

    The method assumes that the pdb and dcd files exist and are
    readable.
    '''
    value = 0
    try:
        pdbmol = system.Molecule(0)
        dcdmol = system.Molecule(1)

        pdbmol.read_pdb(pdbfile, fastread=True)
        natoms_pdb = pdbmol.natoms()

        dcdfile = dcdmol.open_dcd_read(dcdfile)
        natoms_dcd = dcdfile[1]

        if(natoms_pdb == natoms_dcd):
            value = 1
    except:
        value = 0

    return value


def get_pdb_stats(filename, variables):
    value = 0
    try:
        a = system.Molecule(0)
        a.read_pdb(filename, fastread=True)
        result = []
        try:
            for i in xrange(len(variables)):
                if(variables[i] == 'atom'):
                    result.append(a.atom())
                elif(variables[i] == 'index'):
                    result.append(a.index())
                elif(variables[i] == 'name'):
                    result.append(a.name())
                elif(variables[i] == 'loc'):
                    result.append(a.loc())
                elif(variables[i] == 'resname'):
                    result.append(a.resname())
                elif(variables[i] == 'chain'):
                    result.append(a.chain())
                elif(variables[i] == 'resid'):
                    result.append(a.resid())
                elif(variables[i] == 'rescode'):
                    result.append(a.rescode())
                elif(variables[i] == 'x'):
                    result.append(coor[0, :, 0]())
                elif(variables[i] == 'y'):
                    result.append(coor[0, :, 1]())
                elif(variables[i] == 'z'):
                    result.append(coor[0, :, 2]())
                elif(variables[i] == 'occupancy'):
                    result.append(a.occupancy())
                elif(variables[i] == 'beta'):
                    result.append(a.beta())
                elif(variables[i] == 'segname'):
                    result.append(a.segname())
                elif(variables[i] == 'element'):
                    result.append(a.element())
                elif(variables[i] == 'charge'):
                    result.append(a.charge())
                elif(variables[i] == 'moltype'):
                    result.append(a.moltype())
            value = 1

        except:
            value = 0
            result = None
    except:
        value = 0
        result = None

    return value, result


def get_pdb_complex_stats(filename, segname, variables):
    value = 0
    try:
        o = system.Molecule(0)
        o.read_pdb(filename, fastread=True)
        seg_filter = 'segname[i] == "' + segname.strip() + '"'
        error, seg_mask = o.get_subset_mask(seg_filter)
        a = system.Molecule(1)
        error = o.copy_molecule_using_mask(a, seg_mask, 0)
        result = []
        try:
            for i in xrange(len(variables)):
                if(variables[i] == 'atom'):
                    result.append(a.atom())
                elif(variables[i] == 'index'):
                    result.append(a.index())
                elif(variables[i] == 'name'):
                    result.append(a.name())
                elif(variables[i] == 'loc'):
                    result.append(a.loc())
                elif(variables[i] == 'resname'):
                    result.append(a.resname())
                elif(variables[i] == 'chain'):
                    result.append(a.chain())
                elif(variables[i] == 'resid'):
                    result.append(a.resid())
                elif(variables[i] == 'rescode'):
                    result.append(a.rescode())
                elif(variables[i] == 'x'):
                    result.append(coor[0, :, 0]())
                elif(variables[i] == 'y'):
                    result.append(coor[0, :, 1]())
                elif(variables[i] == 'z'):
                    result.append(coor[0, :, 2]())
                elif(variables[i] == 'occupancy'):
                    result.append(a.occupancy())
                elif(variables[i] == 'beta'):
                    result.append(a.beta())
                elif(variables[i] == 'segname'):
                    result.append(a.segname())
                elif(variables[i] == 'element'):
                    result.append(a.element())
                elif(variables[i] == 'charge'):
                    result.append(a.charge())
                elif(variables[i] == 'moltype'):
                    result.append(a.moltype())
            value = 1

        except:
            value = 0
            result = None
    except:
        value = 0
        result = None

    return value, result
