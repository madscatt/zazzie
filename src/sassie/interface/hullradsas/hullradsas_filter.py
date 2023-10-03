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

def check_hullradsas(variables, eflag, monflag, **kwargs):

    runname = variables['runname'][0]
    path = variables['path'][0]
    pdbfile = variables['pdbfile'][0]

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
                     pdbfile[3:] + ', is not a valid pdb file')
        return error

    m1 = sasmol.SasMol(0)
    error = m1.read_pdb(pdbfile, fastread=True)

    if(len(error) > 0):
        error.append(result)
        return error

    return error
