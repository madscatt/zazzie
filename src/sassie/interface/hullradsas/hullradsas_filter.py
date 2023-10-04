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
import sassie.interface.input_filter as input_filter

def check_hullradsas(variables, **kwargs):

    run_name = variables['run_name'][0]
    pdbfile = variables['pdbfile'][0]

    error = []
    error = input_filter.check_name(run_name)
    if(error != []):
        return error

    print('pdbfile = ',pdbfile)

    ev, value = input_filter.check_pdb_dcd(pdbfile, 'pdb')

    print('ev = ',ev)
    print('value = ',value)

    if(ev == 0):
        error.append('input pdb file, ' + pdbfile[3:] + ', does not exist')
        return error
    if(value == 0):
        error.append('input pdb file, ' +
                     #pdbfile[3:] + ', is not a valid pdb file')
                     pdbfile + ', is not a valid pdb file')
        return error

    m1 = system.SasMol(0)
    error = m1.read_pdb(pdbfile, fastread=True)

    if(len(error) > 0):
        error.append(result)
        return error

    return error
