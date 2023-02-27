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
import string
import sys

#       PREPEND_NAMD_INPUT
#
#       5/26/2015       --      initial coding                  :       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    PREPEND_NAMD_INPUT is the function to write a generic input
    file for the program NAMD (version 2.9).

    This module is called by ENERGY_MINIMIZATION.PY

    The module reads in a template input file and adds local
    file information to complete the input file

'''


def prepend_namd_input(inputfilename, pdbfile, psffile, odcd, parmfile, external_input_file, velocity_restart_file, extended_system_restart_file):

    extfile = open(external_input_file, 'r').readlines()

    outfile = open(inputfilename, 'w')

    coordst = 'coordinates             ' + pdbfile
    psfst = 'structure                 ' + psffile
    topst = 'parameters                ' + parmfile
    outst = 'outputname                junk'
    dcdst = 'DCDfile                  ' + odcd

    outfile.write('%s\n' % (coordst))
    outfile.write('%s\n' % (psfst))
    outfile.write('%s\n' % ('paratypecharmm          on'))
    outfile.write('%s\n' % (topst))
    outfile.write('%s\n' % (outst))
    outfile.write('%s\n' % (dcdst))

    overwrite_file_list = ['coordinates', 'structure',
                           'paratypecharmm', 'parameters', 'outputname', 'DCDfile']

    if velocity_restart_file != "False":
        velrestart = 'velocities         ' + velocity_restart_file
        outfile.write('%s\n' % (velrestart))
        overwrite_file_list.append('velocities')

    if extended_system_restart_file != "False":
        extendedrestart = 'extendedSystem         ' + extended_system_restart_file
        outfile.write('%s\n' % (extendedrestart))
        overwrite_file_list.append('extendedSystem')

    for line in extfile:
        lin = string.split(line)
        try:
            if lin[0] not in overwrite_file_list:
                # try:
                if True:
                    path, head = os.path.split(lin[1])
                    if (len(path) > 0):
                        outfile.write('%s\t%s\n' % lin[0], head)
                    else:
                        outfile.write('%s' % (line))
                # except:
                else:
                    outfile.write('%s' % (line))
        except:
            pass

    outfile.close()

    return

if __name__ == '__main__':

    prepend_namd_input()
