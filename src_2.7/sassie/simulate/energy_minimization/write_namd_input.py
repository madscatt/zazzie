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

#       WRITE_NAMD_INPUT
#
#       10/1/2007       --      initial coding                  :       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        WRITE_NAMD_INPUT is the function to write a generic input
	file for the program NAMD (version 2.5).

	This module is called by MINIMIZE.PY

	Writes input file for vacuum minimization @ 300K only.

	REFERENCE:

    	J. C. Phillips et al.
    	Journal of Computational Chemistry  26  1781-1802  (2005)

'''


def write_namd_input(inputfilename, nsteps, dcdfreq, pdb, psf, odcd, toppar, md, mdsteps, dielect, temperature):
    '''
    WRITE_NAMD_INPUT is the function to write a generic input
    file for the program NAMD (version 2.5).


    INPUT:  variable descriptions:

            nsteps:			number of minimization steps
            dcdfreq:		frequency to write dcd file to disk
            pdb:			name of pdb file with coordinates to minimize
            psf:			path and name of psf file
            odcd:			name of final minimized structure dcd file
            toppar:			path and name of topology file		 

            md:			flag to run md before minimization (0==no, 1==yes)
            mdsteps:		number of 1 fs timesteps
            dielect:		solvent dielectric constant
            temperature:		temperature


    OUTPUT:
            inputfilename:          name of input file that function will write

    '''

    outfile = open(inputfilename, 'w')

    coordst = 'coordinates             ' + pdb
    psfst = 'structure                 ' + psf
    outst = 'outputname                junk'
    dcdnst = 'DCDfile                  ' + odcd
    dcdfreq = 'DCDfreq		 ' + dcdfreq
    dielectst = 'dielectric		' + str(dielect)
    minst = 'minimize			' + nsteps

    if(md != 0):

        timest = 'timestep                1.0'
        stepst = 'stepspercycle           20'

        elecst = 'fullElectFrequency      2'
        vdwst = 'nonbondedFreq           2'

        tempst = 'temperature		' + str(temperature)
        mdstepsst = 'run	' + str(mdsteps)

        lan1st = 'langevin                on'
        lan2st = 'langevinDamping         1'
        lan3st = 'langevinTemp            ' + str(temperature)
        lan4st = 'langevinHydrogen        on'

    else:
        tempst = 'temperature 		' + str(300.0)

    outfile.write('%s\n' % (coordst))
    outfile.write('%s\n' % (psfst))
    outfile.write('%s\n' % (tempst))
    #outfile.write('%s\n' % ('temperature             300'))
    outfile.write('%s\n' % ('paratypecharmm          on'))
    parmfile = string.split(toppar, ',')
    for this_parmfile in parmfile:
        topst = 'parameters       ' + this_parmfile
        outfile.write('%s\n' % (topst))

    outfile.write('%s\n' % ('exclude                 scaled1-4'))
    outfile.write('%s\n' % ('1-4scaling              1.0'))
    outfile.write('%s\n' % ('switching               on'))
    outfile.write('%s\n' % ('switchdist              8'))
    outfile.write('%s\n' % ('cutoff                  10'))
    outfile.write('%s\n' % ('pairlistdist            12'))
    outfile.write('%s\n' % ('margin                  2'))
    outfile.write('%s\n' % ('rigidBonds              all'))
    outfile.write('%s\n' % ('rigidTolerance          0.00001'))
    outfile.write('%s\n' % ('rigidIterations         500'))
    outfile.write('%s\n' % ('outputenergies          10'))
    outfile.write('%s\n' % (outst))
    outfile.write('%s\n' % ('binaryoutput            no'))
    outfile.write('%s\n' % (dcdnst))
    outfile.write('%s\n' % (dcdfreq))
    outfile.write('%s\n' % (dielectst))

    if(md != 0):
        outfile.write('%s\n' % (timest))
        outfile.write('%s\n' % (stepst))
        outfile.write('%s\n' % (elecst))
        outfile.write('%s\n' % (vdwst))

        outfile.write('%s\n' % (lan1st))
        outfile.write('%s\n' % (lan2st))
        outfile.write('%s\n' % (lan3st))
        outfile.write('%s\n' % (lan4st))
        outfile.write('%s\n' % (minst))
        outfile.write('%s\n' % (mdstepsst))

    if(md == 0 or md == 2):
        outfile.write('%s\n' % (minst))

    outfile.close()

    return

if __name__ == '__main__':

    write_namd_input()
