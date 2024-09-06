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
import string
import sys
import readpsf
import readparam
import filterparam
import initialize


def init(psffilepath,psffilename,parmfilepath,parmfilename,simple_flag):

	segments,atoms,charge,mass,bonds,angles,dihedrals,impropers=readpsf.getpsf(psffilepath,psffilename)

	pbonds,pangles,pdihedrals,pimpropers,pnonbond=readparam.getparms(parmfilepath,parmfilename)

	vdwparam,simple_vdwparam=initialize.getvdw(atoms,pnonbond)

	if(not simple_flag):

		'>>> calculating exclusion and one-four lists'
		exclusionlist,onefourlist=filterparam.getexclusion(atoms,bonds)

		angleparam=initialize.getangles(atoms,angles,pangles)

		print 'dihedral parameters were NOT determined'
		#dihedralparam=initialize.getdihedrals(atoms,dihedrals,pdihedrals)

	else:
		exclusionlist = [] ; onefourlist = [] ; angleparam = []
	
	dihedralparam = []

	return charge,exclusionlist,onefourlist,vdwparam,simple_vdwparam,angleparam,dihedralparam



