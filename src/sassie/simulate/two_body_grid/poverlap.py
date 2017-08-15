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
import numpy
import scipy.spatial.distance
import scipy
from scipy.weave import converters

sys.path.append('./')

import foverlap

#	POVERLAP
#
#	12/22/09	--	initial coding 				:	jc	
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **

def aaoverlap(coor1,coor2,cut,natom1,natom2):

###	OPEN	Only works for same ATOM overlap, but it works

	check=0
	diff2=(coor1-coor2)*(coor1-coor2)
	dist=numpy.sqrt(diff2.sum(axis=1))
 	if (dist < cut).any():
		check=1

#	n.sqrt(n.sum((a-b)**2,axis=1))
	
	return check

def aboverlap(coor1,coor2,natoms1,natoms2,cut):

###	OPEN	Slow loop: convert to C extension

	check=0
	for i in range(natoms1):
		x1=coor1[i][0] ; y1=coor1[i][1] ; z1=coor1[i][2]
		for j in range(i,natoms2):
			x2=coor2[j][0] ; y2=coor2[j][1] ; z2=coor2[j][2]
			diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
			dist=numpy.sqrt(diff2)
			if(dist<cut):
				check=1
				return check
#
#	timing when accepting all moves and writing data to PDB file
#time =  10.0075511932 : avg =  12.7689574824
#ntrials =  324
#
#	timing when accepting all moves without writing to any file (PDB or DCD)
#time =  1.33096003532 : avg =  4.36706425526
#ntrials =  324
#
#	timing when accepting all moves without writing to any file (PDB or DCD) OR transferring coords
#time =  0.935441017151 : avg =  3.99341740432
#ntrials =  324
#
#
	return check

def saboverlap(coor1,coor2,cut):	
	
	check=0
	dist=scipy.spatial.distance.cdist(coor1,coor2,'euclidean')

#	sh=ipshell()
 	if (dist < cut).any():
		check=1
#
###	OPEN	Re-test spatial.distance.cdist once DCD (or memory/DCD) paradigm is working
###	OPEN	Why doesn't (dist<cut).any() work here?
###	DONE	Compare this result to the C/Fortran subroutine
#	
#	timing when accepting all moves and writing data to PDB file
#time =  9.23265004158 : avg =  8.91779266463
#ntrials =  324
#
#	timing when accepting all moves without writing to any file (PDB or DCD)
#time =  0.556995868683 : avg =  0.536804216879
#ntrials =  324
#
#	timing when accepting all moves without writing to any file (PDB or DCD) OR transferring coords
#time =  0.152242183685 : avg =  0.149324116883
#ntrials =  324
# 
	return check

def waboverlap(coor1,coor2,cut):

	check=0
	natoms1=len(coor1)
	natoms2=len(coor2)
#
###	OPEN	Weave does not compile ... wants some ARCH setting
#
	code="""
		int check;
		double x1,y1,z1,x2,y2,z2,diff2,dist;
		check=0;
		for(int i=1;i<natoms1;i++){
			x1=coor1[i][0];
			y1=coor1[i][1];
			z1=coor1[i][2];
			for(int j=1;j<natoms2;j++){
				x2=coor2[i][0];
				y2=coor2[i][1];
				z2=coor2[i][2];
				diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);	
				dist=sqrt(diff2);
				if(dist<cut){
					check=1;
					break;
				}
			}
			if(check==1){
				break;
			}
		}

		return_val=check;

		"""
	check=scipy.weave.inline(code,['coor1','coor2','natoms1','natoms2','cut'],type_converters=converters.blitz,compiler='gcc')

	return check

def faboverlap(coor1,coor2,cut):
	'''
		f2py -c foverlap.f -m foverlap

	'''	
###	OPEN	Write doc-string for fortran f2py stuff
	check=0

	check=foverlap.foverlap(coor1,coor2,cut)

#	timing when accepting all moves and writing data to PDB file
#time =  9.15971398354 : avg =  8.87105663617
#ntrials =  324
#
#	timing when accepting all moves without writing to any file (PDB or DCD)
#time =  0.465914011002 : avg =  0.448736552839
#ntrials =  324
#
#	timing when accepting all moves without writing to any file (PDB or DCD) OR transferring coords
#time =  0.0609328746796 : avg =  0.0605104057877
#ntrials =  324
#
	return check


