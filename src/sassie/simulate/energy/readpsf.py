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
import locale

def getpsf(psffile):
	print('\nreading protein structure file : ',psffile)
	segments=[] ; atoms=[] ; charge=[] ; mass=[] ; bonds=[] ; angles=[] ; dihedrals=[] ; impropers=[] 
	
	infile=open(psffile,'r').readlines()
	nlines=len(infile)

#       1 FALA 1    ALA  CAY  CT3   -0.270000       12.0110           0

	st=infile[2].split() ; num_remarks=locale.atoi(st[0]) ; print('remarks = ',num_remarks)
	offset1=2+num_remarks+2
	st=infile[offset1].split() ; num_atoms=locale.atoi(st[0]) ; print('natoms = ',num_atoms)
	offset2=offset1+num_atoms+2
	for i in range(offset1+1,offset1+1+num_atoms):
		tal=infile[i].split()
		segments.append(tal[1]) ; atoms.append([tal[4],tal[5]]) ; charge.append(tal[6]) ; mass.append(tal[7])

	st=infile[offset2].split() ; num_bonds=locale.atoi(st[0]) ; print('nbonds = ',num_bonds)
	numbondlines=num_bonds/4 
	if(num_bonds%(numbondlines*4)!=0):
		leftover=1 
	else:
		leftover=0 
#	print 'numbondlines = ',numbondlines+leftover 
	for i in range(offset2+1,offset2+1+numbondlines+leftover):
		tbl=infile[i].split() 
		for j in range(0,len(tbl),2):
			bonds.append([tbl[j],tbl[j+1]])

#	print bonds

	offset3=offset2+numbondlines+2+leftover
	st=infile[offset3].split() ; num_angles=locale.atoi(st[0]) ; print('nangles = ',num_angles)
	numanglelines=num_angles/3 
	if(num_angles%(numanglelines*3)!=0):
		leftover=1 
	else:
		leftover=0
	print('numanglelines = ',numanglelines+leftover)

	for i in range(offset3+1,offset3+1+numanglelines+leftover):
		tal=infile[i].split()
		for j in range(0,len(tal),3):
			angles.append([tal[j],tal[j+1],tal[j+2]])

#	print angles
	
	offset4=offset3+numanglelines+2+leftover	
	st=infile[offset4].split() ; num_dihedrals=locale.atoi(st[0]) ; print('dihedrals = ',num_dihedrals)
	numdihedrallines=num_dihedrals/2
	if(numdihedrallines>0):
		if(num_dihedrals%(numdihedrallines*2)!=0):
			leftover=1
		else:
			leftover=0
	else:
		leftover=0
	offset5=offset4+numdihedrallines+2+leftover	
	
	for i in range(offset4+1,offset4+1+numdihedrallines+leftover):
		tdl=infile[i].split()
		for j in range(0,len(tdl),4):
			dihedrals.append([tdl[j],tdl[j+1],tdl[j+2],tdl[j+3]])
#	print dihedrals


	st=infile[offset5].split() ; 
	if(len(st)>0):
		num_impropers=locale.atoi(st[0]) ; print('impropers = ',num_impropers)
	else:
		num_impropers=0
	numimproperlines=num_impropers/2
	if(numimproperlines>0):
		if(num_impropers%(numimproperlines*2)!=0):
			leftover=1
		else:
			leftover=0
	else:
		leftover=0
	offset6=offset5+numimproperlines+2+leftover

	for i in range(offset5+1,offset5+1+numimproperlines+leftover):
		til=infile[i].split()
		for j in range(0,len(til),4):
			impropers.append([til[j],til[j+1],til[j+2],til[j+3]])
#	print impropers
	print()

	return segments,atoms,charge,mass,bonds,angles,dihedrals,impropers

