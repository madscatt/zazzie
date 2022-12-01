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
import os,locale

def getparms(parameter_file):
	print('\nreading parameters from file : ',parameter_file)
	pbonds=[] ; pangles=[] ; pdihedrals=[] ; pimpropers=[] ; pnonbond=[]

	infile=open(parameter_file,'r').readlines()
	nlines=len(infile)

	for i in range(nlines):
		lin=infile[i].split()
		if(len(lin)>0):
			if(lin[0]=='BONDS'):
				bondoffset=i
			if(lin[0]=='ANGLES'):
				angleoffset=i
			if(lin[0]=='DIHEDRALS'):
				dihedraloffset=i
			if(lin[0]=='IMPROPER'):
				improperoffset=i
			if(lin[0]=='NONBONDED'):
				nonbondoffset=i
	
	for i in range(bondoffset,angleoffset):
		tbl=infile[i].split()
		if(len(tbl)>0 and tbl[0][0]!="!"):
			if(len(tbl)>3):
				pbonds.append([tbl[0],tbl[1],tbl[2],tbl[3]])

	for i in range(angleoffset,dihedraloffset):
		tal=infile[i].split()
		if(len(tal)>0 and tal[0][0]!="!"):
			if(len(tal)>4):
				good=0
				for j in range(len(tal)):
					if(tal[j][0]=='!'):
						good=j
						break
					if(j==len(tal)-1):
						good=j+1
						break
				loc=[]
				for j in range(0,good):
					loc.append(tal[j])
				if(len(loc)>4): pangles.append(loc)

	for i in range(dihedraloffset,improperoffset):
		tdl=infile[i].split()
		if(len(tdl)>0 and tdl[0][0]!="!"):
			if(len(tdl)>6):
				pdihedrals.append([tdl[0],tdl[1],tdl[2],tdl[3],tdl[4],tdl[5],tdl[6]])

	for i in range(improperoffset,nonbondoffset):
		til=infile[i].split()
		if(len(til)>0 and til[0][0]!="!"):
			if(len(til)>6):
				pimpropers.append([til[0],til[1],til[2],til[3],til[4],til[5],til[6]])

	for i in range(nonbondoffset+8,nlines):
		tnl=infile[i].split()
		if(len(tnl)>0 and tnl[0][0]!="!"):
			if(len(tnl)>3):
				good=0
				for j in range(len(tnl)):
					if(tnl[j][0]=='!'):
						good=j
						break
					if(j==(len(tnl)-1)):
						good=j+1
						break
				loc=[]
				for j in range(0,good):
					loc.append(tnl[j])
				if(len(loc)>3): pnonbond.append(loc)
	print('found ',len(pbonds),' bond parameters')
	print('found ',len(pangles),' angle parameters')
	print('found ',len(pdihedrals),' dihedral parameters')
	print('found ',len(pimpropers),' improper parameters')
	print('found ',len(pnonbond),' non-bonded parameters\n')
		
	return pbonds,pangles,pdihedrals,pimpropers,pnonbond

