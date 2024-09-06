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
import sys,string,locale,numpy,random

def boltz(vdi,vdf,vvdwi,vvdwf,veli,velf,vangi,vangf,vdihedi,vdihedf,nonbondscale,beta):

	vi=vdi+nonbondscale*(vvdwi+veli+vangi+vdihedi)
	vf=vdf+nonbondscale*(vvdwf+velf+vangf+vdihedf)
	
	check=1
	if(vf<vi):
		check=0
	else:
        	kcal2kj=4.184 # 1 kcal = 4.184 kJ
        	vikjmol=vi*kcal2kj
        	vfkjmol=vf*kcal2kj
        	kjmol2j=1.660540E-21
        	vij=vikjmol*kjmol2j
        	vfj=vfkjmol*kjmol2j
        	earg=(-beta*(vfj-vij))  
        	r=random.random()
		if(abs(earg)>500.0):
			arg=0.0
        	else:
			arg=numpy.exp(earg)
        	if(r<arg):
        		check=0
			print
			print 'vdi = ',vdi
			print 'vvdwi = ',vvdwi
			print 'veli = ',veli
			print 'vi = ',vi
			print 'vdihedi = ',vdihedi
			print
			print 'vdf = ',vdf
			print 'vvdwf = ',vvdwf
			print 'velf = ',velf
			print 'vf = ',vf
			print 'vdihedf = ',vdihedf
			print
			print 'vf - vi = ',vf-vi
			print '-beta*(vfj - vij) = ',-beta*(vfj-vij)
			print 'r = ',r,' arg = ',arg,' earg = ',earg	
			print

	return check

def getswitch(rij):

	ctonnb=10.0
	ctofnb=12.0
	ron=ctonnb
	roff=ctofnb

	if(rij<=ron):
		switch=1.0
	elif(rij>ron and rij<=roff):
		arg1=(roff-rij)*(roff-rij)
		arg2=(roff+2.0*rij-3*ron)
		arg3=(roff-ron)*(roff-ron)*(roff-ron)
		switch=(arg1*arg2)/arg3
	elif(rij>roff):
		switch=0.0

	return switch

def el(coor,charge,exclusionlist):

#	print '\ncalculating electrostatic energy'

	natoms=len(charge)

	elenergy=0.0

#	qi*qj/(e1*rij)

	e1=1.0

	qconv=1.620177E-19  		# eV to C 
	qconv2=qconv*qconv		# q*q --> C^2
	eps=8.854187816E-12		# epsilon --> C^2/(J*m)
	ang2m=1E-10			# angstrom to m
	jtokjpmol=6.022137E+20		# J to KJ/mol
	kjtokcal=1.0/4.184		# KJ to kcal

	#qconv2=1.0
	pi=numpy.pi
	conv=(qconv2/(4.0*pi*eps*ang2m))*jtokjpmol*kjtokcal
	#print 'conv = ',conv
	conv=332
	conv=332.4683
	np=0
	
	stupid = True

	for i in xrange(natoms):
		qi=locale.atof(charge[i])
		x1=coor[i][0] ; y1=coor[i][1] ; z1=coor[i][2]
		for j in xrange(i+1,natoms):
			#if str(j+1) not in exclusionlist[i][1]:
			if stupid: 
				qj=locale.atof(charge[j])
				x2=coor[j][0] ; y2=coor[j][1] ; z2=coor[j][2]
				dx2=(x1-x2)*(x1-x2)	
				dy2=(y1-y2)*(y1-y2)	
				dz2=(z1-z2)*(z1-z2)	
				rij=numpy.sqrt(dx2+dy2+dz2)
				switch=getswitch(rij)
				vij=((qi*qj)/(e1*rij))*switch	
				#vij=(qi*qj*qconv2)/(e1*eps*rij*ang2m)
				#vij=vij*jtokjpmol*kjtokcal
				elenergy=elenergy+vij
				np=np+1
			else:
				pass
				#print 'atom ',i+1,' and atom ',j+1,' not counted in vdw calculation'
	elenergy=elenergy*conv

	return elenergy


def vdw(coor,vdwparam,exclusionlist,onefourlist,calc_flag):
	
#	print '\ncalculating vdw energy'

#NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
#cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5 
#                !adm jr., 5/08/91, suggested cutoff scheme
#!
#!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
#!
#!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
#!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
#!
#!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
#!

	natoms=len(vdwparam)

	vdwenergy=0.0

	offlag=1
	#offlag=0
	
	for i in xrange(natoms):
		flag1=0
		if(len(vdwparam[i])>3):
			flag1=1
			thisonefour=onefourlist[i][1]
		epsi=vdwparam[i][1]
		rmi=vdwparam[i][2]
		x1=coor[i][0] ; y1=coor[i][1] ; z1=coor[i][2]
		for j in xrange(i+1,natoms):
			flag2=0
			#if str(j+1) not in exclusionlist[i][1]:
			if(True):
				if(flag1==1 and len(vdwparam[j])>3):
					flag2=1
				if(calc_flag and offlag==1 and flag1==1 and flag2==1 and (str(j+1) in onefourlist[i][1])):
					epsi=vdwparam[i][4]
					rmi=vdwparam[i][5]
					epsj=vdwparam[j][4]
					rmj=vdwparam[j][5]
				else:
					epsj=vdwparam[j][1]
					rmj=vdwparam[j][2]

				epsilon=numpy.sqrt(epsi*epsj)
				rminij=rmi+rmj
				x2=coor[j][0] ; y2=coor[j][1] ; z2=coor[j][2]
				dx2=(x1-x2)*(x1-x2)	
				dy2=(y1-y2)*(y1-y2)	
				dz2=(z1-z2)*(z1-z2)	
				rij=numpy.sqrt(dx2+dy2+dz2)
				switch=getswitch(rij)
				rij6=(rminij/rij)**6
				rij12=rij6*rij6
				vij=epsilon*(rij12-2.0*rij6)*switch
				vdwenergy=vdwenergy+vij
			else:
				pass
				#print 'atom ',i+1,' and atom ',j+1,' not counted in vdw calcuation'
	#print
	return vdwenergy

def calcnb(coor,charge,exclusionlist,vdwparam,onefourlist):

	vvdw=vdw(coor,vdwparam,exclusionlist,onefourlist)

	vel=el(coor,charge,exclusionlist)

	return vvdw,vel

