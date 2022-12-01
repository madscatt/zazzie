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
import string,locale
import numpy

def getdihedrals(atoms,dihedrals,pdihedrals):
	dihedralparam=[]

#DIHEDRALS
#
#V(dihedral) = Kchi(1 + cos(n(chi) - delta))
#
#Kchi: kcal/mole
#n: multiplicity
#delta: degrees
#
	notfound=0 ; nfa=[]
        nparam=len(pdihedrals) ; ndihedrals=len(dihedrals)
        for i in xrange(ndihedrals):
		tdihedral=dihedrals[i]
                lex=[]

		tatomnum1=locale.atoi(dihedrals[i][0])-1
		tatomnum2=locale.atoi(dihedrals[i][1])-1
		tatomnum3=locale.atoi(dihedrals[i][2])-1
		tatomnum4=locale.atoi(dihedrals[i][3])-1

		atmname1=atoms[tatomnum1][1]
		atmname2=atoms[tatomnum2][1]
		atmname3=atoms[tatomnum3][1]
		atmname4=atoms[tatomnum4][1]

		'''
		if(i==0):
			print 		
			print 'tdihedral = ',dihedrals[0]
			print 'pdihedrals[0] = ',pdihedrals[0]
			print
			print 'atms = ',atmname1,atmname2,atmname3,atmname4
		'''

#pdihedrals[0] =  ['C', 'CT1', 'NH1', 'C', '0.2000', '1', '180.00']

#H    NH1  C    CT3
		found=0
                for j in xrange(nparam):
                        if(pdihedrals[j][0]=='X' and pdihedrals[j][3]=='X'): 
				if ( (pdihedrals[j][1]==atmname2 and pdihedrals[j][2]==atmname3) or (pdihedrals[j][2]==atmname2 and pdihedrals[j][1]==atmname3) ):
					loc=[]
					loc.append(str(tatomnum1))
					loc.append(str(tatomnum2))
					loc.append(str(tatomnum3))
					loc.append(str(tatomnum4))
					lex.append(loc)	
					loc0=[]	
					for k in xrange(len(pdihedrals[j])):
						if(k>3):
							loc0.append(pdihedrals[j][k])
					lex.append(loc0)	
					dihedralparam.append(lex)
					found=1
					break
                        if(pdihedrals[j][1]=='X' and pdihedrals[j][2]=='X'): 
				if ( (pdihedrals[j][0]==atmname2 and pdihedrals[j][3]==atmname3) or (pdihedrals[j][3]==atmname2 and pdihedrals[j][0]==atmname3) ):
					loc=[]
					loc.append(str(tatomnum1))
					loc.append(str(tatomnum2))
					loc.append(str(tatomnum3))
					loc.append(str(tatomnum4))
					lex.append(loc)	
					loc0=[]	
					for k in xrange(len(pdihedrals[j])):
						if(k>3):
							loc0.append(pdihedrals[j][k])
					lex.append(loc0)	
					dihedralparam.append(lex)
					found=1
					break
                        if(pdihedrals[j][0]==atmname1 and pdihedrals[j][3]==atmname4): 
				if(pdihedrals[j][1]==atmname2 and pdihedrals[j][2]==atmname3):
					loc=[]
					loc.append(str(tatomnum1))
					loc.append(str(tatomnum2))
					loc.append(str(tatomnum3))
					loc.append(str(tatomnum4))
					lex.append(loc)	
					loc0=[]	
					for k in xrange(len(pdihedrals[j])):
						if(k>3):
							loc0.append(pdihedrals[j][k])
					lex.append(loc0)	
					dihedralparam.append(lex)
					found=1
					break
                        if(pdihedrals[j][0]==atmname4 and pdihedrals[j][3]==atmname1): 
				if(pdihedrals[j][1]==atmname3 and pdihedrals[j][2]==atmname2):	
					loc=[]
					loc.append(str(tatomnum1))
					loc.append(str(tatomnum2))
					loc.append(str(tatomnum3))
					loc.append(str(tatomnum4))
					lex.append(loc)	
					loc0=[]	
					for k in xrange(len(pdihedrals[j])):
						if(k>3):
							loc0.append(pdihedrals[j][k])
					lex.append(loc0)	
					dihedralparam.append(lex)
					found=1
					break
		if(found==0):
			this=[atmname1,atmname2,atmname3,atmname4]
			if this not in nfa:
				nfa.append(this)
				notfound=notfound+1
			print 'not found tdihedral = ',dihedrals[i]
			print 'atms = ',atmname1,atmname2,atmname3,atmname4
			#print 'pdihedrals[j] = ',pdihedrals[j]



	print 'assigned ',len(dihedralparam),' dihedral parameters'
	print 'missed ',notfound,' dihedrals parameters'

	for j in xrange(len(nfa)):
		print 'no dihedral parameter found for ',nfa[j]

	if('5' in dihedralparam[0][0]):
		print 'yipee'
		print dihedralparam[0][0]
		print dihedralparam[0][1]

	return dihedralparam

def getangles(atoms,angles,pangles):
        angleparam=[]

#ANGLES
#
#V(angle) = Ktheta(Theta - Theta0)**2
#
#V(Urey-Bradley) = Kub(S - S0)**2
#
#Ktheta: kcal/mole/rad**2
#Theta0: degrees
#Kub: kcal/mole/A**2 (Urey-Bradley)
#S0: A
#
        nparam=len(pangles) ; nangles=len(angles)
        for i in xrange(nangles):
		tangle=angles[i]
                lex=[]

		tatomnum1=locale.atoi(angles[i][0])-1
		tatomnum2=locale.atoi(angles[i][1])-1
		tatomnum3=locale.atoi(angles[i][2])-1

		atmname1=atoms[tatomnum1][1]
		atmname2=atoms[tatomnum2][1]
		atmname3=atoms[tatomnum3][1]

                for j in xrange(nparam):
                        if(pangles[j][0]==atmname1 and pangles[j][2]==atmname3): 
				if(pangles[j][1]==atmname2):
					loc=[]
					loc.append(str(tatomnum1))
					loc.append(str(tatomnum2))
					loc.append(str(tatomnum3))
					lex.append(loc)	
					loc0=[]	
					for k in xrange(len(pangles[j])):
						if(k>2):
							loc0.append(pangles[j][k])
					lex.append(loc0)	
					angleparam.append(lex)
					break
                        elif(pangles[j][0]==atmname3 and pangles[j][2]==atmname1): 
				if(pangles[j][1]==atmname2):	
					loc=[]
					loc.append(str(tatomnum1))
					loc.append(str(tatomnum2))
					loc.append(str(tatomnum3))
					lex.append(loc)	
					loc0=[]	
					for k in xrange(len(pangles[j])):
						if(k>2):
							loc0.append(pangles[j][k])
					lex.append(loc0)	
					angleparam.append(lex)
					
					break
	print 'assigned ',len(angleparam),' angle parameters'

        return angleparam


def getvdw(atoms,pnonbond):
	vdwparam=[]
	simple_vdwparam = []	
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

#['CAY', 'CT3']
#['C', '0.000000', '-0.110000', '2.000000']


	nparam=len(pnonbond) ; natoms=len(atoms)	
	for i in xrange(natoms):
		lex=[] ; simple=[]
		for j in xrange(nparam):
			if(pnonbond[j][0]==atoms[i][1]):
	#			print pnonbond[j]  
				for k in xrange(1,len(pnonbond[j])):
					this_value = locale.atof(pnonbond[j][k])				
					lex.append(this_value)
					if(k == 2 or k == 3):
						simple.append(this_value)	
		#vdwparam.append([atoms[i][1],lex])
		vdwparam.append(lex)
		simple_vdwparam.append(simple)

	simple_vdw = numpy.array(simple_vdwparam,numpy.float32)

	return vdwparam,simple_vdw

