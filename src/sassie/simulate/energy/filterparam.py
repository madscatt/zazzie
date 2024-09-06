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

def getexclusion(atoms,bonds):
	print '\ndetermining 1-2/1-3 exclusion list and 1-4 pair list from bonding'
	exclusionlist=[] ; temp=[] ; temp2=[] ; blist=[] ; bondlist=[]
	natoms=len(atoms)
	nbonds=len(bonds)
	
	for i in xrange(nbonds):
		tatom=bonds[i][0]
		batom=bonds[i][1]
		if(i==0):
			lex=[]
			latm=tatom
			lex.append(batom)
			first=1
		elif(tatom!=latm):
			lex=[]
			lex.append(batom)
			first=1
		if(first==1):
			for j in xrange(i+1,nbonds):
				tatom2=bonds[j][0]
				batom2=bonds[j][1]
				if(tatom2==tatom):
					lex.append(batom2)
			temp.append([tatom,lex])	
			first=0
		latm=tatom
	
# now go through list and complete primary bonds

	for i in xrange(len(temp)):
		tatom=temp[i][0]
		lb=[]
		for j in xrange(len(temp)):
			satom=temp[j][0]
			nb=len(temp[j][1])
			for k in xrange(nb):
				if(temp[j][1][k]==tatom):
					lb.append(satom)
					#print tatom,satom
		lex=[]
		for j in xrange(len(temp[i][1])):
			lex.append(temp[i][1][j])
		for j in xrange(len(lb)):
			lex.append(lb[j])
		
		temp2.append([tatom,lex])		

# now go through list and find all single bond partners
	
	for i in xrange(len(temp2)):
		blist.append(temp2[i][0])

	for i in xrange(natoms):
		atom=str(i+1) ; lex=[]
		if atom in blist:
			for j in xrange(len(temp2)):
				if(atom==temp2[j][0]):
					for k in xrange(len(temp2[j][1])):
						lex.append(temp2[j][1][k])			
		else:
			for j in xrange(len(temp2)):
				for k in xrange(len(temp2[j][1])):
					if(temp2[j][1][k]==atom):
						lex.append(temp2[j][0])

		bondlist.append([atom,lex]) 	### 1-2 bonds 

# now create the exclusion list (excluding self)   THIS IS 1-2 & 1-3

	for i in xrange(len(bondlist)):
		tatom=bondlist[i][0]
		lex=[]
		for j in xrange(len(bondlist[i][1])):
			lex.append(bondlist[i][1][j])

		for j in xrange(len(bondlist[i][1])):
			batom=bondlist[i][1][j]
			for k in xrange(len(bondlist)):
				satom=bondlist[k][0]	
				if(satom==batom):
					for kk in xrange(len(bondlist[k][1])):
						if(bondlist[k][1][kk]!=tatom):
							lex.append(bondlist[k][1][kk])
		exclusionlist.append([tatom,lex])	
	
	nex=0
	for i in xrange(len(exclusionlist)):
		lin=len(exclusionlist[i][1])
		nex=nex+lin	
	print 'number of atom exclusions (1-2 & 1-3 pairs) = ',nex/2
	
	#print

# now get the 1-3 list (difference of 1-2 and 1-2 plus 1-3)

	onethree=[]
	for i in xrange(len(exclusionlist)):
		lex=[]
		for j in xrange(len(exclusionlist[i][1])):
			tatom=exclusionlist[i][1][j]
			if tatom not in bondlist[i][1]:
				lex.append(tatom)
		onethree.append([exclusionlist[i][0],lex])


### then put the 1-2 for each atom in the 1-3 list (this is the overcomplete 1-4) 

	overonefour=[]
	for i in xrange(len(onethree)):
		lex=[]
		for j in xrange(len(onethree[i][1])):
			tatom=onethree[i][1][j]
			# now find the index and the bonded atoms for tatom in the 1-2 list
			for k in xrange(len(bondlist)):
				batom=bondlist[k][0]
				if(tatom==batom):
					for m in xrange(len(bondlist[k][1])):
						if(bondlist[k][1][m] not in lex):
							lex.append(bondlist[k][1][m])
		overonefour.append([exclusionlist[i][0],lex])
					
### then remove the atoms in the original 1-2 list for that atom to get the 1-4 list

	onefourlist=[] ; n14=0
	for i in xrange(len(overonefour)):
		lex=[]
		tatom=overonefour[i][0]
		for k in xrange(len(overonefour[i][1])):
			if(overonefour[i][1][k] not in bondlist[i][1]):
				lex.append(overonefour[i][1][k])
				n14=n14+1
		onefourlist.append([exclusionlist[i][0],lex])

	print 'number of 1-4 pairs = ',n14

	return exclusionlist,onefourlist
