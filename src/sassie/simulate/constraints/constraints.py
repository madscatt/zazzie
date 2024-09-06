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
import sys
import string
import locale
import numpy
import sasmol.sasmol as sasmol

'''
  	Format for constraints between two object definitions:

		(1)		:	(2)	
 
 	SEG RESID NAME 	: SEG RESID NAME : DISTANCE : TYPE1 : TYPE2
 
 	Where SEG is the name of the segment (up to 4 characters)
 	RESID is the residue number (up to 9999) [[ see below for ranges ]]
 	NAME is the name of the atom (up to 4 characters)
 
 	DISTANCE is in angstroms (int or float is okay)

	TYPE1 & TYPE2 define the interaction between (1) and (2) as being
	between atoms (ATM) or center of mass (COM) of each selection.  One
	needs to provide a definition of each objects.  Thus the possible
	values are ATM ATM, ATM COM, COM ATM, COM COM.
 
 	The actual names (SEG, NAME) and residue numbers (RESID)
 	are from the PDB file that you are using.
 
 	SEG can be from the same or different segments
 
 	You can define as many pair constraint definitions that you
 	need.
 
 
 	Example values:
 
 	SEG RESID : SEG RESID NAME : DISTANCE : COM : ATM
 	SEG : SEG RESID NAME : DISTANCE : ATM : ATM
 	SEG : SEG : DISTANCE : ATM : ATM
 	SEG RESID : SEG : DISTANCE : ATM : ATM
 	
 	You can mix and match the granularity of each component
 	of a distance constraint.
 
###	A range of RESID values can be entered as well.  For
	example to cover residues 33 to 37 in a segment that
 	is within 27 angstroms of another segment one could
 	write
 	
 	SEG1 33-37 : SEG2 : 27.0 : ATM : ATM
 
 	Or, if residues 21-28 in one segment are within 8 angstroms
 	of alpha carbons of residues 1321-1359 the protein 
 	backbone of another segment one could write
 
 	SEG1 21-28 : SEG2 1321-1359 CA : 8.0 : ATM : ATM

	With the ATM ATM selection ALL atoms in residues 21-28 in
	SEG1 have to be within ALL CA atoms in residues 1321-1359.

	If one only needed the center of masses of these selections
	to be used one would write

 	SEG1 21-28 : SEG2 1321-1359 CA : 8.0 : COM : COM

	With the COM COM selection the center of mass of residues 21-28 in
	SEG1 have to be within the center of mass of ALL CA atoms 
	in residues 1321-1359.


 	You can list independent pair constraints in multiple
 	entries, for example, you can merely list the previous two
 	examples in the file (use as many lines as needed).
 
 	SEG1 33-37 : SEG2 : 27.0 : ATM : COM
 	SEG1 21-28 : SEG2 1321-1359 CA : 8.0 : COM : COM

	. . . etc.
 
	NOTE: you can not have a RESID or RESID & NAME w/o SEG
		you can not have a NAME w/o SEG & RESID etc.
	
 	The key delineators are the hypen between resid numbers
 	and the colons between selection 1, selection 2, distance,
 	and the type definitions.

'''

def get_index(seq, f):
	ret = None
	for i in xrange(len(seq)):
		if(seq[i] == f):
			ret=i
	return ret
	
def read_constraints(m1,confile,filter_flag):

	error = []


	input = open(confile,'r').readlines()
	number_of_constraints = len(input)
	
	segname=m1.segname()
	resid=m1.resid()
	name=m1.name()

	element1 = []
	element2 = []
	sdistance = []
	type_array = []

	if(filter_flag == 1):
		unique_seg = []
		unique_name = []
		for i in xrange(len(segname)):
			tsegname = segname[i]
			tname = name[i]
			if(tsegname not in unique_seg):
				unique_seg.append(tsegname)
			if(tname not in unique_name):
				unique_name.append(tname)	

		resid_seg = []

		for i in xrange(len(unique_seg)):
			tsegname = unique_seg[i]
			local_res = []
			for j in xrange(len(resid)):	
				tresid = resid[j]	
				tseg = segname[j]
				if(tseg == tsegname):
					local_res.append(tresid)

			resid_seg.append(local_res)

		print 'unique segnames = ',unique_seg
	
	for i in xrange(number_of_constraints):
		print 'i = ',i
		lin = string.split(input[i])
		seg1 = [] ; resid1 = [] ; name1 = []
		seg2 = [] ; resid2 = [] ; name2 = []
		loc_type_array = []
		colon_count = 0 ; seg_count = 0 ; resid_count = 0 ; type_count = 0
		for j in xrange(len(lin)):
			print 'lin[j] = ',lin[j]
			if(lin[j] == ":"):
				if(colon_count == 0):
					element1.append([seg1,resid1,name1])
					seg_count = 0 ; resid_count = 0
					if(filter_flag == 1):
						for seg in seg1:
							if(seg not in unique_seg):
								error.append(" : LINE "+str(i+1)+" segment "+seg+" listed in constraint file is not in your PDB file")	
								return error
						for n in name1:
							if(n not in unique_name):
								error.append(" : LINE "+str(i+1)+" atom name "+n+" listed in constraint file is not in your PDB file")	
								return error
				elif(colon_count == 1):
					element2.append([seg2,resid2,name2])
					seg_count = 0 ; resid_count = 0
					if(filter_flag == 1):
						for seg in seg2:
							if(seg not in unique_seg):
								error.append(" : LINE "+str(i+1)+" segment "+seg+" listed in constraint file is not in your PDB file")	
								return error
						for n in name2:
							if(n not in unique_name):
								error.append(" : LINE "+str(i+1)+" atom name "+n+" listed in constraint file is not in your PDB file")	
								return error
				colon_count += 1
			elif(colon_count == 0):
				if(seg_count == 0):
					tseg = lin[j]
					seg1.append(tseg)
					seg_count += 1
				elif(resid_count == 0):
					tresid = lin[j]	
					resid1.append(tresid)
					resid_count += 1
				else:
					tname = lin[j]
					name1.append(tname)
			elif(colon_count == 1):
				if(seg_count == 0):
					tseg = lin[j]
					seg2.append(tseg)
					seg_count += 1
				elif(resid_count == 0):
					tresid = lin[j]	
					resid2.append(tresid)
					resid_count += 1
				else:
					tname = lin[j]
					name2.append(tname)
					seg_count = 0 ; resid_count = 0
			elif(colon_count == 2):
				if(filter_flag == 1):
					try:
						local_distance = locale.atof(lin[j])
						if(local_distance <= 0.0):
							error.append(" : LINE "+str(i+1)+" distance value is not appropriate: "+lin[j])
							return error
						sdistance.append(local_distance)
					except:
						error.append(" : LINE "+str(i+1)+" no distance specified or error in line: "+lin[j])
						return error
				else:
					sdistance.append(locale.atof(lin[j]))
						
			elif(colon_count == 3):
				type1 = lin[j]
				print 'type1 = ',type1
				if(filter_flag == 1):
					if (type1 != "ATM" and type1 != "COM"):
						error.append(" : LINE "+str(i+1)+" TYPE1 is not valid (ATM OR COM): "+lin[j])
						return error
				loc_type_array.append(type1)
				type_count += 1
		
			elif(colon_count == 4):
				type2 = lin[j]
				print 'type2 = ',type2
				if(filter_flag == 1):
					if (type2 != "ATM" and type2 != "COM"):
						error.append(" : LINE "+str(i+1)+" TYPE2 is not valid (ATM OR COM): "+lin[j])
						return error
				loc_type_array.append(type2)
				type_array.append(loc_type_array)
				type_count += 1

		if(filter_flag == 1 and type_count != 2):
			error.append(" : LINE "+str(i+1)+" Two type definitions are required for each constraint (ATM OR COM)")
			return error
	
	
	constraint_basis1_array = []
	constraint_basis2_array = []

	for i in xrange(number_of_constraints):
		st1 = 'segname[i] == "'+str(element1[i][0][0])+'"'
		st2 = 'segname[i] == "'+str(element2[i][0][0])+'"'
		if(element1[i][1] != []):
			st = element1[i][1][0]
			if(filter_flag == 1):
				tseg = element1[i][0][0]
				local_resid_seg_index = get_index(unique_seg,tseg)
			if(st.find("-")!= -1):
				first = st[:st.find("-")]
				second = st[st.find("-")+1:]
				if(filter_flag == 1):
					if(locale.atoi(second[:]) <= locale.atoi(first[:])):
						error.append(' : resid values in constraint file for constraint '+str(i)+' are incorrect: second value is equal or less than first')
						error.append(st)
						return error
					if(locale.atoi(second[:]) not in resid_seg[local_resid_seg_index]): 
						error.append(' : resid '+str(second[:])+' is not in segment '+tseg)			
						return error
					if(locale.atoi(first[:]) not in resid_seg[local_resid_seg_index]): 
						error.append(' : resid '+str(first[:])+' is not in segment '+tseg)			
						return error

				st1 = st1 + ' and (resid[i] >= '+first[:]+' and resid[i] <= '+second[:]+')'
			else:
				st1 = st1 + ' and resid[i] == '+st[:]
		if(element2[i][1] != []):
			st = element2[i][1][0]
			if(filter_flag == 1):
				tseg = element2[i][0][0]
				local_resid_seg_index = get_index(unique_seg,tseg)
			if(st.find("-")!= -1):
				first = st[:st.find("-")]
				second = st[st.find("-")+1:]
				if(filter_flag == 1):
					if(locale.atoi(second[:]) <= locale.atoi(first[:])):
						error.append(' : resid values in constraint file for constraint '+str(i)+' are incorrect: second value is equal or less than first')
						error.append(st)
						return error
					if(locale.atoi(second[:]) not in resid_seg[local_resid_seg_index]): 
						error.append(' : resid '+str(second[:])+' is not in segment '+tseg)			
						return error
					if(locale.atoi(first[:]) not in resid_seg[local_resid_seg_index]): 
						error.append(' : resid '+str(first[:])+' is not in segment '+tseg)			
						return error
				st2 = st2 + ' and (resid[i] >= '+first[:]+' and resid[i] <= '+second[:]+')'
			else:
				st2 = st2 + ' and resid[i] == '+st[:]
		if(element1[i][2] != []):
			st1 = st1 + ' and name[i] == "'+str(element1[i][2][0])+'"'
		if(element2[i][2] != []):
			st2 = st2 + ' and name[i] == "'+str(element2[i][2][0])+'"'

		constraint_basis1_array.append(st1)
		constraint_basis2_array.append(st2)

	print
	
	if(filter_flag == 1):
		return error
	else:
		return error, constraint_basis1_array, constraint_basis2_array, sdistance, type_array


def check_constraints(m1,mask_a_array,mask_b_array,distance_array,type_array):

	check = 0
	frame = 0

	number_of_constraints = len(distance_array)

	for i in xrange(number_of_constraints):

		type_a = type_array[i][0]
		type_b = type_array[i][1]

		this_distance = distance_array[i]

		if(type_a == 'COM'):

			sub_m1a=sasmol.SasMol(2)
			error = m1.copy_molecule_using_mask(sub_m1a,mask_a_array[i],0)
			coor_a = sub_m1a.calccom(frame)
		else:
			sub_m1a=sasmol.SasMol(2)
			error = m1.copy_molecule_using_mask(sub_m1a,mask_a_array[i],0)
			coor_a = sub_m1a.coor()[0:]
			
		if(type_b == 'COM'):
			sub_m1b=sasmol.SasMol(3)
			error = m1.copy_molecule_using_mask(sub_m1b,mask_b_array[i],0)
			coor_b = sub_m1b.calccom(frame)
		else:
			sub_m1b=sasmol.SasMol(3)
			error = m1.copy_molecule_using_mask(sub_m1b,mask_b_array[i],0)
			coor_b = sub_m1b.coor()[0:]

		if(type_a == 'COM' and type_b == 'COM'):
			delta_distance = numpy.sqrt(numpy.sum((coor_a - coor_b)**2.0))
			if(delta_distance > this_distance):
				check = 1
				return check

		elif(type_a == 'ATM' and type_b == 'ATM'):
			natoms_a = sub_m1a.natoms()	
			natoms_b = sub_m1b.natoms()	

			for j in xrange(natoms_a):
				this_coor_a = coor_a[0][j]
				for k in xrange(natoms_b):
					this_coor_b = coor_b[0][k]
					delta_distance = numpy.sqrt(numpy.sum((this_coor_a - this_coor_b)**2.0))
					if(delta_distance > this_distance):
						check = 1
						return check

		elif(type_a == 'ATM' and type_b == 'COM'):
			natoms_a = sub_m1a.natoms()	
			
			for j in xrange(natoms_a):
				this_coor_a = coor_a[0][j]
				delta_distance = numpy.sqrt(numpy.sum((this_coor_a - coor_b)**2.0))
				if(delta_distance > this_distance):
					check = 1
					return check

		elif(type_a == 'COM' and type_b == 'ATM'):
			natoms_b = sub_m1b.natoms()	
			
			for j in xrange(natoms_b):
				this_coor_b = coor_b[0][j]
				delta_distance = numpy.sqrt(numpy.sum((coor_a - this_coor_b)**2.0))
				if(delta_distance > this_distance):
					check = 1
					return check

	return check


if __name__=='__main__':

	m1=sasmol.SasMol(0)
	m1.read_pdb('complex0.pdb')
	confile = 'constraints.txt'  
	confile = 'int_constraints.txt'  
	filter_flag = 0

	if(filter_flag == 1):
		error = read_constraints(m1,confile,filter_flag)
		if(error != []):
			print 'CONSTRAINT ERROR: '+str(error[0])
		else:
			print 'NO ERRORS FOUND'
	else:
		error,cb1a,cb2a,sd,type_array = read_constraints(m1,confile,filter_flag)
		for i in xrange(len(cb1a)):
			print cb1a[i]
			print cb2a[i]
			print sd[i]
			print type_array[i]
			print
