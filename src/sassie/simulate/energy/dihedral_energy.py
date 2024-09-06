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
import os,sys,string,locale,random
import numpy,math

#       dihedral_energy.py:  deal with dihedral energies
#
#	01/10/07	--	initial coding						: 	jc
#	08/29/08	--	re-write, added terminal parameter options	:	jc
#	01/12/11	--	added sasmol support					:	jc
#

def calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object):
	search=1 ; vdi=1E10 ; vdf=1E10
	pi=numpy.pi
	angle_to_radians=pi/180.0
	thetarad=theta*angle_to_radians
	ithetarad=itheta*angle_to_radians
	ai=angle_index
	number_of_parm = len(parm[ai])/3

	iv = 0.0 ; v = 0.0
	for i in xrange(number_of_parm):
		offset = i*3
		k_ang = parm[ai][0+offset]
		n_ang = parm[ai][1+offset]
		delta_ang = (parm[ai][2+offset])*angle_to_radians

		iv += k_ang*(1.0+math.cos(n_ang*ithetarad-delta_ang))
		v += k_ang*(1.0+math.cos(n_ang*thetarad-delta_ang))

	if(nonbondflag==0):
		if(v<iv):
			search=0
			vdi=iv
			vdf=v	
		else:
			kcal2kj=4.184 # 1 kcal = 4.184 kJ
			vkjmol=v*kcal2kj
			ivkjmol=iv*kcal2kj
			kjmol2j=1.660540E-21
			vj=vkjmol*kjmol2j
			ivj=ivkjmol*kjmol2j
			earg=(-beta*(vj-ivj))	
			if(seed_object != -1):
				r = seed_object.rand()
			else:
				r=random.random()
			arg=numpy.exp(earg)
			if(r<arg):
				search=0
				vdi=iv
				vdf=v	
	
		return search,vdi,vdf
	else:
		vdi=iv
		vdf=v	
		return vdi,vdf

def getphi(p,t):

	if(t=="PRO"):
		if(p=="NULLN"):
			type=4	# NULLN-PRO
		elif(p=="MAN"):
			type=99	# MAN-PRO
		else:
			type=3	# PRO-PRO, GLY-PRO, X-PRO
	elif(t=="GLY"):
		if(p=="NULLN"):
			type=6	# NULLN-GLY
		elif(p=="MAN"):
			type=98	# MAN-GLY
		else:
			type=5	# GLY-GLY, PRO-GLY, X-GLY
	else:
		if(p=="NULLN"):
			type=2	# NULLN-X
		elif(p=="MAN"):
			type=97	# MAN-X
		else:
			type=1	# X-X

	return type	

def getpsi(t,n):

	if(t=="PRO"):
		if(n=="PRO"):
			type=3	# PRO-PRO
		elif(n=="NULLC"):
			type=5	# PRO-NULLC
		elif(n=="MAN"):
			type=99	# PRO-MAN
		else:
			type=4	# PRO-GLY, PRO-X
	elif(t=="GLY"):
		if(n=="PRO"):
			type=7	# GLY-PRO
		elif(n=="NULLC"):
			type=8	# GLY-NULLC
		elif(n=="MAN"):
			type=98	# GLY-MAN
		else:
			type=6	# GLY-GLY, GLY-X
	else:
		if(n=="NULLC"):	
			type=2	# X-NULLC
		elif(n=="MAN"):
			type=97	# X-MAN
		elif(n=="PRO"):
			type=9	# X-PRO
		else:
			type=1	# X-X

	return type

def getparamphi(type):
	k=0.0 ; n=0.0; delta=0.0
	k2=0.0 ; n2=0.0 ; delta2=0.0
	if(type==1):
		k=0.2; n=1.0; delta=180.0	# C-NH1-CT1-C
	elif(type==2):
		k=1E90; n=0.0 ; delta=180.0	# _-NH3-CT1-C
	elif(type==3):
		k=0.8 ; n=3.0 ; delta=0.0	# C-N-CP1-C
	elif(type==4):
		k=1E90 ; n=0.0 ; delta=0.0	# _-NP-CP1-C
	elif(type==5):
		k=0.2 ; n=1.0 ; delta=180.0	# C-NH1-CT2-C
	elif(type==6):
		k=1E90 ; n=0.0 ; delta=0.0	# _-NH3-CT2-C
	elif(type==97):
		man=1				# MAN-X ---> placeholder
	elif(type==98):
		man=1				# MAN-GLY ---> placeholder
	elif(type==99):
		man=1				# MAN-PRO ---> placeholder
	else:
		print "WRONG PHI TYPE: ",type,":: BUG !!!\n\n"
		print 3/0

	return k,n,delta,k2,n2,delta2

def getparampsi(type):
	k=0.0 ; n=0.0; delta=0.0
	k2=0.0 ; n2=0.0 ; delta2=0.0
	if(type==1):
		k=0.6 ; n=1.0 ; delta=0.0	# NH1-CT1-C-NH1
	elif(type==2):
		k=1E90 ; n=0.0 ; delta=0.0	# NH1-CT1-C-_
	elif(type==3):
		k=0.3 ; n=1.0 ; delta=0.0
		k2=-0.3 ; n2=4.0 ; delta2=0.0	# N-CP1-C-N
	elif(type==4):
		k=0.3 ; n=1.0 ; delta=0.0
		k2=-0.3 ; n2=4.0 ; delta2=0.0	# N-CP1-C-NH1
	elif(type==5):
		k=1E90 ; n=0.0 ; delta=0.0	# N-CP1-C-_
	elif(type==6):
		k=0.6 ; n=1.0 ; delta=0.0	# NH1-CT2-C-NH1
	elif(type==7):
		k=0.4 ; n=1.0 ; delta=0.0	# NH1-CT2-C-N
	elif(type==8):
		k=1E90 ; n=0.0 ; delta=0.0	# NH1-CT2-C-_
	elif(type==9):
		k=0.4 ; n=1.0 ; delta=0.0	# NH1-CT1-C-N
	elif(type==97):
		man=1				# X-MAN ---> placeholder
	elif(type==98):
		man=1				# GLY-MAN ---> placeholder
	elif(type==99):
		man=1				# PRO-MAN ---> placeholder
	else:
		print "WRONG PSI TYPE: ",type,":: BUG !!!\n\n"
		print 3/0

	return k,n,delta,k2,n2,delta2

def getrnaparm(null_alpha,null_epsilon,null_eta):

	tbeta = [0.2, 1, 120.0]
	tgamma = [0.2,4,180.0,0.8,3,180.0,0.4,2,0.0,2.5,1,180.0]
	tdelta = [0.2,4,0.0,0.8,3,180.0]
	tepsilon = [0.6,5,0.0,0.2,4,0.0,0.0,3,180.0,0.4,2,0.0,1.9,1,180.0]	
	
	if(null_alpha == 1):
		talpha = [1E90,0,180.0,1E90,2,180.0,1E90,3,180.0,1E90,0,0.0]
	else:
		talpha = [1.20,1,180.0,0.1,2,180.0,0.1,3,180.0,0.0,6,0.0]
	if(null_epsilon == 1):
		tepsilon = [1E90,5,0.0,1E90,4,0.0,1E90,3,180.0,1E90,2,0.0,1E90,1,180.0]	
	else:
		tepsilon = [0.6,5,0.0,0.2,4,0.0,0.0,3,180.0,0.4,2,0.0,1.9,1,180.0]	
	if(null_eta == 1):
		teta = [1E90,1,180.0,1E90,2,180.0,1E90,3,180.0,1E90,6,0.0]
	else:
		teta = [1.20,1,180.0,0.10,2,180.0,0.1,3,180.0,0.0,6,0.0]
	
	return talpha,tbeta,tgamma,tdelta,tepsilon,teta

def rna_initialization(resalpha,resbeta,resgamma,resdelta,resepsilon,reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput):

	first_resid = first_last_resid[0]
	last_resid = first_last_resid[1]

	for i in xrange(numranges):

		first=rlow[i]		# i.e. 23
		last=rlow[i]+rnum[i]	# i.e. 28=23+5 (5 in this range)

		zero=first-1		# i.e. 22=23-1 for init phi(1)
		totres=last-zero	# 28-22 = 6 : 22(0),23,24,25,26,27,28,29(0)

		alpha = [] ; beta = [] ; gamma = [] ; delta = [] ; epsilon = [] ; eta = []
		
		for j in xrange(totres):
			null_alpha = 0 ; null_epsilon = 0 ; null_eta = 0
			thisres=first+j
			if(thisres==first_resid):
				null_alpha = 1
			if((thisres-1)==last_resid-1):
				null_eta = 1
				null_epsilon = 1

			talpha,tbeta,tgamma,tdelta,tepsilon,teta = getrnaparm(null_alpha,null_epsilon,null_eta)
		
			alpha.append(talpha)
			beta.append(tbeta)
			gamma.append(tgamma)
			delta.append(tdelta)
			epsilon.append(tepsilon)
			eta.append(teta)
	
		resalpha.append([i,alpha])	
		resbeta.append([i,beta])
		resgamma.append([i,gamma])
		resdelta.append([i,delta])
		resepsilon.append([i,epsilon])
		reseta.append([i,eta])

	return


def protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput):
#

	first_resid = first_last_resid[0]
	last_resid = first_last_resid[1]

#  loop over the residues in each range and
#    identify type of Psi and Phi angle for 
#    each residue.  Then assign the dihedral
#    energy constants into a energy_pair array for
#    each residue.

#	rlow[] = array of starting amino acid residues
#	rnum[] = array of number of amino acids in each range

#	testoutput=open("test.txt",'w')
	for i in xrange(numranges):

		first=rlow[i]		# i.e. 23
		last=rlow[i]+rnum[i]	# i.e. 28=23+5 (5 in this range)

		zero=first-1		# i.e. 22=23-1 for init phi(1)
		totres=last-zero	# 28-22 = 6 : 22(0),23,24,25,26,27,28,29(0)

		#                        the underlined 6:  _________________
		# 6 residues in range: 8 residue info. needed
		#
		# [22(phi),23(psi)], [23(phi),24(psi)], [24(phi),25(psi)], 
		# [25(phi),26(psi)], [26(phi),27(psi)], [27(phi),28(psi)],
		# [28(phi),29(psi)]
		# where 22(phi) and 29(psi) are NOT rotated
		#
		#txtOutput.insert(END,"first = "+str(first)+" last = "+str(last)+"\n")
		phi=[] ; psi=[]
		for j in xrange(totres):
			thisres=first+j
			residue_offset = thisres - first_resid
			# first == 1-based count up
			# res 23 == index 0,1,2,...,22 for array
			# so, thisres-1 == THIS RESIDUE
			# so, thisres-2 == THE PREVIOUS RESIDUE
			# so, thisres == THE NEXT RESIDUE
			# PHI for THIS RESIDUE [previous,this] === [thisres-2,thisres-1]
			# PSI for THIS RESIDUE [this,next] === [thisres-1,thisres]
			#if(thisres==1):
			if(thisres==first_resid):
				resp_phi="NULLN"
				rest_phi=resname[residue_offset] 
				rest_psi=resname[residue_offset]
				resn_psi=resname[residue_offset+1]
			elif((thisres-1)==last_resid-1):
				resp_phi=resname[residue_offset-1]
				rest_phi=resname[residue_offset] 
				rest_psi=resname[residue_offset]
				resn_psi="NULLC"
			else:
				resp_phi=resname[residue_offset-1]  
				rest_phi=resname[residue_offset] 
				rest_psi=resname[residue_offset]
				resn_psi=resname[residue_offset+1]

			type_phi=getphi(resp_phi,rest_phi)
			type_psi=getpsi(rest_psi,resn_psi)
			
			phi_k,phi_n,phi_d,phi2_k,phi2_n,phi2_d=getparamphi(type_phi)	
			psi_k,psi_n,psi_d,psi2_k,psi2_n,psi2_d=getparampsi(type_psi)
			
			phi.append([phi_k,phi_n,phi_d,phi2_k,phi2_n,phi2_d])
			psi.append([psi_k,psi_n,psi_d,psi2_k,psi2_n,psi2_d])
	
#			testoutput.write("%i: %s\t%s\t%s\t%s %i\t%s %i\n" % (thisres,resname[thisres-2],resname[thisres-1],resname[thisres]," type_phi =",type_phi," type_psi =",type_psi))	
#		testoutput.write("\n")
		resphi.append([i,phi])	
		respsi.append([i,psi])

	# respsi is an array of parameters [(i,(x))], where x = [(k,n,d),(k,n,d),...,(k,n,d)] for each residue each range (i)

#	testoutput.close()


#	outphi=open("junk_phi","w")
#	outpsi=open("junk_psi","w")
#	for i in xrange(numranges):
#		first=rlow[i]		# i.e. 23
#		nres=rnum[i]+1
#		this_phi=resphi[i] ; this_psi=respsi[i]
#		for j in xrange(nres):
#			thisres=first+j
#			outphi.write("%i: %f\t%f\t%f\t%f\t%f\t%f\n" % (thisres,this_phi[1][j][0],this_phi[1][j][1],this_phi[1][j][2],this_phi[1][j][3],this_phi[1][j][4],this_phi[1][j][5]))
#			outpsi.write("%i: %f\t%f\t%f\t%f\t%f\t%f\n" % (thisres,this_psi[1][j][0],this_psi[1][j][1],this_psi[1][j][2],this_psi[1][j][3],this_psi[1][j][4],this_psi[1][j][5]))
#		outpsi.write("\n")
#		outphi.write("\n")
#	outphi.close()	
#	outpsi.close()	
#	print 3/0

	return 

def dihedral(coor,dihedralparam):
	dihedralenergy=0.0
	deg2rad=numpy.pi/180.0

	for i in xrange(len(dihedralparam)):
		nps=1	
		tquad=dihedralparam[i][0]	
		atm1=locale.atoi(tquad[0])	
		atm2=locale.atoi(tquad[1])	
		atm3=locale.atoi(tquad[2])	
		atm4=locale.atoi(tquad[3])	
		ldp=len(dihedralparam[i])
		coor1=coor[atm1] ; coor2=coor[atm2] ; coor3=coor[atm3] ; coor4=coor[atm4]
		xi=deg2rad*sascalc.vdihed(coor1,coor2,coor3,coor4)
		k1xi=locale.atof(dihedralparam[i][1][0])	
		n1xi=locale.atof(dihedralparam[i][1][1])	
		t1xi=deg2rad*locale.atof(dihedralparam[i][1][2])	
		if(ldp>2):
			k2xi=locale.atof(dihedralparam[i][2][0])	
			n2xi=locale.atof(dihedralparam[i][2][1])	
			t2xi=deg2rad*locale.atof(dihedralparam[i][2][2])	
			nps=2
			print 'k1xi = ',k1xi
			print 'n1xi = ',n1xi
			print 't1xi = ',t1xi
			print 3/0
		else:
			dihedralenergy=dihedralenergy+k1xi*(1.0+numpy.cos(n1xi*xi-t1xi))

	return dihedralenergy


 
