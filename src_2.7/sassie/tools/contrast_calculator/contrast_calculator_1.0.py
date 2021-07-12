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
import os,sys,string,locale,time,math,platform
import numpy as np
import sassie.sasmol.sasmol as sasmol
from scipy import stats

#-------------------
#for formula parser
from re import findall
#-------------------
import Gnuplot,Gnuplot.PlotItems, Gnuplot.funcutils  #for plotting

#       CONTRAST
#
#	11/21/2012	--	broadly adapted from center.py :	jc/ks
#	1/3/2013        --	ready for beta testing; no FASTA
#				or solvent functionality yet		ks
#	4/1/2013	--	added plotting capability		sk
#	5/1/2013	--	works for any number of protein
#				dna or rna components with different
#				amounts of deuteration			sk
#	5/15/2013	--	solvent functionality added 
#				code follows that of Whitten et al.
#				J. Appl. Cryst. (2008). 41, 222-226	sk
#       6/18/2014       --      now works for any chemical formula.
#                               Inputs are the chemical formula,
#                               mass density, number of exchangeable
#                               hydrogens and fraction of exchangeable
#                               hydrogens that exchange                 sk 


'''
        CONTRAST is the module that calculates SLD, I(0) and contrast v. %D2O 
	for a given molecule or complex based only on sequence.  It will also
	read the sequence from a PDB file, if available.

        This module is called from Contrast Calculator from the main
        GUI through the graphical_contrast.py script.

'''

# for plotting
def wait(str=None, prompt='Plot will clear in 10 seconds ...\n'):
	if str is not None:
		print str

	try:
		if(platform.system() == "Linux"):
			import curses
			stdscr = curses.initscr()
			stdscr.addstr('press a key to continue')
			c = stdscr.getch()
			curses.endwin()
		#	raw_input(prompt)		
	except:
		time.sleep(1)


def print_failure(message,txtOutput):

	txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
	txtOutput.put(message)

	return

def unpack_variables(variables):

	runname=variables['runname'][0]
	inpath=variables['inpath'][0]
	outfile=variables['outfile'][0]
	solute_conc=variables['solute_conc'][0]
	d2ostep=variables['d2ostep'][0]
	fexchp=variables['fexchp'][0]
	fexchn=variables['fexchn'][0]
	plotflag=variables['plotflag'][0] 	#for plotting

	return runname,inpath,outfile,solute_conc,d2ostep,fexchp,fexchn,plotflag

def unpack_ivariables(ivariables):

	seqfiles=[]
	numunits=[]
	fracdeut=[]
	moltype=[]
	isFasta=[]

	numfiles = len(ivariables)
#	print 'numfiles = ', numfiles

	for i in range(numfiles):
		seqfiles.append(ivariables[i][0])
		numunits.append(ivariables[i][1])
		fracdeut.append(ivariables[i][2])
		moltype.append(ivariables[i][3])
		isFasta.append(ivariables[i][4])

	return seqfiles,numunits,fracdeut,moltype,isFasta,numfiles

def unpack_solvvariables(solvvariables):

	solv_comp=[]		#chemical formula (see below)
	solv_conc=[]		#concentration (M)

	numsolv = len(solvvariables)
#	print 'numsolv = ', numsolv
#	for i in range(numsolv):
#		print solvvariables[i]

	for i in range(numsolv):
	    if solvvariables[i][0] == "":
	    	chemical_string = ''
		numsolv = 0
	    else:		    
		chemical_string = ''
		for key,value in solvvariables[i][0].items():
	    		chemical_string += key+str(value)
#	    print 'solvent_string = ',chemical_string
	    solv_comp.append(chemical_string)
#	    solv_comp.append(solvvariables[i][0])
            solv_conc.append(solvvariables[i][1])

	return solv_comp, solv_conc, numsolv

def unpack_chemvariables(chemical_variables):

    	chem_comp=[]       	 #chemical formula (see below)
    	hexchem=[]          	 #number of exchangeable hydrogens
    	chem_density=[]    	 #mass density
    	fexchchem=[]		 #fraction of exchangeable hydrogens that actually do exchange

    	numformulas = len(chemical_variables)
#	print 'numformulas = ', numformulas
#	for i in range(numformulas):
#		print chemical_variables[i]

    	for i in range(numformulas):
	    if chemical_variables[i][0] == "":
	    	chemical_string = ''
		numformulas = 0
	    else:		
		chemical_string = ''
		for key,value in chemical_variables[i][0].items():
			chemical_string += key+str(value)
#	    print 'chemical_string = ',chemical_string
	    chem_comp.append(chemical_string)
            hexchem.append(chemical_variables[i][1])
	    fexchchem.append(chemical_variables[i][2])
            chem_density.append(chemical_variables[i][3])

    	return chem_comp,hexchem,fexchchem,chem_density,numformulas


# this should, eventually, be a sasmol getter and live elsewhere

def sequence(residlist,resnamelist):
	marker = 0
	seq = []
	for i in range(len(residlist)):
		if marker != residlist[i]:
			seq.append(resnamelist[i])
			marker = residlist[i]
	return seq

def FASTA_sequence(infile):
	seq = []
	o = open(infile)
	for line in o:
        	if line[0]!=">":
			for j in line:
				j = j.upper()
				if j in dna_sl() or j in rna_sl() or j in protein_sl():
					seq.append(j)
	return seq


# Scattering properties for elements, proteins and nucleotides; eventually needed in sasproperties!
# Scattering lengths are in units of 10^-12 cm

def element_sl():

#element name : [MW, vol Ang^3, eSL Ang, nSL Ang]
#approx volumes are per Whitten 2008

	atomic_scattering = {
        'D'   : [2.014, 5.15, 0.282, 0.6671],
        'H'   : [1.008, 5.15, 0.282, -0.3741],
        'C'   : [12.01, 16.44, 1.692, 0.6646],
        'N'   : [14.01, 2.49, 1.974, 0.9360],
        'O'   : [16.00, 9.130, 2.256, 0.5803],
        'Na'  : [22.99, 4.45, 3.102, 0.3630],
        'Mg'  : [24.31, 1.56, 3.384, 0.5375],
        'K'   : [39.10, 11.01, 5.358, 0.3670],
        'Ca'  : [40.08, 4.19, 5.640, 0.4700],
        'Cl'  : [35.45, 24.84, 4.794, 0.9577],
        'Br'  : [79.90, 31.54, 9.870, 0.6795],
        'I'   : [126.9, 44.6, 14.946, 0.5280],
        'P'   : [30.97, 3.37, 4.230, 0.5130],
        'S'   : [32.07, 26.09, 4.512, 0.2847],
        'Fe'  : [55.85, 7.99, 7.332, 0.9450],
        'Co'  : [58.93, 7.99, 7.614, 0.2490],
        'Ni'  : [58.69, 8.18, 7.896, 1.0300],
        'Cu'  : [63.55, 8.78, 8.178, 0.7718],
        'Zn'  : [65.39, 9.85, 8.460, 0.5680],
        'Au'  : [196.97, 14.75, 22.278, 0.7630]}


	return atomic_scattering

def nucleotide_sl():

# nucleotide name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang, #exchngH, #totalH]
# vol = MW*(0.586 cm^3/g)/NA

# electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
# molecule formulae are for nucleotides as embedded in a chain, per Jacrot 1976

	residue_scattering = {
		'DA'  : [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
		'DT'  : [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12], #303
		'DG'  : [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11], #328
		'DC'  : [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11], #288
		'U'   : [305.0, 296.9, 44.2, 9.28,  11.36, 2, 10],
		'A'   : [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
		'G'   : [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
		'C'   : [304.0, 295.9, 44.2, 8.68,  10.77, 2, 11]}
		
	return residue_scattering

def dna_sl():

# nucleotide name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang, #exchngH, #totalH]
# vol = MW*(0.56 cm^3/g)/NA, where partial specific volume of 0.56 is per Hearst 1962, JMB 4, 415-417
# psv depends on pH and salt, so there is a range between ~0.55-0.59

# electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
# molecule formulae are for nucleotides as embedded in a chain, per Jacrot 1976

	residue_scattering = {
		'DA'  : [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
		'DT'  : [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12], #303
		'DG'  : [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11], #328
		'DC'  : [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11],
		'ADE'  : [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
		'THY'  : [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12], #303
		'GUA'  : [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11], #328
		'CYT'  : [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11],
		'A'  : [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
		'T'  : [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12], #303
		'G'  : [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11], #328
		'C'  : [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11]} #288

	return residue_scattering

def rna_sl():

# nucleotide name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang, #exchngH, #totalH]
# vol = MW*(0.55 cm^3/g)/NA, where the partial specific volume of 0.55 is from Chien, et al. 2004,
# Biochemistry 43, 1950-1962.  psv depends on pH and salt and whether RNA is ss or ds.  This value
# is at the high end of the range for ssRNA, ~0.47-0.55.  psv is larger for dsRNA, ~0.57.

# electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
# molecule formulae are for nucleotides as embedded in a chain, per Jacrot 1976

	residue_scattering = {
		'U'   : [305.0, 296.9, 44.2, 9.28,  11.36, 2, 10],
		'A'   : [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
		'G'   : [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
		'C'   : [304.0, 295.9, 44.2, 8.68,  10.77, 2, 11],
		'URA'   : [305.0, 296.9, 44.2, 9.28,  11.36, 2, 10],
		'ADE'   : [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
		'GUA'   : [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
		'CYT'   : [304.0, 295.9, 44.2, 8.68,  10.77, 2, 11]}

	return residue_scattering


def protein_sl():

# aa residue name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang, #exchngH, #totalH]

#electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
#molecule formulae are for residues embedded in a chain, per Jacrot 1976

	residue_scattering = {
		'ALA' : [71.1,  88.6,  10.7, 1.645, 2.686, 1, 5],
		'ARG' : [156.2, 173.4, 24.0, 3.466, 9.714, 6, 13],
		'ASP' : [115.1, 111.1, 16.7, 3.845, 4.886, 1, 4],
		'ASN' : [114.1, 114.1, 16.9, 3.456, 6.580, 3, 6],
		'CYS' : [103.1, 108.5, 15.2, 1.930, 4.013, 1, 5],
		'GLU' : [129.1, 138.4, 18.9, 3.762, 4.803, 1, 6],
		'GLN' : [128.1, 143.8, 19.2, 3.373, 6.497, 3, 8],
		'GLY' : [57.1,   60.1,  8.5, 1.728, 2.769, 1, 3],
		'HSD' : [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
		'HIS' : [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
		'HSE' : [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
		'HSP' : [138.2, 153.2, 20.2, 4.771, 6.521, 3, 7],
		'ILE' : [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
		'LEU' : [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
		'LYS' : [128.2, 168.6, 20.1, 1.586, 5.752, 4, 13],
		'MET' : [131.2, 162.9, 19.8, 1.763, 2.805, 1, 9],
		'PHE' : [147.2, 189.9, 22.0, 4.139, 5.180, 1, 9],
		'PRO' : [97.1,  112.7, 14.7, 2.227, 2.227, 0, 7],
		'SER' : [87.1,  89.0,  13.0, 2.225, 4.308, 2, 5],
		'THR' : [101.1, 116.1, 15.3, 2.142, 4.224, 2, 7],
		'TRP' : [186.2, 227.8, 27.7, 6.035, 8.118, 2, 10],
		'TYR' : [163.2, 193.6, 24.3, 4.719, 6.802, 2, 9],
		'VAL' : [99.1,  140.0, 15.3, 1.478, 2.520, 1, 9],
		'A' : [71.1,  88.6,  10.7, 1.645, 2.686, 1, 5],
		'R' : [156.2, 173.4, 24.0, 3.466, 9.714, 6, 13],
		'D' : [115.1, 111.1, 16.7, 3.845, 4.886, 1, 4],
		'N' : [114.1, 114.1, 16.9, 3.456, 6.580, 3, 6],
		'C' : [103.1, 108.5, 15.2, 1.930, 4.013, 1, 5],
		'E' : [129.1, 138.4, 18.9, 3.762, 4.803, 1, 6],
		'Q' : [128.1, 143.8, 19.2, 3.373, 6.497, 3, 8],
		'G' : [57.1,   60.1,  8.5, 1.728, 2.769, 1, 3],
		'H' : [138.2, 153.2, 20.2, 4.771, 6.521, 3, 7],
		'I' : [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
		'L' : [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
		'K' : [128.2, 168.6, 20.1, 1.586, 5.752, 4, 13],
		'M' : [131.2, 162.9, 19.8, 1.763, 2.805, 1, 9],
		'F' : [147.2, 189.9, 22.0, 4.139, 5.180, 1, 9],
		'P' : [97.1,  112.7, 14.7, 2.227, 2.227, 0, 7],
		'S' : [87.1,  89.0,  13.0, 2.225, 4.308, 2, 5],
		'T' : [101.1, 116.1, 15.3, 2.142, 4.224, 2, 7],
		'W' : [186.2, 227.8, 27.7, 6.035, 8.118, 2, 10],
		'Y' : [163.2, 193.6, 24.3, 4.719, 6.802, 2, 9],
		'V' : [99.1,  140.0, 15.3, 1.478, 2.520, 1, 9]}

	return residue_scattering

#main method

def contrast(variables,ivariables,solvvariables,chemical_variables,txtOutput):

	runname,inpath,outfile,solute_conc,d2ostep,fexchp,fexchn,plotflag = unpack_variables(variables)
	seqfiles,numunits,fracdeut,moltype,isFasta,numfiles = unpack_ivariables(ivariables)
	solv_comp,solv_conc,numsolv = unpack_solvvariables(solvvariables)
    	chem_comp,hexchem,fexchchem,chem_density,numformulas = unpack_chemvariables(chemical_variables)

	print 'numfiles, numsolv, numformulas :', numfiles,numsolv,numformulas

	contpath=runname+'/contrast_calculator/'
	direxist=os.path.exists(contpath)
	if(direxist==0):
		os.system('mkdir -p '+contpath)

	prefix = outfile.find('.')
	if prefix == -1:
		coreoutfile=outfile
	else:
    		coreoutfile = outfile[0:prefix]

	izerofile=coreoutfile+'_izero.txt'
	izfile=open(contpath+izerofile,'w')

	contrastfile=coreoutfile+'_contrast.txt'
	contfile=open(contpath+contrastfile,'w')

	scatlendenfile=coreoutfile+'_sld.txt'
	sldfile=open(contpath+scatlendenfile,'w')

	#ttxt=time.ctime()
	ttxt=time.asctime( time.gmtime( time.time() ) ) 
	st=''.join(['=' for x in xrange(60)])
	txtOutput.put("\n%s \n" %(st))
	txtOutput.put("DATA FROM RUN: %s \n\n" %(ttxt))

	na = 6.023e23


#	rh2o = -0.562	  10^10 cm^-2 calculated below
#	rd2o = 6.4
#	rxh2o = 9.41
#	bd = 0.6671	  10^-12 cm  read from element table
#	bh = -0.3742

	vdna = 0.56       #using an average DNA value for both DNA and RNA, since the RNA value is only slightly smaller
	vprot = 0.730

	fexchprot = fexchp  #fraction of exchangeable H that are assumed to exchange
	fexchdna = fexchn   #now obtained from graphical contrast screen as a variable (9/2013)

#	print 'fraction of exchangeables'
#	print fexchp, fexchn


	working_conc = solute_conc/1000		#convert mg/ml to g/ml


	protein_dict = protein_sl()
	nuc_dict = nucleotide_sl()
	dna_dict = dna_sl()
	rna_dict = rna_sl()
	atomic_dict = element_sl()

#water properties used to calculate bh2o,bd2o,bxh2o,rh2o,rd2o,rxh2o

	hstats = atomic_dict['H']
#	print hstats
	dstats = atomic_dict['D']
#	print dstats
	ostats = atomic_dict['O']
#	print ostats

	bh = hstats[3]				#in units of 10^-12 cm
	bd = dstats[3]
	bo = ostats[3]
	bxh = hstats[2]
	bxo = ostats[2]

	bh2o = bh*2+bo
	bd2o = bd*2+bo
	bxh2o = bxh*2+bxo
#	print bh, bd, bxh, bh2o, bd2o, bxh2o

	h2o_conc = 1.00e3/(hstats[0]*2+ostats[0])	#molar conc of water
	d2o_conc = 1.00e3/(dstats[0]*2+ostats[0]) 
	h2o_vol = 1.0e-3/na/h2o_conc*1e30		#Ang^3	(10^-24 cm^3)
	d2o_vol = 1.0e-3/na/d2o_conc*1e30

#	print h2o_conc, d2o_conc, h2o_vol, d2o_vol


#determine if there are additional components in solvent

	if numsolv == 0:
		usesolventinfo = 0
	else:
		usesolventinfo = 1

	elementmatrix=[]

#	print 'numsolv = ', numsolv

	if usesolventinfo == 1:
		if solv_conc[0] == "":			#What if nothing was entered initially?  This takes care of that possibility.
			numsolv = 0
			usesolventinfo = 0
#	print 'numsolv = ', numsolv

#parse the solvent composition to get the total number of each element
#formula can't have parenthesis at this stage

	for m in range(numsolv):
#		print 'm = ', m		
		thissolv_comp = solv_comp[m]

#		print thissolv_comp
#		print 'length = ', len(thissolv_comp)
		
		elements=[]; elementarray=[]; numberarray=[]; 

		elements = [k for k in findall(r'([A-Z][a-z]*)(\d*)',thissolv_comp) if k]
#		print elements
#		print len(elements)
	
		for i in range(len(elements)):
#			print len(elements[i])
 			for j in range(len(elements[i])):
#  				print i,j,elements[i][j]
				if j == 0:
					elementarray.append(elements[i][j])
				if j == 1:
					if elements[i][j] == "":
						numberarray.append('1')
					else:
						numberarray.append(elements[i][j])
#		print elementarray
#		print numberarray
		elementmatrix.append([elementarray,numberarray])

#	print elementmatrix


#arrays for solvent components

	volsolvarray=[]; bxsolvarray=[]; bhsolvarray=[]; mwsolvarray=[]; concsolvarray=[]
	rxsolvarray=[]; rhsolvarray=[]


#calculate solvent properties if solvent contains only water

	if usesolventinfo == 0:

		rxh2o = bxh2o*100/h2o_vol		#in units of 10^10 cm^-2
		rh2o = bh2o*100/h2o_vol
		rd2o = bd2o*100/h2o_vol			#H2O volume is used for both H2O and D2O to get the accepted value of 6.4 for rd2o
#		rd2o2 = bd2o*100/d2o_vol
#		print rh2o, rd2o, rd2o2, rxh2o

		rsx = rxh2o
		rsd2o = rd2o
		rsh2o = rh2o
#		print rsh2o, rsd2o, rsx

#calculate solute properties if usesolventinfo=1
#H-D exchange between water and small solute molecules in the solvent is not considered

	elif usesolventinfo == 1:
		
		rxh2o = bxh2o*100/h2o_vol		#in units of 10^10 cm^-2
		rh2o = bh2o*100/h2o_vol			#for printing out to files

		soluteconc = 0.0

		for i in range(numsolv):
		
			mwsolv = 0.0
			bhsolv = 0.0
			bxsolv = 0.0
			volsolv = 0.0

			for j in range(len(elementmatrix[i][0])):

				element=elementmatrix[i][0][j]
				number=int(elementmatrix[i][1][j])
#				print 'element, number: ', element,number
				
				estats = atomic_dict[element]

				mwsolv = mwsolv + estats[0]*number
				volsolv = volsolv + estats[1]*number
				bhsolv = bhsolv + estats[3]*number
				bxsolv = bxsolv + estats[2]*number
#				print mwsolv,volsolv,bhsolv,bxsolv
		
			concsolvarray.append(float(solv_conc[i]))
			bxsolvarray.append(bxsolv)
			bhsolvarray.append(bhsolv)
			mwsolvarray.append(mwsolv)
			volsolvarray.append(volsolv)

			rxsolvarray.append(bxsolvarray[i]*100/volsolvarray[i])			#10^10 cm^-2
			rhsolvarray.append(bhsolvarray[i]*100/volsolvarray[i])

			soluteconc = soluteconc + (concsolvarray[i]*volsolvarray[i]*na)		#total conc of solutes	

#		print concsolvarray
#		print bxsolvarray
#		print bhsolvarray
#		print mwsolvarray
#		print volsolvarray
#		print rxsolvarray
#		print rhsolvarray
#		print 'solute conc: ', soluteconc	

#		adjust water concentration due to presence of solutes

		newh2o_conc = h2o_conc*(1e27 - soluteconc)/1e27					#1e27 is the number of A^3 per liter

#		print 'new h2o conc: ', newh2o_conc
						
		bxsolvtot = bxh2o*newh2o_conc
		bhsolvtot = bh2o*newh2o_conc
		bdsolvtot = bd2o*newh2o_conc
		volsolvtot = h2o_vol*newh2o_conc

#		print bhsolvtot, bdsolvtot, bxsolvtot, volsolvtot

		for i in range(numsolv):

			bxsolvtot = bxsolvtot + concsolvarray[i]*bxsolvarray[i]
			bhsolvtot = bhsolvtot + concsolvarray[i]*bhsolvarray[i]
			bdsolvtot = bdsolvtot + concsolvarray[i]*bhsolvarray[i]
			volsolvtot = volsolvtot + concsolvarray[i]*volsolvarray[i]

#		print bh, bd, bhsolvtot, bdsolvtot, bxsolvtot, volsolvtot

		rhsolvtot = bhsolvtot*100/volsolvtot						#10^10 cm^-2
		rdsolvtot = bdsolvtot*100/volsolvtot
		rxsolvtot = bxsolvtot*100/volsolvtot		

#		print rhsolvtot, rdsolvtot, rxsolvtot
			
		rsx = rxsolvtot
		rsd2o = rdsolvtot
		rsh2o = rhsolvtot

#		print rsh2o, rsd2o, rsx


#Get the properties of the COMPLEX.

#First, determine if there are any chemical forumlas, which indicates that there are components
#of the complex that are not protein, dna or rna.


	formulamatrix=[]; allformulas=[]; vxchemarray=[]
	mwchemarray=[]; hchemarray=[]; dchemarray=[]; hexchemarray=[]; fexchchemarray=[]
	denschemarray=[]; btotchemarray=[]; bxtotchemarray=[]

#	print 'numformulas =', numformulas


#Each chemical formula is treated as a separate component regardless of the amount of deuteration
#because the mass densities could be different

	if numformulas > 0:
#		print 'CHEMICAL FORMULA'

#parse the chemical formula to get the total number of each element

		for i in range(numformulas):
			thischem_comp = chem_comp[i]
			allformulas.append(thischem_comp)

#			print thischem_comp
#			print 'length = ', len(thischem_comp)
		
			elements=[]; elementarray=[]; numberarray=[]; 

			elements = [k for k in findall(r'([A-Z][a-z]*)(\d*)',thischem_comp) if k]
#			print elements
#			print len(elements)
	
			for i in range(len(elements)):
#				print len(elements[i])
 				for j in range(len(elements[i])):
#  					print i,j,elements[i][j]
					if j == 0:
						elementarray.append(elements[i][j])
					if j == 1:
						if elements[i][j] == "":
							numberarray.append('1')
						else:
							numberarray.append(elements[i][j])
#			print elementarray
#			print numberarray
			formulamatrix.append([elementarray,numberarray])

#		print formulamatrix

#calculate Mw and total x-ray and neutron scattering lengths for each formula
#the calculations assume D is explicitely listed in the chemical formula

		for i in range(numformulas):
		
			mwchem = 0.0
			hchem = 0.0
			bh2chem = 0.0
			bxchem = 0.0
			dchem = 0.0

			for j in range(len(formulamatrix[i][0])):

				elementname=formulamatrix[i][0][j]
				elementnumber=int(formulamatrix[i][1][j])
#				print 'element, number: ', elementname,elementnumber
				
				estats = atomic_dict[elementname]

				mwchem = mwchem + estats[0]*elementnumber
				bh2chem = bh2chem + estats[3]*elementnumber
				bxchem = bxchem + estats[2]*elementnumber
				if elementname == 'H':
					hchem = hchem + elementnumber
				if elementname == 'D':
					dchem = dchem + elementnumber
#			print mwchem,bh2chem,bxchem,hchem,dchem
		
			mwchemarray.append(mwchem)
			hchemarray.append(hchem)
			dchemarray.append(dchem)
			hexchemarray.append(float(hexchem[i]))
			fexchchemarray.append(float(fexchchem[i]))
			denschemarray.append(float(chem_density[i]))
			btotchemarray.append(bh2chem)
			bxtotchemarray.append(bxchem)
	
	
		for i in range(len(mwchemarray)):
			vol = mwchemarray[i]/(denschemarray[i]*na)
			vxchemarray.append(vol)

#		print mwchemarray
#		print hchemarray
#		print dchemarray
#		print hexchemarray
#		print fexchchemarray
#		print denschemarray
#		print btotchemarray
#		print bxtotchemarray
#		print vxchemarray


#Get the properties of protein, dna and/or rna in the complex.
#First, determine how many different molecule types there are and the fract deut for each given molecule type.
#Each different fract deut is one component in the total complex.

	proteinmatrix=[]; dnamatrix=[]; rnamatrix=[]
	mwprotarray=[]; hprotarray=[]; hexprotarray=[]; btotprotarray=[]; bxtotprotarray=[]
	mwdnaarray=[]; hdnaarray=[]; hexdnaarray=[]; btotdnaarray=[]; bxtotdnaarray=[]
	mwrnaarray=[]; hrnaarray=[]; hexrnaarray=[]; btotrnaarray=[]; bxtotrnaarray=[]

#	print 'numfiles =', numfiles

	if numfiles > 0:
		if seqfiles[0] == "":			#In case 0 files was entered initially but "Then Click Here" wasn't clicked
			numfiles = 0

#	print 'numfiles =', numfiles

	for index in range(numfiles):

		thismoltype = moltype[index]
		thisnumunits = int(numunits[index])
		thisfracdeut = float(fracdeut[index])
		thisIsFasta = isFasta[index]  
		thisfilename = seqfiles[index]
#		print thisfilename,thisfracdeut,thisnumunits,thisIsFasta,thismoltype,index

		if thismoltype == "protein":
			proteinmatrix.append([thisfilename,thisfracdeut,thisnumunits,thisIsFasta])
		elif thismoltype == "dna":
			dnamatrix.append([thisfilename,thisfracdeut,thisnumunits,thisIsFasta])
		elif thismoltype == "rna":
			rnamatrix.append([thisfilename,thisfracdeut,thisnumunits,thisIsFasta])
		else:
			print("Molecule type is neither protein nor dna nor rna.")
		
# First calculate the contrast parameters for the protein components.  Sort the protein values by fract deut and bin accordingly.  Then go through each bin and calculate the contrast parameters and store them in arrays.  The size of the arrays will depend on the number of components.

	proteinfdval=[]; dnafdval=[]; rnafdval=[]; allfiles=[]; allnumunits=[]

	if len(proteinmatrix) > 0:
#		print 'PROTEIN'
#		print proteinmatrix
		newproteinmatrix=sorted(proteinmatrix, key=lambda x: x[1])   #sort by fract deut
#		print newproteinmatrix

		current=newproteinmatrix[0][1]
		proteinfdval.append(newproteinmatrix[0][1])
		numvalues=1

		for i in range(1,len(newproteinmatrix)):		#determine how many components there are
			if newproteinmatrix[i][1] != current:
				numvalues = numvalues+1
				proteinfdval.append(newproteinmatrix[i][1])
				current=newproteinmatrix[i][1]
#		print 'numvalues = ', numvalues						#number of components
#		print 'proteinfdval = ', proteinfdval					#value of fract deut  

		proteinbins=[[[] for x in range(len(newproteinmatrix))] for y in range(numvalues)]
			
		for i in range(numvalues):				#sort into bins according to fract deut
			j=0
			k=0
			while j < len(newproteinmatrix):
				if newproteinmatrix[j][1] == proteinfdval[i]:
					proteinbins[i][k] = newproteinmatrix[j]
					j=j+1
					k=k+1
				else:
					j=j+1
					k=0
#		print proteinbins

#calculate Mw, number of hydrogens, number of exchangeable hydrogens and total x-ray and neutron scattering lengths for each protein component

		for i in range(numvalues):
			mwprottot = 0.0
			hprottot = 0.0
			hexprottot = 0.0
			btotprot = 0.0
			bxtotprot = 0.0

			proteinbins[i] = [e for e in proteinbins[i] if e]	#eliminate empty spaces in bins
#			print proteinbins[i]
#			print len(proteinbins[i])

			for j in range(len(proteinbins[i])):
#read the file
				filename=proteinbins[i][j][0]
				nunits=proteinbins[i][j][2]
				allfiles.append(filename)			#for printing to output files
				allnumunits.append(nunits)
#				print 'filename: ', filename

				if proteinbins[i][j][3] != '1':			#not thisIsFasta
					m1 = sasmol.SasMol(0)
					m1.read_pdb(inpath+"/"+filename)
	
					resids = m1.resid()
					resnames = m1.resname()

					seq = sequence(resids,resnames)
				else:
					seq = FASTA_sequence(inpath+"/"+filename)

				mwprot = 0.0
				hprot = 0.0
				hexprot = 0.0
				bh2oprot = 0.0
				bxprot = 0.0
		
				for m in protein_dict.keys():
					nres = seq.count(m)
					aastats = protein_dict[m]
					mwprot = mwprot + proteinbins[i][j][2]*nres*aastats[0]	#thisnumunits
					hprot =hprot + proteinbins[i][j][2]*nres*aastats[6]
					hexprot = hexprot + proteinbins[i][j][2]*nres*aastats[5]
					bh2oprot = bh2oprot + proteinbins[i][j][2]*nres*aastats[3]
					bxprot = bxprot + proteinbins[i][j][2]*nres*aastats[2]

				mwprottot = mwprottot + mwprot + proteinbins[i][j][1]*(hprot-hexprot)
				hprottot = hprottot + hprot
				hexprottot = hexprottot + hexprot
				btotprot = btotprot + bh2oprot + (hprot-hexprot)*proteinbins[i][j][1]*(bd-bh)  #thisfracdeut
				bxtotprot = bxtotprot + bxprot

			mwprotarray.append(mwprottot)
			hprotarray.append(hprottot)
			hexprotarray.append(hexprottot)
			btotprotarray.append(btotprot)
			bxtotprotarray.append(bxtotprot)

#		print mwprotarray
#		print hexprotarray
#		print btotprotarray
#		print bxtotprotarray
	
#Now repeat for the DNA contrast parameters

	if len(dnamatrix) > 0:
#		print 'DNA'
#		print dnamatrix
		newdnamatrix=sorted(dnamatrix, key=lambda x: x[1])   #sort by fract deut
#		print newdnamatrix

		current=newdnamatrix[0][1]
		dnafdval.append(newdnamatrix[0][1])
		numdnavalues=1

		for i in range(1,len(newdnamatrix)):		#determine how many components there are
			if newdnamatrix[i][1] != current:
				numdnavalues = numdnavalues+1
				dnafdval.append(newdnamatrix[i][1])
				current=newdnamatrix[i][1]
#		print 'numdnavalues = ', numdnavalues						#number of components
#		print 'dnafdval = ', dnafdval							#value of fract deut  

		dnabins=[[[] for x in range(len(newdnamatrix))] for y in range(numdnavalues)]
			
		for i in range(numdnavalues):				#sort into bins according to fract deut
			j=0
			k=0
			while j < len(newdnamatrix):
				if newdnamatrix[j][1] == dnafdval[i]:
					dnabins[i][k] = newdnamatrix[j]
					j=j+1
					k=k+1
				else:
					j=j+1
					k=0
#		print dnabins

#calculate Mw, number of hydrogens, number of exchangeable hydrogens and total x-ray and neutron scattering lengths for each DNA component

		for i in range(numdnavalues):
			mwdnatot = 0.0
			hdnatot = 0.0
			hexdnatot = 0.0
			btotdna = 0.0
			bxtotdna = 0.0

			dnabins[i] = [e for e in dnabins[i] if e]	#eliminate empty spaces in bins
#			print dnabins[i]
#			print len(dnabins[i])
			for j in range(len(dnabins[i])):
#read the file
				filename=dnabins[i][j][0]
				nunits=dnabins[i][j][2]
				allfiles.append(filename)
				allnumunits.append(nunits)
#				print 'DNA filename', filename

				if dnabins[i][j][3] != '1':			#not thisIsFasta
					m1 = sasmol.SasMol(0)
					m1.read_pdb(inpath+"/"+filename)
	
					resids = m1.resid()
					resnames = m1.resname()

					seq = sequence(resids,resnames)
				else:
					seq = FASTA_sequence(inpath+"/"+filename)

				mwdna = 0.0
				hdna = 0.0
				hexdna = 0.0
				bh2odna = 0.0
				bxdna = 0.0
		
				for m in dna_dict.keys():
					nres = seq.count(m)
					aastats = dna_dict[m]
					mwdna = mwdna + dnabins[i][j][2]*nres*aastats[0]	#thisnumunits
					hdna =hdna + dnabins[i][j][2]*nres*aastats[6]
					hexdna = hexdna + dnabins[i][j][2]*nres*aastats[5]
					bh2odna = bh2odna + dnabins[i][j][2]*nres*aastats[3]
					bxdna = bxdna + dnabins[i][j][2]*nres*aastats[2]

				mwdnatot = mwdnatot + mwdna + dnabins[i][j][1]*(hdna-hexdna)
				hdnatot = hdnatot + hdna
				hexdnatot = hexdnatot + hexdna
				btotdna = btotdna + bh2odna + (hdna-hexdna)*dnabins[i][j][1]*(bd-bh)  #thisfracdeut
				bxtotdna = bxtotdna + bxdna

			mwdnaarray.append(mwdnatot)
			hdnaarray.append(hdnatot)
			hexdnaarray.append(hexdnatot)
			btotdnaarray.append(btotdna)
			bxtotdnaarray.append(bxtotdna)

#		print mwdnaarray
#		print hdnaarray
#		print hexdnaarray
#		print btotdnaarray
#		print bxtotdnaarray

#Now repeat for RNA contrast parameters

	if len(rnamatrix) > 0:
#		print 'RNA'
#		print rnamatrix
		newrnamatrix=sorted(rnamatrix, key=lambda x: x[1])   #sort by fract deut
#		print newrnamatrix

		current=newrnamatrix[0][1]
		rnafdval.append(newrnamatrix[0][1])
		numrnavalues=1

		for i in range(1,len(newrnamatrix)):		#determine how many components there are
			if newrnamatrix[i][1] != current:
				numrnavalues = numrnavalues+1
				rnafdval.append(newrnamatrix[i][1])
				current=newrnamatrix[i][1]
#		print 'numrnavalues = ', numrnavalues						#number of components
#		print 'rnafdval = ', rnafdval							#value of fract deut  

		rnabins=[[[] for x in range(len(newrnamatrix))] for y in range(numrnavalues)]
			
		for i in range(numrnavalues):				#sort into bins according to fract deut
			j=0
			k=0
			while j < len(newrnamatrix):
				if newrnamatrix[j][1] == rnafdval[i]:
					rnabins[i][k] = newrnamatrix[j]
					j=j+1
					k=k+1
				else:
					j=j+1
					k=0
#		print rnabins

#calculate Mw, number of hydrogens, number of exchangeable hydrogens and total x-ray and neutron scattering lengths for each RNA component

		for i in range(numrnavalues):
			mwrnatot = 0.0
			hrnatot = 0.0
			hexrnatot = 0.0
			btotrna = 0.0
			bxtotrna = 0.0

			rnabins[i] = [e for e in rnabins[i] if e]	#eliminate empty spaces in bins
#			print rnabins[i]
#			print len(rnabins[i])
			for j in range(len(rnabins[i])):
#read the file
				filename=rnabins[i][j][0]
				nunits=rnabins[i][j][2]
				allfiles.append(filename)
				allnumunits.append(nunits)
#				print 'RNA filename', filename

				if rnabins[i][j][3] != '1':			#not thisIsFasta
					m1 = sasmol.SasMol(0)
					m1.read_pdb(inpath+"/"+filename)
	
					resids = m1.resid()
					resnames = m1.resname()

					seq = sequence(resids,resnames)
				else:
					seq = FASTA_sequence(inpath+"/"+filename)

				mwrna = 0.0
				hrna = 0.0
				hexrna = 0.0
				bh2orna = 0.0
				bxrna = 0.0
		
				for m in rna_dict.keys():
					nres = seq.count(m)
					aastats = rna_dict[m]
					mwrna = mwrna + rnabins[i][j][2]*nres*aastats[0]	#thisnumunits
					hrna =hrna + rnabins[i][j][2]*nres*aastats[6]
					hexrna = hexrna + rnabins[i][j][2]*nres*aastats[5]
					bh2orna = bh2orna + rnabins[i][j][2]*nres*aastats[3]
					bxrna = bxrna + rnabins[i][j][2]*nres*aastats[2]

				mwrnatot = mwrnatot + mwrna + rnabins[i][j][1]*(hrna-hexrna)
				hrnatot = hrnatot + hrna
				hexrnatot = hexrnatot + hexrna
				btotrna = btotrna + bh2orna + (hrna-hexrna)*rnabins[i][j][1]*(bd-bh)  #thisfracdeut
				bxtotrna = bxtotrna + bxrna

			mwrnaarray.append(mwrnatot)
			hrnaarray.append(hrnatot)
			hexrnaarray.append(hexrnatot)
			btotrnaarray.append(btotrna)
			bxtotrnaarray.append(bxtotrna)

#		print mwrnaarray
#		print hrnaarray
#		print hexrnaarray
#		print btotrnaarray
#		print bxtotrnaarray


#print date and list of filenames and chemical formulas used to output files

	izfile.write('#Date: '+ttxt+'\n\n')
	contfile.write('#Date: '+ttxt+'\n\n')
	sldfile.write('#Date: '+ttxt+'\n\n')

	for i in range(len(allfiles)):
		if i == 0:
			files = str(allfiles[i])+' ('+str(allnumunits[i])+')'
		else:
			files = files+', '+str(allfiles[i])+' ('+str(allnumunits[i])+')'
	if len(allfiles) > 0:
#		print 'Files used: ', files
		izfile.write('#Files used: '+files+'\n\n')
		contfile.write('#Files used: '+files+'\n\n')
		sldfile.write('#Files used: '+files+'\n\n')

	for i in range(len(allformulas)):
		if i == 0:
			formulas = str(i+1)+'. '+str(allformulas[i])
		else:
			formulas = formulas+', '+str(i+1)+'. '+str(allformulas[i])
	if len(allformulas) > 0:
#		print 'Chemical formulas used: ', formulas
		izfile.write('#Chemical formulas used: '+formulas+'\n\n')
		contfile.write('#Chemical formulas used: '+formulas+'\n\n')
		sldfile.write('#Chemical formulas used: '+formulas+'\n\n')


#print description of the solvent components and their parameters
#Solvent component, molar conc, volume, mass, bhtot, bxtot, rxtot, rhtot

	izfile.write('#Solvent Components:\n')
	contfile.write('#Solvent Components:\n')
	sldfile.write('#Solvent Components:\n')

	izfile.write('#Component, molar conc, volume (A^3), Mw (kDA), x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
	contfile.write('#Component, molar conc, volume (A^3), Mw (kDA), x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
	sldfile.write('#Component, molar conc, volume (A^3), Mw (kDA), x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
 
#print water values

	solvstring = 'H2O'
	concstring = '%6.2f' %(h2o_conc)
	volstring = '%5.1f' %(h2o_vol)
	massstring = ' 0.018'
	bxstring = '%7.3f' %(bxh2o)
	bhstring = '%7.3f' %(bh2o)
	rxstring = '%7.3f' %(rxh2o)
	rhstring = '%7.3f' %(rh2o)

	izfile.write('#\t'+solvstring+'\t'+concstring+'\t'+volstring+'\t'+massstring+'\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')
	contfile.write('#\t'+solvstring+'\t'+concstring+'\t'+volstring+'\t'+massstring+'\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')
	sldfile.write('#\t'+solvstring+'\t'+concstring+'\t'+volstring+'\t'+massstring+'\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')


#print solute values

	if usesolventinfo == 1:

		for i in range(numsolv):
			solvstring = solv_comp[i]
			concstring = '%6.2f' %(concsolvarray[i])
			volstring = '%5.1f' %(volsolvarray[i])
			massstring = '%6.3f' %(mwsolvarray[i]/1000.0)
			bxstring = '%7.3f' %(bxsolvarray[i])
			bhstring = '%7.3f' %(bhsolvarray[i])
			rxstring = '%7.3f' %(rxsolvarray[i])
			rhstring = '%7.3f' %(rhsolvarray[i])
		
			izfile.write('#\t'+solvstring+'\t'+concstring+'\t'+volstring+'\t'+massstring+'\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')
			contfile.write('#\t'+solvstring+'\t'+concstring+'\t'+volstring+'\t'+massstring+'\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')
			sldfile.write('#\t'+solvstring+'\t'+concstring+'\t'+volstring+'\t'+massstring+'\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')

#print total SL and SLD values if there are non-water solvent components

		izfile.write('#Totals: x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
		contfile.write('#Totals: x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
		sldfile.write('#Totals: x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')		
		
		bxstring = '%7.3f' %(bxsolvtot)
		bhstring = '%7.3f' %(bhsolvtot)
		rxstring = '%7.3f' %(rsx)
		rhstring = '%7.3f' %(rsh2o)			

		izfile.write('#\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')
		contfile.write('#\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')
		sldfile.write('#\t'+bxstring+'\t'+bhstring+'\t'+rxstring+'\t'+rhstring+'\n')


#calculate x-ray properties of the system:  combine all protein components, all NA (DNA+RNA) and molecule components to find the fraction of proteins, nucleic acids and other molecules in the complex. 

	izfile.write('\n#Complex concentration:  '+str(solute_conc)+' mg/ml\n')

#Mw and x-ray scattering length

#	print 'X-RAY'
 
	mwxprot=0
	bxprot=0
	mwxdna=0
	bxdna=0
	mwxrna=0
	bxrna=0
	mwxchem=0

	fracx=[]; volx=[]; pvolx=[]; bx=[]; rx=[]; rxbar=[]
	chemxsldtitle=[]; chemxcontrasttitle=[]; chemxmwtitle=[]

	for i in range(len(mwprotarray)):
		mwxprot = mwxprot + mwprotarray[i]/1.0e3
		bxprot = bxprot + bxtotprotarray[i]
#	print mwxprot, bxprot
	
	for i in range(len(mwdnaarray)):
		mwxdna = mwxdna + mwdnaarray[i]/1.0e3
		bxdna = bxdna + bxtotdnaarray[i]
#	print mwxdna, bxdna

	for i in range(len(mwrnaarray)):
		mwxrna = mwxrna + mwrnaarray[i]/1.0e3
		bxrna = bxrna + bxtotrnaarray[i]
#	print mwxrna, bxrna

	bxtotprot = bxprot
	bxtotdna = bxdna + bxrna
	bx.append(bxtotprot)
	bx.append(bxtotdna)
	vxtotdna = vdna*(mwxdna + mwxrna)*1.0e3/na
	vxtotprot = vprot*mwxprot*1.0e3/na
	volx.append(vxtotprot)
	volx.append(vxtotdna)
	pvolx.append(vprot)
	pvolx.append(vdna)

#	print 'mwchemarray length = ', len(mwchemarray)

	for i in range(len(mwchemarray)):
		mwxchem = mwxchem + mwchemarray[i]/1.0e3
#	print mwxprot, mwxdna, mwxrna, mwxchem

	mwxcomp = mwxprot + mwxdna + mwxrna + mwxchem
#	print mwxcomp

	fxprot = mwxprot/mwxcomp
	fxdna = (mwxdna + mwxrna)/mwxcomp
	fracx.append(fxprot)
	fracx.append(fxdna)

	for i in range(len(mwchemarray)):
		fraction = mwchemarray[i]/(1.0e3*mwxcomp)
		volume = vxchemarray[i]
		partial_vol = 1.0/denschemarray[i]
		fracx.append(fraction)
		volx.append(volume)
		pvolx.append(partial_vol)
		bx.append(bxtotchemarray[i])
		chemxsldtitle.append('Molecule '+str(i+1))
		chemxcontrasttitle.append('Molecule '+str(i+1))
		chemxmwtitle.append('Molecule '+str(i+1)+' Mw')

#	print fracx
#	print volx, pvolx
#	print bx

#x-ray SLD, contrast and I(0)

	rxdna = 0.0
	rxdnabar = 0.0
	rxprot = 0.0
	rxprotbar = 0.0
	izeq = 0.0
	contrx = 0.0
	sldx = 0.0

	bxtotprot = bxtotprot*1.0e-12
	bxtotdna = bxtotdna*1.0e-12

	if (mwxdna + mwxrna) > 0:
		rxdna = (bxtotdna/vxtotdna)/1.0e10
		rxdnabar = (rxdna - rsx)
#	print rxdna,rxdnabar

	if mwxprot > 0:
		rxprot = (bxtotprot/vxtotprot)/1.0e10
		rxprotbar = (rxprot - rsx)
#	print rxprot,rxprotbar

	rx.append(rxprot)
	rx.append(rxdna)
	rxbar.append(rxprotbar)
	rxbar.append(rxdnabar)

	if mwxchem > 0:
		for i in range(2,len(mwchemarray)+2):
			b = bx[i]*1.0e-12
			v = volx[i]
			rhox = (b/v)/1.0e10
			rx.append(rhox)
			rhobarx = (rhox - rsx)
			rxbar.append(rhobarx)
	
#	print rx, rsx
#	print rxbar


	for i in range(len(mwchemarray)+2):
		izeq = izeq + fracx[i]*pvolx[i]*rxbar[i]*1.0e10
		sldx = sldx + fracx[i]*rx[i]
		contrx = contrx + fracx[i]*rxbar[i]
#	print sldx, contrx		


	izerox = working_conc*(mwxcomp)*1.0e3/na*(izeq)**2
#	print izerox


#	write information to files

	protxsldtitle = 'Protein'
	dnaxsldtitle = 'NA'
	protxcontrasttitle = 'Protein'
	dnaxcontrasttitle = 'NA'
	protxmwtitle = 'Protein Mw'
	dnaxmwtitle =  'NA Mw'

	
	if mwxprot > 0:
		xsldtitle = protxsldtitle
		xcontrasttitle = protxcontrasttitle
		xmwtitle = protxmwtitle
		if (mwxdna + mwxrna) > 0:
			xsldtitle = xsldtitle+', '+dnaxsldtitle
			xcontrasttitle = xcontrasttitle+', '+dnaxcontrasttitle
			xmwtitle = xmwtitle+', '+dnaxmwtitle
		if mwxchem > 0:
			for i in range(len(mwchemarray)):
				xsldtitle = xsldtitle+', '+chemxsldtitle[i]
				xcontrasttitle = xcontrasttitle+', '+chemxcontrasttitle[i]
				xmwtitle = xmwtitle+', '+chemxmwtitle[i]	
	elif (mwxdna + mwxrna) > 0:
		xsldtitle = dnaxsldtitle
		xcontrasttitle = dnaxcontrasttitle
		xmwtitle = dnaxmwtitle
		if mwxchem > 0:
			for i in range(len(mwchemarray)):
				xldtitle = xsldtitle+', '+chemxsldtitle[i]
				xcontrasttitle = xcontrasttitle+', '+chemxcontrasttitle[i]
				xmwtitle = xmwtitle+', '+chemxmwtitle[i]	
	elif mwxchem > 0:
		xsldtitle = chemxsldtitle[0]
		xcontrasttitle = chemxcontrasttitle[0]
		xmwtitle = chemxmwtitle[0]		
		for i in range(1,len(mwchemarray)):
			xsldtitle = xsldtitle+', '+chemxsldtitle[i]
			xcontrasttitle = xcontrasttitle+', '+chemxcontrasttitle[i]
			xmwtitle = xmwtitle+', '+chemxmwtitle[i]	

#	print xsldtitle
#	print xcontrasttitle
#	print xmwtitle
	

	izfile.write('\n#XRAY I(0):\n')
	sldfile.write('\n#XRAY SLD (10^10 cm^-2):\n')
	contfile.write('\n#XRAY Contrast (10^10 cm^-2):\n')


	izfile.write('#'+xmwtitle+', Complex Mw (kDa), I(0) (cm^-1)\n')
	sldfile.write('#'+xsldtitle+', Complex, Solvent\n')
	contfile.write('#'+xcontrasttitle+', Complex\n')

	xsldtable=[]; xcontrasttable=[]; xmwtable=[]	

	mwxna = mwxdna + mwxrna

	if mwxprot > 0:
		xsldtable.append(rx[0])
		xcontrasttable.append(rxbar[0])
		xmwtable.append(mwxprot)

		if mwxna > 0:
			xsldtable.append(rx[1])
			xcontrasttable.append(rxbar[1])
			xmwtable.append(mwxna)	
		if mwxchem > 0:
			for j in range(2,len(mwchemarray)+2):
				xsldtable.append(rx[j])
				xcontrasttable.append(rxbar[j])
				xmwtable.append(mwchemarray[j-2]/1.0e3)	
	elif mwxna > 0:
		xsldtable.append(rx[1])
		xcontrasttable.append(rxbar[1])
		xmwtable.append(mwxna)
		if mwxchem > 0:
			for j in range(2,len(mwchemarray)+2):
				xsldtable.append(rx[j])
				xcontrasttable.append(rxbar[j])
				xmwtable.append(mwchemarray[j-2]/1.0e3)	
	elif mwxchem > 0:
		for j in range(2,len(mwchemarray)+2):
			xsldtable.append(rx[j])
			xcontrasttable.append(rxbar[j])
			xmwtable.append(mwchemarray[j-2]/1.0e3)


#	print xmwtable
#	print xcontrasttable
#	print xsldtable

	xmwvalstring='\t'.join(['%8.3f' % (x) for x in xmwtable])
	xmwcompstring='%8.3f' % (mwxcomp)
	xizerostring='%7.3f' % (izerox)
	izfile.write('#'+xmwvalstring+'\t'+xmwcompstring+'\t'+xizerostring+'\n')
	xsldvalstring='\t'.join(['%7.3f' % (x) for x in xsldtable])
	xsldcompstring='%7.3f' % (sldx)
	xrsstring='%7.3f' % (rsx)
	sldfile.write('#'+xsldvalstring+'\t'+xsldcompstring+'\t'+xrsstring+'\n')
	xcontrastvalstring='\t'.join(['%7.3f' % (x) for x in xcontrasttable])
	xcontrastcompstring='%7.3f' % (contrx)
	contfile.write('#'+xcontrastvalstring+'\t'+xcontrastcompstring+'\n')
	


#calculate neutron properties of the system:  number of components for each molecule type depends on proteinfdval, dnafdval and rnafdval.  Combine DNA and RNA results for same fract deut.  (It is unlikely that there will be both DNA and RNA in a complex, but both will be picked up if they exist.)

	izfile.write('\n#NEUTRONS:\n')
	sldfile.write('\n#NEUTRON SLDs:\n')
	contfile.write('\n#NEUTRON Contrast:\n') 

#determine the number of fd2o values

#	print 'NEUTRON'

	fd2o = []; rs = []
	d2ovals = (100/float(d2ostep))+1
#	print d2ovals

	for i in range(int(d2ovals)):
		newfd2o = (i*(float(d2ostep)))/100
		fd2o.append(newfd2o)

	fd2otitle = 'frac D2O'
#	print fd2otitle
#	print fd2o
#	print len(fd2o)

#determine the solvent from the initial values calculated for the solvent + any non-water components

	for i in range(len(fd2o)):
		newrs = (rsh2o + (rsd2o-rsh2o)*fd2o[i])*1.0e10
		rs.append(newrs)
#	print rs
		

#calculate chemical formula parameters
#column headings for SLD and contrast will depend on the number of chemical formulas

	
	chemsld=[[[] for x in range(len(mwchemarray))] for y in range(len(fd2o))]		
	chemcontrast=[[[] for x in range(len(mwchemarray))] for y in range(len(fd2o))]		
	chemmw=[[[] for x in range(len(mwchemarray))] for y in range(len(fd2o))]
	chemsldtitle=[]; chemcontrasttitle=[]; chemmwtitle=[]; chemmatchtitle=[]; chemmatchpoint=[]

	for i in range(len(mwchemarray)):

#		print 'CHEMICAL FORMULA'

		working_bchem = 0.0
		working_mwchem = 0.0
		working_vtotchem = 0.0
		working_rchem = 0.0
		working_rchembar = 0.0
		temparray=[]

		chemsldtitle.append('Molecule '+str(i+1))
		chemcontrasttitle.append('Molecule '+str(i+1))
		chemmwtitle.append('Molecule '+str(i+1)+' Mw')
		chemmatchtitle.append('Molecule '+str(i+1)+' Match Point')


		for j in range(len(fd2o)):
			
			working_bchem = (btotchemarray[i] + hexchemarray[i]*fexchchemarray[i]*fd2o[j]*(bd-bh))*1.0e-12
			working_mwchem = mwchemarray[i] + fd2o[j]*hexchemarray[i]*fexchchemarray[i]
			working_vtotchem = (working_mwchem)/(denschemarray[i]*na)			
			working_rchem = working_bchem/working_vtotchem
			working_rchembar = working_rchem - rs[j]
			chemmw[j][i] = working_mwchem/1.0e3
			chemsld[j][i] = working_rchem/1.0e10
			chemcontrast[j][i] = working_rchembar/1.0e10
			temparray.append(working_rchembar)

#calculate protein match point and print to screen and files

		x=np.array(fd2o)
		y=np.array(temparray)
		slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x, y)
		x_intercept = -intercept/slope
#		print slope, intercept, x_intercept
		matchpoint=x_intercept*100
		chemmatchpoint.append(matchpoint)
		matchstring='%.2f' % (matchpoint)
		txtOutput.put(chemmatchtitle[i]+': '+matchstring+" %D2O\n\n")
		izfile.write('#'+chemmatchtitle[i]+': '+matchstring+" %D2O\n")
		sldfile.write('#'+chemmatchtitle[i]+': '+matchstring+" %D2O\n")
		contfile.write('#'+chemmatchtitle[i]+': '+matchstring+" %D2O\n")

#	print chemmwtitle
#	print chemmw
#	print chemsldtitle
#	print chemsld
#	print chemcontrasttitle
#	print chemcontrast
#	print chemmatchtitle
#	print chemmatchpoint	


#column headings for SLD and contrast will depend on fract deut
#calculate protein parameters
	
	protsld=[[[] for x in range(len(mwprotarray))] for y in range(len(fd2o))]		
	protcontrast=[[[] for x in range(len(mwprotarray))] for y in range(len(fd2o))]		
	protmw=[[[] for x in range(len(mwprotarray))] for y in range(len(fd2o))]
	protsldtitle=[]; protcontrasttitle=[]; protmwtitle=[]; protmatchtitle=[]; protmatchpoint=[]

	for i in range(len(mwprotarray)):

#		print 'PROTEIN'

		working_bprot = 0.0
		working_mwprot = 0.0
		working_vtotprot = 0.0
		working_rprot = 0.0
		working_rprotbar = 0.0
		temparray=[]

		percentdeut = 100*proteinfdval[i]
		if percentdeut == 0:
			protsldtitle.append('Protein')
			protcontrasttitle.append('Protein')
			protmwtitle.append('Protein Mw')
			protmatchtitle.append('Protein Match Point')
		else:
			protsldtitle.append(str(percentdeut)+' %D Protein')
			protcontrasttitle.append(str(percentdeut)+' %D Protein')
			protmwtitle.append(str(percentdeut)+' %D Protein Mw')
			protmatchtitle.append(str(percentdeut)+' %D Protein Match Point')

		for j in range(len(fd2o)):
			
			working_bprot = (btotprotarray[i] + hexprotarray[i]*fexchprot*fd2o[j]*(bd-bh))*1.0e-12
			working_mwprot = mwprotarray[i] + fd2o[j]*hexprotarray[i]*fexchprot
			working_vtotprot = (vprot*working_mwprot)/na			
			working_rprot = working_bprot/working_vtotprot
			working_rprotbar = working_rprot - rs[j]
			protmw[j][i] = working_mwprot/1.0e3
			protsld[j][i] = working_rprot/1.0e10
			protcontrast[j][i] = working_rprotbar/1.0e10
			temparray.append(working_rprotbar)

#calculate protein match point and print to screen and files

		x=np.array(fd2o)
		y=np.array(temparray)
		slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x, y)
		x_intercept = -intercept/slope
#		print slope, intercept, x_intercept
		matchpoint=x_intercept*100
		protmatchpoint.append(matchpoint)
		matchstring='%.2f' % (matchpoint)
		txtOutput.put(protmatchtitle[i]+': '+matchstring+" %D2O\n\n")
		izfile.write('#'+protmatchtitle[i]+': '+matchstring+" %D2O\n")
		sldfile.write('#'+protmatchtitle[i]+': '+matchstring+" %D2O\n")
		contfile.write('#'+protmatchtitle[i]+': '+matchstring+" %D2O\n")

#	print protmwtitle
#	print protmw
#	print protsldtitle
#	print protsld
#	print protcontrasttitle
#	print protcontrast
#	print protmatchtitle
#	print protmatchpoint	
	

#calculate dna params
	
	dnasld=[[[] for x in range(len(mwdnaarray))] for y in range(len(fd2o))]		
	dnacontrast=[[[] for x in range(len(mwdnaarray))] for y in range(len(fd2o))]		
	dnamw=[[[] for x in range(len(mwdnaarray))] for y in range(len(fd2o))]
	dnasldtitle=[]; dnacontrasttitle=[]; dnamwtitle=[]; dnamatchtitle=[]; dnamatchpoint=[] 

	for i in range(len(mwdnaarray)):

#		print 'DNA'

		working_bdna = 0.0
		working_mwdna = 0.0
		working_vtotdna = 0.0
		working_rdna = 0.0
		working_rdnabar = 0.0
		tempdnaarray=[]

		dnapercentdeut = 100*dnafdval[i]
		if dnapercentdeut == 0:
			dnasldtitle.append('DNA')
			dnacontrasttitle.append('DNA')
			dnamwtitle.append('DNA Mw')
			dnamatchtitle.append('DNA Match Point')
		else:
			dnasldtitle.append(str(dnapercentdeut)+' %D DNA')
			dnacontrasttitle.append(str(dnapercentdeut)+' %D DNA')
			dnamwtitle.append(str(dnapercentdeut)+' %D DNA Mw')
			dnamatchtitle.append(str(dnapercentdeut)+' %D DNA Match Point')

		for j in range(len(fd2o)):
			
			working_bdna = (btotdnaarray[i] + hexdnaarray[i]*fexchdna*fd2o[j]*(bd-bh))*1.0e-12
			working_mwdna = mwdnaarray[i] + fd2o[j]*hexdnaarray[i]*fexchdna
			working_vtotdna = (vdna*working_mwdna)/na			
			working_rdna = working_bdna/working_vtotdna
			working_rdnabar = working_rdna - rs[j]
			dnamw[j][i] = working_mwdna/1.0e3
			dnasld[j][i] = working_rdna/1.0e10
			dnacontrast[j][i] = working_rdnabar/1.0e10
			tempdnaarray.append(working_rdnabar)

#calculate DNA match point

		x=np.array(fd2o)
		y=np.array(tempdnaarray)
		slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x, y)
		x_intercept = -intercept/slope
#		print slope, intercept, x_intercept
		matchpoint=x_intercept*100
		dnamatchpoint.append(matchpoint)
		matchstring='%.2f' % (matchpoint)
		txtOutput.put(dnamatchtitle[i]+': '+matchstring+" %D2O\n\n")
		izfile.write('#'+dnamatchtitle[i]+': '+matchstring+" %D2O\n")
		sldfile.write('#'+dnamatchtitle[i]+': '+matchstring+" %D2O\n")
		contfile.write('#'+dnamatchtitle[i]+': '+matchstring+" %D2O\n")

#	print dnamwtitle
#	print dnamw
#	print dnasldtitle
#	print dnasld
#	print dnacontrasttitle
#	print dnacontrast
#	print dnamatchtitle
#	print dnamatchpoint		

#calculate rna params


	rnasld=[[[] for x in range(len(mwrnaarray))] for y in range(len(fd2o))]		
	rnacontrast=[[[] for x in range(len(mwrnaarray))] for y in range(len(fd2o))]		
	rnamw=[[[] for x in range(len(mwrnaarray))] for y in range(len(fd2o))]
	rnasldtitle=[]; rnacontrasttitle=[]; rnamwtitle=[]; rnamatchtitle=[]; rnamatchpoint=[] 

	for i in range(len(mwrnaarray)):

#		print 'RNA'

		working_brna = 0.0
		working_mwrna = 0.0
		working_vtotrna = 0.0
		working_rrna = 0.0
		working_rrnabar = 0.0
		temprnaarray=[]

		rnapercentdeut = 100*rnafdval[i]
		if rnapercentdeut == 0:
			rnasldtitle.append('RNA')
			rnacontrasttitle.append('RNA')
			rnamwtitle.append('RNA Mw')
			rnamatchtitle.append('RNA Match Point')
		else:
			rnasldtitle.append(str(rnapercentdeut)+' %D RNA')
			rnacontrasttitle.append(str(rnapercentdeut)+' %D RNA')
			rnamwtitle.append(str(rnapercentdeut)+' %D RNA Mw')
			rnamatchtitle.append(str(rnapercentdeut)+' %D RNA Match Point')

		for j in range(len(fd2o)):
			
			working_brna = (btotrnaarray[i] + hexrnaarray[i]*fexchdna*fd2o[j]*(bd-bh))*1.0e-12
			working_mwrna = mwrnaarray[i] + fd2o[j]*hexrnaarray[i]*fexchdna
			working_vtotrna = (vdna*working_mwrna)/na			
			working_rrna = working_brna/working_vtotrna
			working_rrnabar = working_rrna - rs[j]
			rnamw[j][i] = working_mwrna/1.0e3
			rnasld[j][i] = working_rrna/1.0e10
			rnacontrast[j][i] = working_rrnabar/1.0e10
			temprnaarray.append(working_rrnabar)

#calculate DNA match point

		x=np.array(fd2o)
		y=np.array(temprnaarray)
		slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x, y)
		x_intercept = -intercept/slope
#		print slope, intercept, x_intercept
		matchpoint=x_intercept*100
		rnamatchpoint.append(matchpoint)
		matchstring='%.2f' % (matchpoint)
		txtOutput.put(rnamatchtitle[i]+': '+matchstring+" %D2O\n\n")
		izfile.write('#'+rnamatchtitle[i]+': '+matchstring+" %D2O\n")
		sldfile.write('#'+rnamatchtitle[i]+': '+matchstring+" %D2O\n")
		contfile.write('#'+rnamatchtitle[i]+': '+matchstring+" %D2O\n")

#	print rnamwtitle
#	print rnamw
#	print rnasldtitle
#	print rnasld
#	print rnacontrasttitle
#	print rnacontrast
#	print rnamatchtitle
#	print rnamatchpoint		


#calculate complex parameters
#first calculate mwcomp and then use the result to calculate fdna, fprot, fchem, sldcomp, contrastcomp, izero as func of fd2o

	mwtotprot=[]; mwtotdna=[]; mwtotrna=[]; mwtotchem=[]; mwcomp=[]; sldcomp=[]; contrastcomp=[]

#	calculate total Mw of the complex as a func of fd2o

#	print 'COMPLEX'

	for j in range(len(fd2o)):

		working_mwtotprot = 0.0
		working_mwtotdna = 0.0
		working_mwtotrna = 0.0
		working_mwtotchem = 0.0
		working_mwcomp = 0.0

		for i in range(len(mwprotarray)):
			working_mwtotprot = working_mwtotprot + protmw[j][i]
		mwtotprot.append(working_mwtotprot)

		for i in range(len(mwdnaarray)):
			working_mwtotdna = working_mwtotdna + dnamw[j][i]
		mwtotdna.append(working_mwtotdna)

		for i in range(len(mwrnaarray)):
			working_mwtotrna = working_mwtotrna + rnamw[j][i]
		mwtotrna.append(working_mwtotrna)	

		for i in range(len(mwchemarray)):
			working_mwtotchem = working_mwtotchem + chemmw[j][i]
		mwtotchem.append(working_mwtotchem)	

		working_mwcomp = mwtotprot[j]+ mwtotdna[j] + mwtotrna[j] + mwtotchem[j]
		mwcomp.append(working_mwcomp)
	
#	print mwtotprot
#	print mwtotdna
#	print mwtotrna
#	print mwtotchem
#	print mwcomp

#calculate Mw fraction vs Mw total for each protein, DNA, RNA and molecule component

	fprot=[[[] for x in range(len(mwprotarray))] for y in range(len(fd2o))]
	fdna=[[[] for x in range(len(mwdnaarray))] for y in range(len(fd2o))]
	frna=[[[] for x in range(len(mwrnaarray))] for y in range(len(fd2o))]
	fchem=[[[] for x in range(len(mwchemarray))] for y in range(len(fd2o))]

	for j in range(len(fd2o)):

		for i in range(len(mwprotarray)):
			fprot[j][i] = protmw[j][i]/mwcomp[j]
		
		for i in range(len(mwdnaarray)):
 			fdna[j][i] = dnamw[j][i]/mwcomp[j]

		for i in range(len(mwrnaarray)):
			frna[j][i] = rnamw[j][i]/mwcomp[j]

		for i in range(len(mwchemarray)):
			fchem[j][i] = chemmw[j][i]/mwcomp[j]

#	print fprot
#	print fdna
#	print frna
#	print fchem


#first calculate parameters for protein, DNA, RNA and molecule component and then add them up for each fd2o

	sldtotprot=[]; sldtotdna=[]; sldtotrna=[]; sldtotchem=[] 
	sldcomp=[]; izerocomp=[]; sqrtizerocomp=[]; contrastcomp=[]
	contrasttotprot=[]; contrasttotdna=[]; contrasttotrna=[]; contrasttotchem=[]

	for j in range(len(fd2o)):	
		
		working_rtotprot = 0.0
		working_rbartotprot = 0.0
		working_rtotdna = 0.0
		working_rbartotdna = 0.0
		working_rtotrna = 0.0
		working_rbartotrna = 0.0
		working_rtotchem = 0.0
		working_rbartotchem = 0.0
		working_vchem = 0.0
		working_rcomp = 0.0
		working_rbarcomp = 0.0
		working_izero = 0.0
		working_sqrti = 0.0

		rs[j] = rs[j]/1.0e10		#for printing out

		for i in range(len(mwprotarray)):
			working_rtotprot = working_rtotprot + fprot[j][i]*protsld[j][i]
			working_rbartotprot = working_rbartotprot + fprot[j][i]*protcontrast[j][i]
		sldtotprot.append(working_rtotprot)
		contrasttotprot.append(working_rbartotprot)
		
		for i in range(len(mwdnaarray)):
			working_rtotdna = working_rtotdna + fdna[j][i]*dnasld[j][i]
			working_rbartotdna = working_rbartotdna + fdna[j][i]*dnacontrast[j][i]
		sldtotdna.append(working_rtotdna)
		contrasttotdna.append(working_rbartotdna)
		
		for i in range(len(mwrnaarray)):
			working_rtotrna = working_rtotrna + frna[j][i]*rnasld[j][i]
			working_rbartotrna = working_rbartotrna + frna[j][i]*rnacontrast[j][i]
		sldtotrna.append(working_rtotrna)
		contrasttotrna.append(working_rbartotrna)

		for i in range(len(mwchemarray)):
			working_rtotchem = working_rtotchem + fchem[j][i]*chemsld[j][i]
			working_rbartotchem = working_rbartotchem + fchem[j][i]*chemcontrast[j][i]
			working_vchem = 1.0/denschemarray[i]
		sldtotchem.append(working_rtotchem)
		contrasttotchem.append(working_rbartotchem)

		working_rcomp = sldtotprot[j] + sldtotdna[j] + sldtotrna[j] + sldtotchem[j]
		working_rbarcomp = contrasttotprot[j] + contrasttotdna[j] + contrasttotrna[j] + contrasttotchem[j]
		sldcomp.append(working_rcomp)
		contrastcomp.append(working_rbarcomp)

		working_izero = working_conc*mwcomp[j]*1.0e3/na*(vprot*contrasttotprot[j]*1.0e10 + vdna*(contrasttotdna[j] + contrasttotrna[j])*1.0e10 + working_vchem*contrasttotchem[j]*1.0e10)**2

		working_sqrti = math.sqrt(working_izero)

		if contrastcomp[j] < 0.0:
			working_sqrti = -working_sqrti

		izerocomp.append(working_izero)
		sqrtizerocomp.append(working_sqrti)


#	print 'SLD'
#	print sldtotprot
#	print sldtotdna
#	print sldtotrna
#	print sldtotchem
#	print sldcomp
#	print 'Contrast'
#	print contrasttotprot
#	print contrasttotdna
#	print contrasttotrna
#	print contrasttotchem	
#	print contrastcomp
#	print 'izero'
#	print izerocomp
#	print sqrtizerocomp

#calculate Complex match point

	x=np.array(fd2o)
	y=np.array(contrastcomp)
	slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x, y)
	x_intercept = -intercept/slope
#	print slope, intercept, x_intercept
	matchpoint=x_intercept*100
	matchstring='%.2f' % (matchpoint)
	txtOutput.put('Complex Match Point: '+matchstring+" %D2O\n\n")
	izfile.write('#Complex Match Point: '+matchstring+" %D2O\n")
	sldfile.write('#Complex Match Point: '+matchstring+" %D2O\n")
	contfile.write('#Complex Match Point: '+matchstring+" %D2O\n")


#	print out fraction of exchangeable protein, nucleic acid and molecule hydrogens that actually do exchange to output files

	fexchp_string='%.2f' % (fexchp)
	fexchn_string='%.2f' % (fexchn)
	izfile.write('\n#Fraction of exchanged protein hydrogens: '+fexchp_string+'\n')
	izfile.write('#Fraction of exchanged nucleic acid hydrogens: '+fexchn_string+'\n')
	sldfile.write('\n#Fraction of exchanged protein hydrogens: '+fexchp_string+'\n')
	sldfile.write('#Fraction of exchanged nucleic acid hydrogens: '+fexchn_string+'\n')
	contfile.write('\n#Fraction of exchanged protein hydrogens: '+fexchp_string+'\n')
	contfile.write('#Fraction of exchanged nucleic acid hydrogens: '+fexchn_string+'\n')

	for i in range(len(mwchemarray)):
		fexchc = fexchchemarray[i]
		fexchc_string='%.2f' % (fexchc)
		izfile.write('#Fraction of exchanged molecule '+str(i+1)+' hydrogens: '+fexchc_string+'\n')
		sldfile.write('#Fraction of exchanged molecule '+str(i+1)+' hydrogens: '+fexchc_string+'\n')
		contfile.write('#Fraction of exchanged molecule '+str(i+1)+' hydrogens: '+fexchc_string+'\n')


#	combine protein, dna, rna and molecule matrices to create contrast, sld and izero tables to print to files
	
	if len(mwprotarray) > 0:
		sldtitle = protsldtitle[0]
		contrasttitle = protcontrasttitle[0]
		mwtitle=protmwtitle[0]
		for i in range(1,len(mwprotarray)):
			sldtitle = sldtitle+', '+protsldtitle[i]
			contrasttitle = contrasttitle+', '+protcontrasttitle[i]
			mwtitle = mwtitle+', '+protmwtitle[i]
		if len(mwdnaarray) > 0:
			for i in range(len(mwdnaarray)):
				sldtitle = sldtitle+', '+dnasldtitle[i]
				contrasttitle = contrasttitle+', '+dnacontrasttitle[i]
				mwtitle = mwtitle+', '+dnamwtitle[i]
		if len(mwrnaarray) > 0:
			for i in range(len(mwrnaarray)):
				sldtitle = sldtitle+', '+rnasldtitle[i]
				contrasttitle = contrasttitle+', '+rnacontrasttitle[i]
				mwtitle = mwtitle+', '+rnamwtitle[i]
		if len(mwchemarray) > 0:
			for i in range(len(mwchemarray)):
				sldtitle = sldtitle+', '+chemsldtitle[i]
				contrasttitle = contrasttitle+', '+chemcontrasttitle[i]
				mwtitle = mwtitle+', '+chemmwtitle[i]

	elif len(mwdnaarray) > 0:
		sldtitle = dnasldtitle[0]
		contrasttitle = dnacontrasttitle[0]
		mwtitle=dnamwtitle[0]
		for i in range(1,len(mwdnaarray)):
			sldtitle = sldtitle+', '+dnasldtitle[i]
			contrasttitle = contrasttitle+', '+dnacontrasttitle[i]
			mwtitle = mwtitle+', '+dnamwtitle[i]
		if len(mwrnaarray) > 0:
			for i in range(len(mwrnaarray)):
				sldtitle = sldtitle+', '+rnasldtitle[i]
				contrasttitle = contrasttitle+', '+rnacontrasttitle[i]
				mwtitle = mwtitle+', '+rnamwtitle[i]
		if len(mwchemarray) > 0:
			for i in range(len(mwchemarray)):
				sldtitle = sldtitle+', '+chemsldtitle[i]
				contrasttitle = contrasttitle+', '+chemcontrasttitle[i]
				mwtitle = mwtitle+', '+chemmwtitle[i]

	elif len(mwrnaarray) > 0:
		sldtitle = rnasldtitle[0]
		contrasttitle = rnacontrasttitle[0]
		mwtitle=rnamwtitle[0]
		for i in range(1,len(mwrnaarray)):
			sldtitle = sldtitle+', '+rnasldtitle[i]
			contrasttitle = contrasttitle+', '+rnacontrasttitle[i]
			mwtitle = mwtitle+', '+rnamwtitle[i]
		if len(mwchemarray) > 0:
			for i in range(len(mwchemarray)):
				sldtitle = sldtitle+', '+chemsldtitle[i]
				contrasttitle = contrasttitle+', '+chemcontrasttitle[i]
				mwtitle = mwtitle+', '+chemmwtitle[i]

	elif len(mwchemarray) > 0:
		sldtitle = chemsldtitle[0]
		contrasttitle = chemcontrasttitle[0]
		mwtitle=chemmwtitle[0]
		for i in range(1,len(mwchemarray)):
			sldtitle = sldtitle+', '+chemsldtitle[i]
			contrasttitle = contrasttitle+', '+chemcontrasttitle[i]
			mwtitle = mwtitle+', '+chemmwtitle[i]

#	print sldtitle
#	print contrasttitle
#	print mwtitle
	
	sldfile.write('\n#NEUTRON SLDs (10^10 cm^-2):\n')
	contfile.write('\n#NEUTRON Contrast (10^10 cm^-2):\n') 
	izfile.write('\n# frac D2O, '+mwtitle+', Complex Mw (kDa), I(0) (cm^-1), sqrtI(0)\n')	
	sldfile.write('# frac D2O, '+sldtitle+', Complex, Solvent\n')
	contfile.write('# frac D2O, '+contrasttitle+', Complex\n')

	if len(mwprotarray) > 0:
		sldtable=protsld
		contrasttable=protcontrast
		mwtable=protmw
		for j in range(len(fd2o)):
			if len(mwdnaarray) > 0:
				sldtable[j].extend(dnasld[j])
				contrasttable[j].extend(dnacontrast[j])
				mwtable[j].extend(dnamw[j])	
			if len(mwrnaarray) > 0:
				sldtable[j].extend(rnasld[j])
				contrasttable[j].extend(rnacontrast[j])
				mwtable[j].extend(rnamw[j])
			if len(mwchemarray) > 0:
				sldtable[j].extend(chemsld[j])
				contrasttable[j].extend(chemcontrast[j])
				mwtable[j].extend(chemmw[j])		
	elif len(mwdnaarray) > 0:
		sldtable=dnasld
		contrasttable=dnacontrast
		mwtable=dnamw
		for j in range(len(fd2o)):
			if len(mwrnaarray) > 0:
				sldtable[j].extend(rnasld[j])
				contrasttable[j].extend(rnacontrast[j])
				mwtable[j].extend(rnamw[j])
			if len(mwchemarray) > 0:
				sldtable[j].extend(chemsld[j])
				contrasttable[j].extend(chemcontrast[j])
				mwtable[j].extend(chemmw[j])		
	elif len(mwrnaarray) > 0:
		sldtable=rnasld
		contrasttable=rnacontrast
		mwtable=rnamw
		for j in range(len(fd2o)):
			if len(mwchemarray) > 0:
				sldtable[j].extend(chemsld[j])
				contrasttable[j].extend(chemcontrast[j])
				mwtable[j].extend(chemmw[j])
	elif len(mwchemarray) > 0:
		sldtable=chemsld
		contrasttable=chemcontrast
		mwtable=chemmw	

	for j in range(len(fd2o)):
		fd2ostring='%5.2f' % (fd2o[j])
		mwvalstring='\t'.join(['%8.3f' % (x) for x in mwtable[j]])
		mwcompstring='%8.3f' %(mwcomp[j])
		izerostring='%7.3f' %(izerocomp[j])
		sqrtizerostring='%7.3f' %(sqrtizerocomp[j])
		izfile.write(fd2ostring+'\t'+mwvalstring+'\t'+mwcompstring+'\t'+izerostring+'\t'+sqrtizerostring+'\n')
		sldvalstring='\t'.join(['%7.3f' % (x) for x in sldtable[j]])
		sldcompstring='%7.3f' % (sldcomp[j])
		rsstring='%7.3f' % (rs[j])
		sldfile.write(fd2ostring+'\t'+sldvalstring+'\t'+sldcompstring+'\t'+rsstring+'\n')
		contrastvalstring='\t'.join(['%7.3f' % (x) for x in contrasttable[j]])
		contrastcompstring='%7.3f' % (contrastcomp[j])
		contfile.write(fd2ostring+'\t'+contrastvalstring+'\t'+contrastcompstring+'\n')
		
	txtOutput.put("\nFiles "+izerofile+", "+scatlendenfile+" and "+contrastfile+" written to ./"+contpath+".")

#plot only if there are two components or less

	if (len(mwprotarray) + len(mwdnaarray) + len(mwrnaarray) + len(mwchemarray)) > 2:
		plotflag = 0

	if(plotflag == 1):
		graph = Gnuplot.Gnuplot(debug=1)
		graph.clear()
		graph('set title "SqrtI(0) vs D2O Fraction"')
		graph('set zeroaxis')
		graph('set xtics 0,.05,1')
		graph.xlabel('D2O Fraction')
		graph.ylabel('sqrtI(0)')
		graph2 = Gnuplot.Gnuplot(debug=1)
		graph2.clear()
		graph2('set title "Scattering Length Density vs D2O Fraction"')
		graph2('set zeroaxis')
		graph2('set xtics 0,.05,1')
		graph2.xlabel('D2O Fraction')
		graph2.ylabel('SLD x 10^10 cm^-2')
		graph3 = Gnuplot.Gnuplot(debug=1)
		graph3.clear()
		graph3('set title "Contrast vs D2O Fraction"')
		graph3('set zeroaxis')
		graph3('set xtics 0,.05,1')
		graph3.xlabel('D2O Fraction')
		graph3.ylabel('Contrast x 10^10 cm^-2')
		graph4 = Gnuplot.Gnuplot(debug=1)
		graph4.clear()
		graph4('set title "I(0) vs D2O Fraction"')
		graph4('set xtics 0,.05,1')
		graph4.xlabel('D2O Fraction')
		graph4.ylabel('I(0) cm^-1')
	
#	arrays for plots
		izerocoords=[]; sizerocoords=[]; rscoords=[]; rcomp1coords=[]; rcomp2coords=[]
		ccomp1coords=[]; ccomp2coords=[]; rcomplexcoords=[]; ccomplexcoords=[]
			
		for j in range(len(fd2o)):
			izerocoords.append([fd2o[j],izerocomp[j]])
			sizerocoords.append([fd2o[j],sqrtizerocomp[j]])
			rscoords.append([fd2o[j],rs[j]])
			rcomplexcoords.append([fd2o[j],sldcomp[j]])
			ccomplexcoords.append([fd2o[j],contrastcomp[j]])

			if len(mwprotarray) == 2:
				rcomp1coords.append([fd2o[j],protsld[j][0]])
				ccomp1coords.append([fd2o[j],protcontrast[j][0]])
				comp1contrasttitle = protcontrasttitle[0]
				comp1sldtitle = protsldtitle[0]
				rcomp2coords.append([fd2o[j],protsld[j][1]])
				ccomp2coords.append([fd2o[j],protcontrast[j][1]])
				comp2contrasttitle = protcontrasttitle[1]
				comp2sldtitle = protsldtitle[1]
				continue
			elif len(mwprotarray) == 1:
				rcomp1coords.append([fd2o[j],protsld[j][0]])
				ccomp1coords.append([fd2o[j],protcontrast[j][0]])
				comp1contrasttitle = protcontrasttitle[0]
				comp1sldtitle = protsldtitle[0]		
				if len(mwdnaarray) == 1:
					rcomp2coords.append([fd2o[j],dnasld[j][0]])
					ccomp2coords.append([fd2o[j],dnacontrast[j][0]])
					comp2contrasttitle = dnacontrasttitle[0]
					comp2sldtitle = dnasldtitle[0]
				elif len(mwrnaarray) == 1:
					rcomp2coords.append([fd2o[j],rnasld[j][0]])
					ccomp2coords.append([fd2o[j],rnacontrast[j][0]])
					comp2contrasttitle = rnacontrasttitle[0]
					comp2sldtitle = rnasldtitle[0]
				elif len(mwchemarray) == 1:
					rcomp2coords.append([fd2o[j],chemsld[j][0]])
					ccomp2coords.append([fd2o[j],chemcontrast[j][0]])
					comp2contrasttitle = chemcontrasttitle[0]
					comp2sldtitle = chemsldtitle[0]
				else:
					rcomp2coords.append([fd2o[j],0.0])
					ccomp2coords.append([fd2o[j],0.0])
					comp2contrasttitle = 'None'
					comp2sldtitle = 'None'
				continue
			elif len(mwdnaarray) == 2:
				rcomp1coords.append([fd2o[j],dnasld[j][0]])
				ccomp1coords.append([fd2o[j],dnacontrast[j][0]])
				comp1contrasttitle = dnacontrasttitle[0]
				comp1sldtitle = dnasldtitle[0]
				rcomp2coords.append([fd2o[j],dnasld[j][1]])
				ccomp2coords.append([fd2o[j],dnacontrast[j][1]])
				comp2contrasttitle = dnacontrasttitle[1]
				comp2sldtitle = dnasldtitle[1]
				continue
			elif len(mwdnaarray) == 1:
				rcomp1coords.append([fd2o[j],dnasld[j][0]])
				ccomp1coords.append([fd2o[j],dnacontrast[j][0]])
				comp1contrasttitle = dnacontrasttitle[0]
				comp1sldtitle = dnasldtitle[0]		
				if len(mwrnaarray) == 1:
					rcomp2coords.append([fd2o[j],rnasld[j][0]])
					ccomp2coords.append([fd2o[j],rnacontrast[j][0]])
					comp2contrasttitle = rnacontrasttitle[0]
					comp2sldtitle = rnasldtitle[0]
				elif len(mwchemarray) == 1:
					rcomp2coords.append([fd2o[j],chemsld[j][0]])
					ccomp2coords.append([fd2o[j],chemcontrast[j][0]])
					comp2contrasttitle = chemcontrasttitle[0]
					comp2sldtitle = chemsldtitle[0]
				else:
					rcomp2coords.append([fd2o[j],0.0])
					ccomp2coords.append([fd2o[j],0.0])
					comp2contrasttitle = 'None'
					comp2sldtitle = 'None'
				continue
			elif len(mwrnaarray) == 2:
				rcomp1coords.append([fd2o[j],rnasld[j][0]])
				ccomp1coords.append([fd2o[j],rnacontrast[j][0]])
				comp1contrasttitle = rnacontrasttitle[0]
				comp1sldtitle = rnasldtitle[0]
				rcomp2coords.append([fd2o[j],rnasld[j][1]])
				ccomp2coords.append([fd2o[j],rnacontrast[j][1]])
				comp2contrasttitle = rnacontrasttitle[1]
				comp2sldtitle = rnasldtitle[1]
			elif len(mwrnaarray) == 1:
				rcomp1coords.append([fd2o[j],rnasld[j][0]])
				ccomp1coords.append([fd2o[j],rnacontrast[j][0]])
				comp1contrasttitle = rnacontrasttitle[0]
				comp1sldtitle = rnasldtitle[0]		
				if len(mwchemarray) == 1:
					rcomp2coords.append([fd2o[j],chemsld[j][0]])
					ccomp2coords.append([fd2o[j],chemcontrast[j][0]])
					comp2contrasttitle = chemcontrasttitle[0]
					comp2sldtitle = chemsldtitle[0]
				else:
					rcomp2coords.append([fd2o[j],0.0])
					ccomp2coords.append([fd2o[j],0.0])
					comp2contrasttitle = 'None'
					comp2sldtitle = 'None'
				continue
			elif len(mwchemarray) == 2:
				rcomp1coords.append([fd2o[j],chemsld[j][0]])
				ccomp1coords.append([fd2o[j],chemcontrast[j][0]])
				comp1contrasttitle = chemcontrasttitle[0]
				comp1sldtitle = chemsldtitle[0]
				rcomp2coords.append([fd2o[j],chemsld[j][1]])
				ccomp2coords.append([fd2o[j],chemcontrast[j][1]])
				comp2contrasttitle = chemcontrasttitle[1]
				comp2sldtitle = chemsldtitle[1]
				continue
			elif len(mwchemarray) == 1:
				rcomp1coords.append([fd2o[j],chemsld[j][0]])
				ccomp1coords.append([fd2o[j],chemcontrast[j][0]])
				comp1contrasttitle = chemcontrasttitle[0]
				comp1sldtitle = chemsldtitle[0]
				rcomp2coords.append([fd2o[j],0.0])
				ccomp2coords.append([fd2o[j],0.0])
				comp2contrasttitle = 'None'
				comp2sldtitle = 'None'	
		

		try:
        		graph.plot(Gnuplot.Data(sizerocoords,using='1:2 w lines lw 2',title=''))
		except:
			message='error trying to plot sqrtI(0) vs D2O Fraction data'
			message+=' :  stopping here'
			print_failure(message,txtOutput)

		try:
        		graph2.plot(Gnuplot.Data(rscoords,using='1:2 w lines lw 2',title='Solvent'),Gnuplot.Data(rcomp1coords,using='1:2 w lines lw 2',title=comp1sldtitle),Gnuplot.Data(rcomp2coords,using='1:2 w lines lw 2',title=comp2sldtitle),Gnuplot.Data(rcomplexcoords,using='1:2 w lines lw 2',title='Complex'))
		except:
			message='error trying plot Scattering Length Density vs D2O Fraction data'
			message+=' :  stopping here'
			print_failure(message,txtOutput)

		try:
        		graph3.plot(Gnuplot.Data(ccomp1coords,using='1:2 w lines lw 2',title=comp1contrasttitle),Gnuplot.Data(ccomp2coords,using='1:2 w lines lw 2',title=comp2contrasttitle),Gnuplot.Data(ccomplexcoords,using='1:2 w lines lw 2',title='Complex'))
		except:
			message='error trying plot Contrast vs D2O Fraction data'
			message+=' :  stopping here'
			print_failure(message,txtOutput)

		try:
        		graph4.plot(Gnuplot.Data(izerocoords,using='1:2 w lp pt 5 lw 2',title=''))
		except:
			message='error trying plot I(0) vs D2O Fraction data'
			message+=' :  stopping here'
			print_failure(message,txtOutput)
	
	else:
		txtOutput.put("\n\nResults are not plotted for complexes with more than 2 components.")


	st=''.join(['=' for x in xrange(60)])
	txtOutput.put('\n'+st)
	fraction_done=1
	report_string='STATUS\t'+str(fraction_done)
	txtOutput.put(report_string)

	time.sleep(2.0)

	return()


if __name__ == "__main__":

    fasta_file = "pai_seq.txt"
	
    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT

    runname = 'run_0'
    inpath = './'
    outfile = 'dum.txt'
    numfile = '2'
    #numfile = '1'
    #numfile = '2'

    solute_conc = "1.0"
    d2ostep = "5"
    numsolv = "1"
    frExchHpro = "0.95"
    frExchHnuc = "1.0"

    seqfiles = ["protein_sequence.txt", "dna_sequence.txt"] 
    #seqfiles = ["pai_seq.txt"]
    #seqfiles = ["pai_seq.txt","vn_seq.txt"]
 
    numunits = ["1", "1"]
    #numunits = ["1"]
    #numunits = ["1", "1"]
  
    fracdeut = ["0", "0"]
    #fracdeut = ["0"]
    #fracdeut = ["0.6","0.0"]
   
    moltype = ["protein", "dna"]
    #moltype = ["protein"]
    #moltype = ["protein","protein"]
    
    isFasta = ["1", "1"]
    #isFasta = ["1"]
    #isFasta = ["1", "1"]

    ivariables = []
    for i in xrange(locale.atoi(numfile)):
        ivariables.append([seqfiles[i], numunits[i], fracdeut[i], moltype[i], isFasta[i]])

    solv_comp = ["NaCl"]
    solv_conc = ["0.15"]

    number_of_chemicals = '2'
    formula_array = ["(C3H4O3)12", "(C3H4O3)12"]
    
    number_exchangeable_hydrogens = ["12", "5"]
    fraction_exchangeable_hydrogens = ["0.95", "0.45"]
    mass_density = ["1.1", "1.3"]

#    chemical_variables = []

    #for i in xrange(locale.atoi(number_of_chemicals)):
    #    this_chemical_formula = chemical_formula[i]
    #    this_number_exchangeable_hydrogens = number_exchangeable_hydrogens[i]
    #    this_fraction_exchangeable_hydrogens = fraction_exchangeable_hydrogens[i]
    #    this_mass_density = mass_density[i]
    #    chemical_variables.append([this_chemical_formula, this_number_exchangeable_hydrogens, this_fraction_exchangeable_hydrogens, this_mass_density])

    path = ''

    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT

    import sys, shutil
    import multiprocessing 
    import sassie.interface.input_filter as input_filter
    import sassie.interface.contrast_calculator_filter as contrast_calculator_filter
   
    svariables = {}

    svariables['runname'] = (str(runname), 'string')
    svariables['inpath'] = (str(inpath), 'string')
    svariables['outfile'] = (str(outfile), 'string')
    svariables['numfiles'] = (str(numfile), 'int')
    svariables['solute_conc'] = (str(solute_conc), 'float')
    svariables['d2ostep'] = (str(d2ostep), 'int')
    svariables['numsolv'] = (str(numsolv), 'int')
    svariables['fexchp'] = (str(frExchHpro), 'float')
    svariables['fexchn'] = (str(frExchHnuc), 'float')
    svariables['plotflag'] = (str('0'), 'int')
    svariables['number_of_chemicals'] = (str(number_of_chemicals), 'int')

    error = []

    error, variables = input_filter.type_check_and_convert(svariables)

    if(len(error) > 0):
        print 'error = ', error
        sys.exit()

    else:
        error = contrast_calculator_filter.check_contrast(variables)

        if(len(error) > 0):
            print 'error = ', error
            sys.exit()

        else:

            solv_variables = []

            if(int(numsolv) > 0):
                error, solv_formula = input_filter.check_and_convert_formula(solv_comp)

                if(len(error) > 0):
                    print 'error = ', str(error)
                    sys.exit() 

                for i in xrange(locale.atoi(numsolv)):
                    # solv_variables.append([solv_comp[i],solv_conc[i]])
                    solv_variables.append([solv_formula[i], solv_conc[i]])

            chemical_variables = []
            
            if(int(number_of_chemicals) > 0):
                error, formulas = input_filter.check_and_convert_formula(formula_array)
                if(len(error) > 0):
                    print 'error = ', str(error)
                    sys.exit() 
           
                for i in xrange(int(number_of_chemicals)):

                    this_chemical_formula = formulas[i]
                    this_number_exchangeable_hydrogens = number_exchangeable_hydrogens[i]
                    this_fraction_exchangeable_hydrogens = fraction_exchangeable_hydrogens[i]
                    this_mass_density = mass_density[i]
                    chemical_variables.append([this_chemical_formula, this_number_exchangeable_hydrogens, this_fraction_exchangeable_hydrogens, this_mass_density])

            runname = variables['runname'][0]

            txtQueue = multiprocessing.JoinableQueue()
            if os.path.exists(runname + '/contrast_calculator'):
                shutil.rmtree(runname + '/contrast_calculator')

            process = multiprocessing.Process(target=contrast, args=( \
                    variables, ivariables, solv_variables, chemical_variables, txtQueue))
            process.start()


    
    
    	
