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
import os,sys,locale,string
import sassie.interface.input_filter as input_filter
import sasmol.sasmol as sasmol
import sassie.tools.contrast_calculator.contrast_helper as contrast_helper
#import contrast_helper as contrast_helper

#def fasta_check(filename):

def check_numfiles(numfiles):

    error = []

    try: 
        intvalue=int(numfiles)
    except ValueError:
        error.append("The number of input files must be an integer.")
        return error
    if (intvalue < 0):
        error.append("The number of input files must be greater than or equal to zero.")
        return error
    return error

def check_numsolvcomp(numSolvComp):

    error = []

    try: 
        intvalue=int(numSolvComp)
    except ValueError:
        error.append("The number of solvent components must be an integer.")
        return error
    if (intvalue < 0):
        error.append("The number of solvent components must be greater than or equal to zero.")
        return error	

    return error

def check_numchemcomp(numformulas):

    error = []

    try: 
        intvalue=int(numformulas)
    except ValueError:
        error.append("The number of additional components must be an integer.")
        return error
    if (intvalue < 0):
        error.append("The number of additional components must be greater than or equal to zero.")
        return error	

    return error
		

def check_ivariables_exist(ivariables,glbls):
    error = []
    if ivariables not in glbls:
        error.append('Please enter at least one input file.')
    return error

def check_ivariables(inpath,ivariables):

    error = []

    allfilenames=[]
    allnumunits=[]
    allfracdeut=[]
    allmoltype=[]
    allIsFasta=[]

    for i in range(len(ivariables)):
        allfilenames.append(ivariables[i][0])
        allnumunits.append(ivariables[i][1])
        allfracdeut.append(ivariables[i][2])
        allmoltype.append(ivariables[i][3])
        allIsFasta.append(ivariables[i][4])

        ev,rv,wv=input_filter.check_permissions(inpath)
        if(not ev or not rv or not wv):
            error.append('Permission error in input file path '+inpath+':  [code = '+str(ev)+str(rv)+str(wv)+']')
            if(ev==False):
                error.append('Path does not exist.')
            elif(rv==False):
                error.append('Read permission not allowed.')
            elif(wv==False):
                error.append('Write permission not allowed.')
            return error


    for i in range(len(allfilenames)):

        filename=""
        filename=allfilenames[i].strip()


        if (not filename):
            error.append('No input file has been specified on line '+str(i+1)+'.')
            return error

        pdbfile=os.path.join(inpath,filename)
        print 'pdbfile: ', pdbfile
        ev,value=input_filter.check_pdb_dcd(pdbfile,'pdb')

        if(ev == 0):
            error.append('Input pdb file "'+pdbfile+'" does not exist.')
            return error

        if allIsFasta[i] != '1':
            if(ev==1 and value ==0):
                error.append('Input file "'+pdbfile+'" is not in valid PDB format. FASTA?')
                return error

            m1=sasmol.SasMol(0)
            m1.read_pdb(pdbfile)

            resname=m1.resnames()

#           necessary catch, with 3-letter codes, for the ambiguous DNA/RNA case
            print sorted(resname)
            if sorted(resname) == ['ADE','CYT','GUA']:
#           if ambiguity, trust the user
                moltypes=[allmoltype[i]]
            else:
                moltypes=list(set(m1.moltype()))		

#            print 'moltypes: ',moltypes
            possibilities=['dna','rna','protein']
            ithtype = allmoltype[i]
#            print 'ithtype: ',ithtype

            if ithtype not in possibilities:
                error.append("Input file "+str(i+1)+": Select a molecule type")
                return error


            matching = [s for s in moltypes if s in possibilities]
#            print 'matching: ', matching
#            print 'len(matching: ',len(matching)

            if len(matching) == 0:
                error.append("Input file "+str(i+1)+" does not contain a valid dna, rna or protein molecule.")
                return error
#            elif len(matching) > 1:  
#                error.append("Input file "+str(i+1)+" must contain a single molecule type (dna, rna, protein).")
#                return error
#            elif (matching[0] != ithtype):
#                error.append("The molecule represented by input file "+str(i+1)+" must match the molecule type (dna, rna, protein) indicated.")
#                return error

#           necessary to catch DNA/RNA files that are classified as 'protein', protein files classified as 'dna' or 'rna',
#           and RNA files that are classified as 'dna' BUT
#           DNA can be classified as 'rna' here and get through because only THY is classified as 'dna' in SasMol
            if ithtype not in matching:
                    error.append('Input file '+str(i+1)+' contains residues that do not match the moltype "'+str(ithtype)+'"')
                    return error

        else:
            fastaseq = contrast_helper.FASTA_sequence(pdbfile)
            if not fastaseq:
                error.append('The FASTA file "'+pdbfile+'" appears to be empty.')
                return error

            if allmoltype[i] == 'dna':
                fastadict = contrast_helper.dna_sl()
            elif allmoltype[i] == 'rna':
                fastadict = contrast_helper.rna_sl()
            elif allmoltype[i] == 'protein':
                fastadict = contrast_helper.protein_sl()
            else:
                error.append('The molecule type for input file '+str(i+1)+' must be DNA, RNA or protein.')
                return error


#           invalid letter codes that aren't protein, rna or dna are discarded in contrast_helper.FASTA_sequence above
#           so they won't be caught here and the user likely won't know that they have been ignored
#           residue letter codes are tested against given molecule type BUT
#           User beware: DNA can be called 'protein' here and get through.
            for k in fastaseq:
                if k not in fastadict:
                    error.append("Input file "+str(i+1)+" contains the non-"+allmoltype[i]+" residue "+k+".")
                    return error 

    for i in range(len(allnumunits)):
        try:
            intvalue = int(allnumunits[i])
        except ValueError:
            error.append('The number of units of each molecule for input file '+str(i+1)+' must be an integer.')
            return error
        else:
            if intvalue < 1:
                error.append('The number of units of each molecule for input file '+str(i+1)+' must be at least 1.')
                return error


    for i in range(len(allfracdeut)):
        try:
            floatvalue = float(allfracdeut[i])
        except ValueError:
            error.append('The deuteration fraction of input file '+str(i+1)+' must be a number.')
            return error

        if floatvalue < 0 or floatvalue > 1:
            error.append('The deuteration fraction of input file '+str(i+1)+' must be between 0.0 and 1.0')
            return error

    return error



def check_chemvariables(chemical_variables):

#    print 'chemvariables = ', chemical_variables

    error = []

#   chemical formula is parsed and checked in input filter
#    chemcomps=[]

    hexchem=[]
    chemdens=[]
    fexchchem=[]

    for i in range(len(chemical_variables)):
#        chemcomps.append(chemical_variables[i][0])
        hexchem.append(chemical_variables[i][1])
        fexchchem.append(chemical_variables[i][2])
        chemdens.append(chemical_variables[i][3])
        
    for i in range(len(hexchem)):
        try:
            intvalue = int(hexchem[i])
        except ValueError:
            error.append('The number of exchangeable hydrogens for additional component '+str(i+1)+' must be an integer.')
            return error
        else:
            if intvalue < 0:
                error.append('The number of exchangeable hydrogens for additional component '+str(i+1)+' must be greater than or equal to zero.')
                return error

    for i in range(len(fexchchem)):
        try:
            floatvalue = float(fexchchem[i])
        except ValueError:
            error.append('The fraction of exchangeable hydrogens for additional component '+str(i+1)+' must be a number.')
            return error
        else:
            if (floatvalue < 0.0 or floatvalue > 1.0):
                error.append('The fraction of exchangeable hydrogens for additional component '+str(i+1)+' must be between 0.0 and 1.0.')
                return error

    for i in range(len(chemdens)):
        try:
            floatvalue = float(chemdens[i])
        except ValueError:
            error.append('The mass density of additional component '+str(i+1)+' must be a number.')
            return error
        else:
            if (floatvalue <= 0.0):
                error.append('The mass density of additional component '+str(i+1)+' must be greater than 0.')
                return error

    return error



def check_solvvariables(solvvariables):

#    print 'solvvariables = ', solvvariables

    error = []

#   solvent composition is checked and parsed in input filter
#    solvcomps=[]
    solvconcs=[]

    for i in range(len(solvvariables)):
#        solvcomps.append(solvvariables[i][0])
        solvconcs.append(solvvariables[i][1])

    for i in range(len(solvconcs)):
        try:
            floatvalue = float(solvconcs[i])
        except ValueError:
            error.append('The concentration of solvent component '+str(i+1)+' must be a number.')
            return error
        else:
            if floatvalue <= 0.0:
                error.append('The concentration of solvent component '+str(i+1)+' must be greater than 0.')
                return error
    return error


def check_contrast(variables):

#    print 'variables = ',variables

    error=[]

    runname=variables['runname'][0]
    inpath=variables['inpath'][0]
    outfile=variables['outfile'][0]
    solute_conc=variables['solute_conc'][0]
    d2ostep=variables['d2ostep'][0]
    fexchp=variables['fexchp'][0]
    fexchn=variables['fexchn'][0]

    numfiles = variables['numfiles'][0]
    numformulas = variables['number_of_chemicals'][0]
    numsolv = variables['numsolv'][0]


#   checked at time of input when pressing Then Click Here
#        if(numsolv < 0):	
#            error.append("Number of non-solvent components must be greater or equal to zero")
#            return error

    if(numfiles + numformulas < 1): 
        error.append("You must enter either one file (PDB or Fasta) or one chemical formula")
        return error

    if (not runname) or (runname == "Enter project name"):
        error.append("Please enter a project name")
        return error

    if (not outfile) or (outfile == "Enter output filename"):
        error.append("Please enter an output filename")
        return error

#   checked in input filter
#    try:
#        solute_conc_float = float(solute_conc)
#    except ValueError:
#        error.append("The solute concentration must be a number.")
#    return error

    if solute_conc <= 0:
        error.append('Solute concentration must be greater than 0.')
        return error

#   checked in input filter
#    try:
#        d2ostep_int = int(d2ostep)
#    except ValueError:
#        error.append('D2O step must be an integer.')
#        return error

    if (d2ostep > 100 or d2ostep <=0):
        error.append('D2O step must be an integer between 1 and 100')
        return error
		
    if(100%d2ostep !=0):
        error.append('D2O step size must divide 100% evenly (1%, 2%, 5%, etc.)')
        return error

#   checked in input filter
#    try:
#        fexchp_float = float(fexchp)
#    except ValueError:
#        error.append("The fraction of exchangeable protein hydrogens must be a number.")
#        return error

    if (fexchp < 0.0 or fexchp > 1.0):
        error.append('The fraction of exchangeable protein hydrogens must be between 0.0 and 1.0.')
        return error

#   checked in input filter
#    try:
#        fexchn_float = float(fexchn)
#    except ValueError:
#        error.append("The fraction of exchangeable nucleic acid hydrogens must be a number.")
#        return error

    if (fexchn < 0.0 or fexchn > 1.0):
        error.append('The fraction of exchangeable nucleic acid hydrogens must be between 0.0 and 1.0.')
        return error
		
    return error

