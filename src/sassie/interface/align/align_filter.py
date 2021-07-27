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
import os,sys
import sassie.interface.input_filter as input_filter

def check_align(variables,**kwargs):


    path=variables['path'][0]
    infile=variables['infile'][0]
    pdbmol1=variables['pdbmol1'][0]
    pdbmol2=variables['pdbmol2'][0]
    basis1=variables['basis1'][0]
    basis2=variables['basis2'][0]
    lowres1=variables['lowres1'][0]
    lowres2=variables['lowres2'][0]
    highres1=variables['highres1'][0]
    highres2=variables['highres2'][0]
    ebasis1=variables['ebasis1'][0]
    ebasis2=variables['ebasis2'][0]

    error=[]

    if 'no_file_check' not in kwargs:

        ev,rv,wv=input_filter.check_permissions(path)

        infile=path+'/'+infile
        pdbmol1=path+'/'+pdbmol1
        pdbmol2=path+'/'+pdbmol2


        if(not ev or not rv or not wv):
            error.append('permission error in input file path '+path+'  [code = '+str(ev)+str(rv)+str(wv)+']')
            if(ev==False):
                error.append('path does not exist')
            elif(rv==False):
                error.append('read permission not allowed')
            elif(wv==False):
                error.append('write permission not allowed')
            return error

    print 'checking pdbmol1 : pdb'

    ev,value=input_filter.check_pdb_dcd(pdbmol1,'pdb')

    if(ev == 0):
        error.append('mol 1 reference pdb file, '+pdbmol1+', does not exist')
        return error
    if(value == 0):
        error.append( 'mol 1 reference pdb file, '+pdbmol1+', is not a valid pdb file')
        return error

    print 'checking pdbmol2 : pdb'
        
    ev,value=input_filter.check_pdb_dcd(pdbmol2,'pdb')

    if(ev == 0):
        error.append('mol 2 reference pdb file, '+pdbmol2+', does not exist')
        return error
    if(value == 0):
        error.append( 'mol 2 reference pdb file, '+pdbmol2+', is not a valid pdb file')
        return error

    ### now check and see if input file is a correct pdb or dcd     

    print 'checking infile : as pdb'
		
    if(infile[-3:] == 'pdb'):
        ev,value=input_filter.check_pdb_dcd(infile,'pdb')
        print 'pdb: ev = ',ev,'\tvalue = ',value
        value = 1	
    elif(infile[-3:] == 'dcd'):
        ev,value=input_filter.check_pdb_dcd(infile,'dcd')
        print 'dcd: ev = ',ev,'\tvalue = ',value
        value = 0	
    else:
        error.append('infile needs to have a "pdb" or "dcd" suffix')
        return error

    if(ev == 1):   # if the file exists
        if(value == 0):         # not a pdb file
            print  'checking infile : as dcd'

            ev,value=input_filter.check_pdb_dcd(infile,'dcd')
            if(value == 1):
                cvalue=input_filter.certify_pdb_dcd(pdbmol2,infile)
                if(cvalue == 0):
                    error.append('mol 2 pdbfile '+pdbmol2+' and dcdfile '+infile+' are not compatible (different number of atoms)')
                    return error
            else:
                error.append('mol 2 input file, '+infile+', is not a valid pdb or dcd file')
                return error
        else:                   # is a pdb file
            locvariables=['name']
            value,result1=input_filter.get_pdb_stats(infile,locvariables)
            value,result2=input_filter.get_pdb_stats(pdbmol2,locvariables)
            if(result1 != result2):
                error.append('mol 2 reference pdb file '+pdbmol2+' and input pdb file '+infile+' are not compatible')

                return error
    else:  # file does not exist
        error.append('input (pdb or dcd) file '+infile+' does not exist')
        return error
	
    locvariables=['resid']
    value,result1=input_filter.get_pdb_stats(pdbmol1,locvariables)

    resid1=map(int,result1[0])
    value,result2=input_filter.get_pdb_stats(pdbmol2,locvariables)
    resid2=map(int,result2[0])

    if(lowres1 not in resid1):
        error.append('mol 1 pdb file does not have low residue amino acid, '+str(lowres1)+', range = '+str(resid1[0])+' : '+str(resid1[-1]))
        return error
    elif(highres1 not in resid1):
        error.append('mol 1 pdb file does not have high residue amino acid, '+str(highres1)+', range = '+str(resid1[0])+' : '+str(resid1[-1]))
        return error
    elif(lowres2 not in resid2):
        error.append('mol 2 pdb file does not have low residue amino acid, '+str(lowres2)+', range = '+str(resid2[0])+' : '+str(resid2[-1]))
        return error
    elif(highres2 not in resid2):
        error.append('mol 2 pdb file does not have high residue amino acid, '+str(highres2)+', range = '+str(resid2[0])+' : '+str(resid2[-1]))
        return error
	
    elif(highres1 - lowres1 < 2):
        error.append('mol 1 alignment basis is too small (less than 3 points) or low residue > high residue')
        return error
    elif(highres2 - lowres2 < 2):
        error.append('mol 2 alignment basis is too small (less than 3 points) or low residue > high residue')
        return error

    if(basis1.strip().upper()==basis2.strip().upper()):	
        num1=0 ; num2=0
        for val in resid1:
            if val == highres1 :
                num1+=1

        for val in resid2:
            if val == highres2 :
                num2+=1
        if(num1 != num2):
            st='mol 1'
            if(num2>num1):	
                st='mol 2'
            #	error.append('mol 1 and mol2 have a different number of basis atoms: multiple high residues in '+st)
            #	return error
	
        diff1=highres1-lowres1
        diff2=highres2-lowres2
        if(diff1 != diff2):
            error.append('mol 1 and mol2 have a different number of basis atoms, check low and high residue input values')
            return error

    else:
        error.append('basis1 = '+basis1+' and basis2 = '+basis2+' do not match')
        return error

    if(ebasis1 == "None" and ebasis2 == "None"):
        return error
    elif(ebasis1 == "" or ebasis1 == " "):
        error.append("extra basis for molecule 1 needs to be either 'None' or an appropriate command")
        return error
    elif(ebasis2 == "" or ebasis2 == " "):
        error.append("extra basis for molecule 2 needs to be either 'None' or an appropriate command")
        return error

		# note that we do not currently have a test to confirm validity for extra basis
	
    return error


