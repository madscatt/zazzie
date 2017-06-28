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
import sasmol.sasmol as sasmol
import sassie.simulate.constraints.constraints as constraints
import sassie.interface.input_filter as input_filter

def check_input_values(variables, **kwargs):

    runname = variables['runname'][0]
    path = variables['path'][0]
    pdbmol1 = variables['pdbmol1'][0]
    pdbmol2 = variables['pdbmol2'][0]
    ofile = variables['ofile'][0]
    accpos = variables['accpos'][0] 
    pos = variables['pos'][0] 
    trans = variables['trans'][0]
    dtrans = variables['dtrans'][0]
    theta = variables['theta'][0]
    dtheta = variables['dtheta'][0]
    basis = variables['basis'][0]
    cutoff = variables['cutoff'][0]
    lowrg = variables['lowrg'][0]
    highrg = variables['highrg'][0]
    zflag = variables['zflag'][0]
    zcutoff	 = variables['zcutoff'][0]
    cflag = variables['cflag'][0]
    confile = variables['confile'][0]
    nexsegments1 = variables['nexsegments1'][0]
    nsegments1 = variables['nsegments1'][0]
    reslow1 = variables['reslow1'][0]
    numcont1 = variables['numcont1'][0]
    nexsegments2 = variables['nexsegments2'][0]
    nsegments2 = variables['nsegments2'][0]
    reslow2 = variables['reslow2'][0]
    numcont2 = variables['numcont2'][0]

    error = input_filter.check_name(runname)
    if(error!=[]):
        return error

    error=[]

    if 'no_file_check' not in kwargs:
        ev,rv,wv=input_filter.check_permissions(path)
        if(not ev or not rv or not wv):
            error.append('permission error in input file path '+path+'  [code = '+str(ev)+str(rv)+str(wv)+']')
            if(ev==False):
                error.append('path does not exist')
            elif(rv==False):
                error.append('read permission not allowed')
            elif(wv==False):
                error.append('write permission not allowed')
            return error

    pdbfile1=path+'/'+pdbmol1
    pdbfile2=path+'/'+pdbmol2

    print 'pdbfile1: ', pdbfile1

    error = input_filter.check_file_exists(pdbfile1)
    if(error!=[]):
        return error

    error = input_filter.check_file_exists(pdbfile2)
    if(error!=[]):
        return error

    ev,value=input_filter.check_pdb_dcd(pdbfile1,'pdb')

#ev == 0 not tested since non-existent file will trigger check_file_exists test above
    if(ev == 0):
        error.append('input pdb file, '+pdbfile1[3:]+', does not exist')
        return error
    if(value == 0):
        error.append( 'input pdb file, '+pdbfile1[3:]+', is not a valid pdb file')
        return error

    ev,value=input_filter.check_pdb_dcd(pdbfile2,'pdb')

#ev == 0 not tested since non-existent file will trigger check_file_exists test above
    if(ev == 0):
        error.append('input pdb file, '+pdbfile2[3:]+', does not exist')
        return error
    if(value == 0):
        error.append( 'input pdb file, '+pdbfile2[3:]+', is not a valid pdb file')
        return error

    m1=sasmol.SasMol(0)
    m1.read_pdb(pdbfile1)
    m2=sasmol.SasMol(1)
    m2.read_pdb(pdbfile2)

    segname1 = m1.segname()
    segname2 = m2.segname()
    segment_names_1 = string.split(nsegments1,',')
    segment_names_2 = string.split(nsegments2,',')

    if(accpos != 0 and accpos != 1):
        error.append( 'accept supplied position needs to be (0==no or 1==yes) : '+str(accpos))
        return error
    elif(len(pos)!=3):
        error.append( 'three float values are required for initial position (x,y,z) : '+str(pos))
        return error
    elif(len(trans)!=3):
        error.append( 'three int values are required for number of x,y,z moves : '+str(trans))
        return error
    elif(trans[0] < 1 or trans[1] < 1 or trans[2] < 1):
        error.append( 'you must specifiy at least ONE translational "move" for each axis : '+str(trans))
        return error
    elif(len(dtrans)!=3):
        error.append( 'three float values are required for dx,dy,dz step sizes : '+str(dtrans))
        return error
    elif(len(theta)!=3):
        error.append( 'three int values are required for theta angular moves : '+str(theta))
        return error
    elif(theta[0] < 1 or theta[1] < 1 or theta[2] < 1):
        error.append( 'you must specifiy at least ONE angular "move" for each axis : '+str(theta))
        return error
    elif(len(dtheta)!=3):
        error.append( 'three float values are required for dtheta (x,y,z) step sizes : '+str(dtheta))
        return error
    elif(basis!='CA'):
        error.append( 'only "CA" is accepted as a basis')
        return error
    elif(cutoff < 1.0):
        error.append( 'use a larger cutoff value, cutoff = '+str(cutoff))
        return error
    elif(zflag != 0 and zflag != 1):
        error.append( 'ERROR in Z coordinate filter selection: zflag == 0 for "no" and 1 for "yes", zflag = '+str(zflag))
        return error
    elif(cflag != 0 and cflag != 1):
        error.append( 'ERROR in atomic constraints selection: cflag == 0 for "no" and 1 for "yes", cflag = '+str(cflag))
        return error
    elif(cflag == 1):
        err = input_filter.check_file_exists(confile)
        if(err != []):
            lerr=['ERROR in constraint filename selection: ']
            lerr.append(err)
            error.append(lerr[0]+err[0])
            return error
        filter_flag = 1
        m3=sasmol.SasMol(2)
        err0 = m3.merge_two_molecules(m1,m2)
        err = constraints.read_constraints(m3,confile,filter_flag)
        if(err != []):
            error.append(err[0])
            return error
    elif(lowrg > highrg):
        error.append( 'low Rg cutoff is larger than high Rg cutoff, lowrg = '+str(lowrg)+' highrg = '+str(highrg))
        return error
    elif(lowrg < 0 or highrg < 0):
        error.append( 'Rg cutoffs need to be >= zero, lowrg = '+str(lowrg)+' highrg = '+str(highrg))
        return error
    elif(nexsegments1 < 0):
        error.append( 'number of excluded segments needs to be >= 0 (mol1) : '+str(nexsegments1))
        return error
    elif(nexsegments2 < 0):
        error.append( 'number of excluded segments needs to be >= 0 (mol2) : '+str(nexsegments2))
        return error
    elif(nexsegments1 > 0 and len(segment_names_1)!=nexsegments1):
        error.append( 'number of segment names does not match number of excluded segments (mol1) : '+str(nsegments1))
        return error
    elif(nexsegments2 > 0 and len(segment_names_2)!=nexsegments2):
        error.append( 'number of segment names does not match number of excluded segments (mol2) : '+str(nsegments2))
        return error
    elif(nexsegments1 > 0 and len(reslow1) != nexsegments1):
        error.append( 'the number of low residue values does not match the number of excluded segments (mol1), lowres1 = '+str(reslow1)+' nexsegments1 = '+str(nexsegments1))
        return error
    elif(nexsegments2 > 0 and len(reslow2) != nexsegments2):
        error.append( 'the number of low residue values does not match the number of excluded segments (mol2), lowres2 = '+str(reslow2)+' nexsegments2 = '+str(nexsegments2))
        return error
    elif(nexsegments1 > 0 and len(numcont1) != nexsegments1):
        error.append( 'the number of contiguous residues does not match the number of excluded segments (mol1), numcont1 = '+str(numcont1)+' nexsegments1 = '+str(nexsegments1))
        return error
    elif(nexsegments2 > 0 and len(numcont2) != nexsegments2):
        error.append( 'the number of contiguous residues does not match the number of excluded segments (mol2), numcont2 = '+str(numcont2)+' nexsegments2 = '+str(nexsegments2))
        return error

    return error




