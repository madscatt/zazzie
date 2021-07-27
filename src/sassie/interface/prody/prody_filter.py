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
import sasmol.sasmol as sasmol

def check_prody(variables,**kwargs):

    runname = variables['runname'][0]
    pdbfile = variables['pdbfile'][0]

    number_modes = variables['number_modes'][0]
    number_conformations_samp = variables['number_conformations_samp'][0]
    number_steps_traverse = variables['number_steps_traverse'][0]
    rmsd_conformations_samp = variables['rmsd_conformations_samp'][0]
    rmsd_traverse = variables['rmsd_traverse'][0]
    advanced_usage = variables['advanced_usage'][0]
    advanced_usage_cmd = variables['advanced_usage_cmd'][0]


    error=[]
    error = input_filter.check_name(runname)

    if(error!=[]):
        return error

    error=input_filter.check_file_exists(pdbfile)
    if(len(error) != 0):
#        error.append('input pdb file, '+pdbfile+', does not exist')        ##check_file_exists returns its own error
        return error
    ev,value=input_filter.check_pdb_dcd(pdbfile,'pdb')
    #ev == 0 not tested since file exists check is performed above
    if(ev == 0):
        error.append('check input pdb file: '+pdbfile)
        return error
    if(value == 0):
        error.append( 'input pdb file, '+pdbfile+', is not a valid pdb file')
        return error

    #NOT tested.  A PDB file with no frames is not a valid PDB file, so test fails above and doesn't get to this point.
    try:
        m1 = sasmol.SasMol(0)
        m1.read_pdb(pdbfile)
        number_of_frames = m1.number_of_frames()
        print '> found '+str(number_of_frames)+' frames in PDB file'
    except:
        error.append('could not open PDB file '+pdbfile+' to check number of frames')
        return error
    if(number_of_frames < 1):
        error.append('PDB file has no frames : '+pdbfile)
        return error

    if number_modes <= 0:
        error.append('number of normal modes needs to be greater than zero')
        return error
    elif number_conformations_samp <= 0:        
        error.append('number of conformations per sample needs to be greater than zero')
        return error
    elif number_steps_traverse <= 0:        
        error.append('number of frames to traverse per mode needs to be greater than zero')
        return error
    elif rmsd_conformations_samp <= 0:        
        error.append('average sampled RMSD for normal mode trajectory needs to be greater than zero')
        return error
    elif rmsd_traverse <= 0:        
        error.append('maximum sampled RMSD for traverse mode trajectory needs to be greater than zero')
        return error

    # advanced options

    # check ProDy command

    # check Rg cutoffs

    '''
    elif(lowrg > highrg):
        error.append( 'low Rg cutoff is larger than high Rg cutoff, lowrg = '+str(lowrg)+' highrg = '+str(highrg))
        return error
    elif(lowrg < 0 or highrg < 0):
        error.append( 'Rg cutoffs need to be >= zero, lowrg = '+str(lowrg)+' highrg = '+str(highrg))
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
        m1=sasmol.SasMol(0)
        m1.read_pdb(pdbfile)
        err = constraints.read_constraints(m1,confile,filter_flag)
        if(err != []):
            error.append(err[0])
            return error

    '''


    return error

