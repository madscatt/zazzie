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
import input_filter
import sasmol.sasmol as sasmol

def check_capriqorn(variables,**kwargs):
   
    error=[]
    
    runname = variables['runname'][0]
    pdbfile = variables['pdbfile'][0]
    dcdfile = variables['dcdfile'][0]

    number_q_values = variables['number_q_values'][0]
    q_max = variables['q_max'][0]

    if number_q_values < 1:
        error.append('number of q-values needs to be greater than zero')
        return error
    elif q_max <= 0:        
        error.append('q-max needs to be greater than zero')
        return error
    
    create_alias_flag = variables['create_alias_flag'][0]

    if not create_alias_flag:
        aliasfile = variables['aliasfile'][0]

        ### OPEN check if aliasfile exists and can be read

        ### OPEN check if aliasfile can be used with supplied pdbfile

    ### OPEN advanced options

    ### OPEN check number of gpu/cpu cores etc. 

    ### OPEN check alias.dat file 

    error = input_filter.check_name(runname)

    if(error!=[]):
        return error

    error=input_filter.check_file_exists(pdbfile)
    if(len(error) != 0):
        error.append('input pdb file, '+pdbfile+', does not exist')
        return error
    ev,value=input_filter.check_pdb_dcd(pdbfile,'pdb')
    if(ev == 0):
        error.append('check input pdb file: '+pdbfile)
        return error
    if(value == 0):
        error.append( 'input pdb file, '+pdbfile+', is not a valid pdb file')
        return error
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

    if dcdfile[-3:] == 'dcd':

        ev,value=input_filter.check_pdb_dcd(dcdfile,'dcd')
        value = 0
        if(ev == 1):   # if the file exists
            if(value == 0):         # not a pdb file
                #  'checking infile : as dcd'

                ev,value=input_filter.check_pdb_dcd(dcdfile,'dcd')
                if(value == 1):
                    cvalue=input_filter.certify_pdb_dcd(pdbfile,dcdfile)
                    if(cvalue == 0):
                        error.append('input pdb file '+pdbfile+' and dcd file '+dcdfile+' are not compatible (different number of atoms)')
                        return error
                else:
                    error.append('dcd input file '+dcdfile+' is not a valid dcd file')
                    return error

        xstfile = variables['xstfile'][0]

        ### OPEN check if xstfile exists and can be read
        
        ### OPEN  check if xstfile can be used with supplied dcdfile

    elif dcdfile[-6:] == 'crdbox':

        pass

        ### OPEN check if crdbox exists and can be read

        ### OPEN check if crdbox can be used with supplied pdbfile

        ### OPEN check if box dimensions are readable and have ~valid values 

    else:
        error.append('infile needs to have a "dcd" or "crdbox" suffix')
        return error


    return error

