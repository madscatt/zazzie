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

def check_apbs(variables,**kwargs):

    runname        = variables['runname'][0]
    infile         = variables['infile'][0]
    pdbfile        = variables['pdbfile'][0]
    temperature    = variables['temperature'][0]
    ph        = variables['ph'][0]
    ion_charge    = variables['ion_charge'][0]
    ion_conc    = variables['ion_conc'][0]
    ion_radius    = variables['ion_radius'][0]
    manual_flag    = variables['manual_flag'][0]
    manual_file    = variables['manual_file'][0]

    error=[]
    error = input_filter.check_name(runname)
    if(error!=[]):
        return error

#   path is not an input variable; this test isn't needed
#    path = './'
#
#    ev,rv,wv=input_filter.check_permissions(path)
#    if(not ev or not rv or not wv):
#        error.append('permission error in input file path '+path+'  [code = '+str(ev)+str(rv)+str(wv)+']')
#        if(ev==False):
#            error.append('path does not exist')
#        elif(rv==False):
#            error.append('read permission not allowed')
#        elif(wv==False):
#            error.append('write permission not allowed')
#        return error

#   pdbfile and infile existence is tested below
#    error = input_filter.check_file_exists(pdbfile)
#    if(error!=[]):
#       return error

#    error = input_filter.check_file_exists(infile)
#    if(error!=[]):
#        return error

    ev,value=input_filter.check_pdb_dcd(pdbfile,'pdb')

    if(ev == 0):
        error.append('reference pdb file, '+pdbfile[3:]+', does not exist')
        return error
    if(value == 0):
        error.append( 'reference pdb file, '+pdbfile[3:]+', is not a valid pdb file')
        return error

### now check and see if input file is a correct pdb or dcd     

    ev,value=input_filter.check_pdb_dcd(infile,'pdb')
#    print 'ev,value: ', ev,value
    if(ev == 1):
        if(value == 0):         # not a pdb file
            ev,value=input_filter.check_pdb_dcd(infile,'dcd')
            if(value == 1):
                cvalue=input_filter.certify_pdb_dcd(pdbfile,infile)
                if(cvalue == 0):
                    error.append('input pdb '+pdbfile[3:]+' and dcd file '+infile[3:]+' are not compatible')
                    return error
            else:
                error.append('input file, '+infile[3:]+', is not a valid pdb or dcd file')
                return error
        else:                   # is a pdb file
            variables=['name']
            value,result1=input_filter.get_pdb_stats(infile,variables)
            value,result2=input_filter.get_pdb_stats(pdbfile,variables)
            if(result1 != result2):
                error.append('reference pdb file '+pdbfile[3:]+' and input pdb file '+infile[3:]+' are not compatible')
                return error
    else:
        error.append('file : '+infile[3:]+' does not exist')
        return error

    if(temperature <= 0):
        error.append('temperature needs to be greater than zero : '+str(temperature))
        return error
    elif(ph < 0):
        error.append('ph needs to be greater than or equal to zero : '+str(ph))
        return error
    elif(ion_conc < 0):
        error.append('ion concentration needs to be greater than or equal to zero : '+str(ion_conc))
        return error
    elif(ion_radius <= 0):
        error.append('ion radius needs to be greater than zero : '+str(ion_radius))
        return error
# feature not yet implemented
#    elif(manual_flag != 0 and manual_flag != 1):
#        error.append('manual flag needs to be 0 or 1 : '+str(manual_flag))
#        return error
#    if(manual_flag == 1):
#        if(len(manual_file)<0):    
#            error.append('incorrect user input filename : '+manual_flag)
#            return error
#
#        error = input_filter.check_name(manual_file)
#        if(error!=[]):
#            return error
#
#        error = input_filter.check_file_exists(manual_file)
#        if(error!=[]):
#            return error

    return error




