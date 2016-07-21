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
import os,sys,string,locale
import input_filter

def check_chi_square(variables,**kwargs):

    runname			= variables['runname'][0]           
    saspath			= variables['saspath'][0]           
    sasintfile		= variables['sasintfile'][0]        
    io			= variables['io'][0]                
    x2highcut		= variables['x2highcut'][0]         
    x2highweightfile	= variables['x2highweightfile'][0]  
    x2lowcut		= variables['x2lowcut'][0]          
    x2lowweightfile		= variables['x2lowweightfile'][0]   
    rghighcut		= variables['rghighcut'][0]         
    rghighweightfile	= variables['rghighweightfile'][0]  
    rglowcut		= variables['rglowcut'][0]          
    rglowweightfile		= variables['rglowweightfile'][0] 
    sastype			= variables['sastype'][0]           

    error=[]
    error = input_filter.check_name(runname)
    if(error!=[]):
        return error
    
    if "no_file_check" not in kwargs:

        ev,rv,wv=input_filter.check_permissions(saspath)

        if(not ev or not rv or not wv):
            error.append('permission error in input file path '+saspath+'  [code = '+str(ev)+str(rv)+str(wv)+']')
            if(ev==False):
                error.append('path does not exist')
            elif(rv==False):
                error.append('read permission not allowed')
            elif(wv==False):
                error.append('write permission not allowed')
            return error

    if(x2highcut < 0.0 or x2lowcut < 0.0):
        error.append('x2 cutoffs need to be > 0, x2highcut = '+str(x2highcut)+' x2lowcut = '+str(x2lowcut))
        return error
    elif(x2highcut < x2lowcut):
        error.append('x2highcut < x2lowcut, x2highcut = '+str(x2highcut)+' x2lowcut = '+str(x2lowcut))
        return error
    elif(rghighcut < 0.0 or rglowcut < 0.0):
        error.append('Rg cutoffs need to be > 0, Rghighcut = '+str(rghighcut)+' Rglowcut = '+str(rglowcut))
        return error
    elif(rghighcut < rglowcut):
        error.append('Rghighcut < Rglowcut, Rghighcut = '+str(rghighcut)+' Rglowcut = '+str(rglowcut))
        return error
    elif(sastype != 1 and sastype != 2 and sastype !=3):
        error.append('sastype needs to be equal to 1, 2, or 3 : sastype = '+str(sastype))
        return error

    error=input_filter.check_file_exists(sasintfile)
    if(len(error) != 0):
        error.append('check input SAS exp. data file path + filename : '+sasintfile)
        return error

### check interpolated sasintfile for three columns (q,I(q),error)
#	nsfit.data 
#	0.000000        0.019000        0.000824

    try:
        infile = open(sasintfile,'r').readlines()
        if(len(infile)<1):
            error.append("no lines in your SAS data file : "+sasintfile)
            return error
        else:
            sum = 0.0
            for i in xrange(len(infile)):
                lin = string.split(infile[i])
                if(lin[0][0] != "#"):
                    this_value = locale.atof(lin[0])
                    this_value = locale.atof(lin[2])
                    this_value = locale.atof(lin[1])
                    sum += abs(this_value)
            if(sum == 0.0):
                error.append("all I(q) values in SAS data file are zero : "+sasintfile)
    except:
        print '>>> ERROR IN SAS DATA FILE'
        error.append("unable to open and read your SAS data file : "+sasintfile)

    return error

