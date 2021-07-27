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
import sassie.interface.input_filter as input_filter

def check_experimental_data(datafile):

	error = []
	inputfile = open(datafile,'r').readlines()
	number_of_lines = len(inputfile)

	z_values = []

	for i in xrange(number_of_lines):
		lin = string.split(inputfile[i])
		if(len(lin)>0):
			if(lin[0][0] != "#"):
				this_value = locale.atof(lin[0])
				if(i == 0):
					z_values.append(this_value)
				elif(this_value in z_values):
					message = 'repeated z-value in experimental data: '+str(this_value)+'\n'
					error.append(message)	
				else:
					z_values.append(this_value)
			
	return error,z_values

def check_sld_mol(variables,**kwargs):

    runname		= variables['runname'][0]          
    path		= variables['path'][0]              
    pdbfile		= variables['pdbfile'][0]           
    dcdfile		= variables['dcdfile'][0]           
    expdatafile	= variables['expdatafile'][0]	
    outputfile	= variables['outputfile'][0]	

    runtype		= variables['runtype'][0]		
    bulk_sld	= variables['bulk_sld'][0]		

    xon		= variables['xon'][0]		
    num_deut_regions= variables['num_deut_regions'][0]	
    deut_low_res	= variables['deut_low_res'][0]	
    deut_high_res	= variables['deut_high_res'][0]	

    dbin		= variables['dbin'][0]		
    width		= variables['width'][0]		

    sldfit		= variables['sldfit'][0]		
    sldoffset	= variables['sldoffset'][0]		

    zfit0		= variables['zfit0'][0]		
    zfitmin		= variables['zfitmin'][0]		
    zfitmax		= variables['zfitmax'][0]		
    zevalmin	= variables['zevalmin'][0]		
    zevalmax	= variables['zevalmax'][0]		

    A0		= variables['A0'][0]		
    Amin		= variables['Amin'][0]		
    Amax 		= variables['Amax'][0]		

    plotflag	= variables['plotflag'][0]		

    error=[]
    error = input_filter.check_name(runname)
    if(error!=[]):
        return error

    if "no_file_check" not in kwargs:
        dcdfile=path+'/'+dcdfile
        pdbfile=path+'/'+pdbfile
	
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

    if(runtype != 0 and runtype != 1):
        error.append('run type entered needs to be either 0 or 1 : '+str(runtype))
        return error
    elif(xon != 0 and xon != 1):
        error.append('scattering type (xon) entered needs to be either 0 or 1 : '+str(xon))
        return error
    elif(plotflag != 0 and plotflag != 1 and plotflag != 2):
        error.append('plotflag needs to be 0, 1, or 2 : '+str(plotflag))
        return error
    elif(num_deut_regions < 0):
        error.append('number of deuterated regions needs to be >= 0 : '+str(num_deut_regions))
        return error
    elif((len(deut_low_res) != num_deut_regions) and num_deut_regions > 0):
        error.append('number of low deuterated values does not match the number of regions: len(deut_low_res) = '+str(len(deut_low_res))+' num_deut_regions = '+str(num_deut_regions))
        return error
    elif((len(deut_high_res) != num_deut_regions) and num_deut_regions > 0):
        error.append('number of high deuterated values does not match the number of regions: len(deut_high_res) = '+str(len(deut_high_res))+' num_deut_regions = '+str(num_deut_regions))
        return error

    error=input_filter.check_file_exists(expdatafile)
    if(len(error)>0):
        return error

    error,z_values = check_experimental_data(expdatafile)
    if(len(error)>0):
        return error

    minz = min(z_values)+sldoffset ; maxz = max(z_values)+sldoffset

    if(zevalmin < minz):
        error.append('minimum evaluation value is less than experimental data + sldoffset: minz = '+str(minz)+' : zevalmin = '+str(zevalmin))
        return error
    elif(zevalmax > maxz):
        error.append('maximum evaluation value is greater than experimental data + sldoffset: maxz = '+str(maxz)+' : zevalmax = '+str(zevalmax))
        return error

    error=input_filter.check_file_exists(dcdfile)
    if(len(error)>0):
        return error

    error=input_filter.check_file_exists(pdbfile)
    if(len(error)>0):
        return error

    ev,value=input_filter.check_pdb_dcd(dcdfile,'pdb')
    if(ev == 1):
        if(value == 0):         # not a pdb file
            ev,value=input_filter.check_pdb_dcd(dcdfile,'dcd')
            if(ev == 1 and value == 1):
                cvalue=input_filter.certify_pdb_dcd(pdbfile,dcdfile)
                if(cvalue == 0):
                    error.append('pdbfile '+pdbfile+' and dcdfile '+dcdfile+' are not compatible (different number of atoms)')
                    return error
            else:
                error.append('input file '+dcdfile+' is not a valid pdb or dcd file or it does not exist')
                return error
        else:                   # is a pdb file
            variables=['name']
            value,result1=input_filter.get_pdb_stats(dcdfile,variables)
            value,result2=input_filter.get_pdb_stats(pdbfile,variables)
            if(result1 != result2):
                error.append('reference pdb file '+pdbfile+' and input pdb file '+dcdfile+' are not compatible')

                return error
#   not tested; non-existent dcd file is detected in check_file_exists above
    else:
        error.append('file '+dcdfile+' does not exist')
        return error

    if(dbin <= 0.0 or dbin > 1.0):
        error.append('bin width needs to be greater than zero and less than 1.0 angstroms: '+str(dbin))
        return error
    elif(width <= 0.0 or width > 5.0):
        error.append('smoothing width needs to be greater than zero and less than 5.0 angstroms: '+str(width))
        return error

    elif(sldfit != 0 and sldfit != 1):
        error.append('fit SLD needs to be either 0 or 1 :'+str(sldfit))
        return error

    elif((A0 <= 0.0 or Amin <= 0.0 or Amax <= 0.0) or (A0 > 1.0 or Amin > 1.0 or Amax > 1.0 )):
        error.append('surface coverage has to be a positive value less than or equal to 1.0: A0: '+str(A0)+' Amin :'+str(Amin)+' Amax: '+str(Amax))
        return error

    elif(Amin >= Amax):
        error.append('surface coverage maximum has to be greater than surface coverage minumum: Amin :'+str(Amin)+' Amax: '+str(Amax))
        return error

    elif(Amin > A0 or Amax < A0):
        error.append('surface coverage must be between Amin and Amax: A0 :'+str(A0)+' Amin: '+str(Amin)+' Amax: '+str(Amax))
        return error  

    elif(zevalmin < zfitmin):
        error.append('error evaluation minimum has to be greater than zfit minimum: zevalmin = '+str(zevalmin)+' zevalmax = '+str(zevalmax))
        return error

    elif(zfitmin > zfitmax):
        error.append('zfit maxiumum has to be greater than zfit minimum: zfitmin = '+str(zfitmin)+' zfitmax = '+str(zfitmax))
        return error

    elif(zfitmin > zfit0 or zfitmax < zfit0):
        error.append('zfit value must be between zfit maxiumum and zfit minimum: zfit0 :'+str(zfit0)+' zfitmin: '+str(zfitmin)+' zfitmax: '+str(zfitmax))

    elif(zevalmin >= zevalmax):
        error.append('error evaluation maximum has to be greater than error evaluation minumum: zevalmin = '+str(zevalmin)+' zevalmax = '+str(zevalmax))
        return error

    return error

