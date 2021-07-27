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

def check_sascalc(variables,**kwargs):

    runname = variables['runname'][0]

    error=[]
    error = input_filter.check_name(runname)

    pdbfile = variables['pdbfile'][0]
    dcdfile = str(variables['dcdfile'][0])
    xon = variables['xon'][0]
    if xon not in ['neutron','xray','neutron_and_xray']:
        error.append("Either neutron or xray input need to be checked")
        return error

    number_contrast_points = variables['number_contrast_points'][0]
    D2O_percentage_array = variables['D2O_percentage_array'][0]
    I0_array = variables['I0_array'][0]
    #
    if xon in ['neutron','neutron_and_xray']:
        number_exH_regions = variables['number_exH_regions'][0]
        exH_basis_string_array = variables['exH_basis_string_array'][0]
        fraction_exH_array = variables['fraction_exH_array'][0]
        #
        number_deuteration_regions = variables['number_deuteration_regions'][0]
        deuterated_basis_string_array = variables['deuterated_basis_string_array'][0]
        fraction_deuterated_array = variables['fraction_deuterated_array'][0]
    #
    xray_number_contrast_points = variables['xray_number_contrast_points'][0]
    xray_D2O_percentage_array = variables['xray_D2O_percentage_array'][0]
    xray_I0_array = variables['xray_I0_array'][0]
        
    number_q_values = variables['number_q_values'][0]
    q_max = variables['q_max'][0]
    number_r_values = variables['number_r_values'][0]
    golden_vector_method_option = variables['golden_vector_method_option'][0]
    if golden_vector_method_option == 'fixed':
        number_golden_vectors = variables['number_golden_vectors'][0]
    elif golden_vector_method_option == 'converge':
        golden_vector_method_converge_tolerance = variables['golden_vector_method_converge_tolerance'][0]

    solvent_volume = variables['solvent_volume'][0]
    VDW_scaling_factor = variables['VDW_scaling_factor'][0]

    if(error!=[]):
        return error

    error=input_filter.check_file_exists(pdbfile)
    if(len(error) != 0):
        #error.append('input pdb file, '+pdbfile+', does not exist')
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

    error=input_filter.check_file_exists(dcdfile)
    if(len(error) != 0):
        #error.append('input dcd file, '+dcdfile+', does not exist')
        return error
    '''
    ev,value=input_filter.check_pdb_dcd(dcdfile,'dcd')
    if(ev == 0):
        error.append('check input dcd filename : '+dcdfile)
        return error
    elif(value == 0):
        ev,value=input_filter.check_pdb_dcd(dcdfile,'pdb')
        if(value ==0):
            error.append( 'input dcd file, '+dcdfile+', is not a valid pdb file')
            return error
        infile_type = 'pdb'
    else:
        infile_type = 'dcd'

    if(infile_type == 'dcd'):
        value=input_filter.certify_pdb_dcd(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and dcd file '+dcdfile+', are not compatiable')
            return error
    elif(infile_type == 'pdb'):
        value=input_filter.certify_pdb_pdb(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and dcd file '+dcdfile+', are not compatiable')
            return error
    '''

    if dcdfile[-3:] == 'dcd':
        infile_type = 'dcd'
        ev,value=input_filter.check_pdb_dcd(dcdfile,'dcd')
        if(ev == 0):
            error.append('check input dcd filename : '+dcdfile)
            return error
        elif(value == 0):
            error.append( 'input dcd file, '+dcdfile+', is not a valid dcd file')
            return error

        value=input_filter.certify_pdb_dcd(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and dcd file '+dcdfile+', are not compatible')
            return error

    elif dcdfile[-3:] == 'pdb':
        infile_type = 'pdb'
        ev,value=input_filter.check_pdb_dcd(dcdfile,'pdb')
        if(ev == 0):
            error.append('check input dcd filename : '+dcdfile)
            return error
        elif(value == 0):
            error.append( 'input dcd file, '+dcdfile+', is not a valid pdb file')
            return error

        value=input_filter.certify_pdb_pdb(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and pdb file '+dcdfile+', are not compatible')
            return error

    try:
        m1 = sasmol.SasMol(0)
        if(infile_type == 'dcd'):
            dcdinputfile = m1.open_dcd_read(dcdfile)
            number_of_frames = dcdinputfile[2]
        elif(infile_type == 'pdb'):
            m1.read_pdb(dcdfile)
            number_of_frames = m1.number_of_frames()
    except:
        error.append('could not open dcd file '+dcdfile+' to check number of frames')
        return error
    if(number_of_frames < 1):
        error.append('dcd file has no frames : '+dcdfile)
        return error

    if xon not in ['xray', 'neutron', 'neutron_and_xray']:
        error.append("scattering source can only be 'xray', 'neutron', or 'neutron_and_xray'; '"+xon+"' is provided")
        return error

    if xon in ['neutron','xray','neutron_and_xray']:
        if number_contrast_points<=0:
            error.append("number of contrast points must be an integer")
            return error
        if len(D2O_percentage_array) != number_contrast_points:
            error.append("len of the D2O_percentage_array does not match number of contrast points")
            return error
        if len(I0_array) != number_contrast_points:
            error.append("len of the I0_array does not match number of contrast points")
            return error
        if len(D2O_percentage_array) != len(set(D2O_percentage_array)):
            error.append("You have duplicate D2O percentage")
            return error
        for D2O_percentage in D2O_percentage_array:
            if not 0.0<=D2O_percentage<=100.0:
                error.append("D2O percentage must be within [0.0, 100.0]")
                return error
        for I0 in I0_array:
            if I0<0.0:
                error.append("I0 value cannot be smaller than 0")
                return error
        if xon in ['neutron','neutron_and_xray']:
            if number_exH_regions<0:
                error.append("number of exH regions must be a non-negative integer")
                return error
            if number_exH_regions and len(fraction_exH_array) != number_exH_regions:
                error.append("len of the fraction_exH_array does not match number of exH regions")
                return error
            if number_deuteration_regions<0:
                error.append("number of deuterated regions must be a non-negative integer")
                return error
            if number_deuteration_regions and len(fraction_deuterated_array) != number_deuteration_regions:
                error.append("len of the fraction_deuterated_array does not match number of deuteration regions")
                return error
            for fraction_exH in fraction_exH_array:
                if not 0.0<=fraction_exH<=1.0:
                    error.append("fraction_exH must be within [0.0, 1.0]")
                    return error
            for fraction_deuterated in fraction_deuterated_array:
                if not 0.0<=fraction_deuterated<=1.0:
                    error.append("fraction_deuterated must be within [0.0, 1.0]")
                    return error

    if number_q_values<=0:
        error.append("number of q values is smaller than or equal to 0")
        return error

    if q_max<=0.0:
        error.append("maximum q value is smaller than or equal to 0")
        return error

    if number_r_values<=0:
        error.append("number of r values is smaller than or equal to 0")
        return error

    if golden_vector_method_option not in ['fixed','converge']:
        error.append("golden vector option: '"+golden_vector_method_option+"' not recognized")
        return

    if golden_vector_method_option=='fixed':
        if number_golden_vectors<=0:
            error.append("number of golden vectors is smaller than or equal to 0")
            return error

    if golden_vector_method_option=='converge':
        if golden_vector_method_converge_tolerance<=0.0:
            error.append("convergence tolerance value for golden vector method is smaller than or equal to 0")
            return error

    if solvent_volume<=0.0:
        error.append("solvent volume is smaller than or equal to 0")
        return error
    if VDW_scaling_factor<=0.0:
        error.append("VDW scaling factor is smaller than or equal to 0")
        return error

    ## NOTE to developers: directed_mc not implemented yet
    #directed_mc = variables['directed_mc'][0]


    return error

