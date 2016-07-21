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
import os,sys,locale,string,glob
import sasmol.sasmol as sasmol
import sassie.simulate.constraints.constraints as constraints
#import input_filter
import sassie.interface.input_filter as input_filter

def check_merge_utilities(variables,**kwargs):
    runname             = variables['runname'][0]
    pdb_file            = variables['pdb_file'][0]
    merge_option        = variables['merge_option'][0]
    number_of_runs      = variables['number_of_runs'][0]
    trajectory_names    = variables['trajectory_names'][0]
    sas_type            = variables['sas_type'][0]
    sas_paths           = variables['sas_paths'][0]
    merge_type_option   = variables['merge_type_option'][0]
    local_value         = variables['local_value'][0]
    output_filename     = variables['output_filename'][0]

    error=[]
    error = input_filter.check_name(runname)
    if(error!=[]):
        return error

    ## @NOTE to ZHL: do I really need to do the following check for options/types?
    if merge_option not in [0,1,2]:
        error.append('merge option needs to be 0, 1, or 2, you entered : '+str(merge_option))
        return error
    if merge_type_option not in [0,1,2]:
        error.append('merge_type option needs to be 0, 1, or 2, you entered : '+str(merge_type_option))
        return error
    if sas_type not in [1,2,3]:
        error.append('sas type needs to be 1, 2, or 3, you entered : '+str(sas_type))
        return error

    ## check for pdb/dcd files
    if merge_option in [0,2]:
        #
        try:
            if (output_filename[0] == '.' or output_filename[-4:] not in ['.pdb','.dcd'] or len(output_filename) < 5):
                error.append('output filename must be greater than four characters long and end with .pdb or .dcd : '+output_filename)
                return error
        except:
            return error

        trajectory_name_list = trajectory_names.split(',') ## @NOTE to ZHL: !
        if len(trajectory_name_list) != number_of_runs:
            error.append('number of trajectory files "%d" does not match number of runs "%d"'%(len(trajectory_name_list),number_of_runs))
            return error
        #
        error=input_filter.check_file_exists(pdb_file)
        if(len(error) != 0):
            #error.append('input pdb file, '+pdb_file+', does not exist')
            return error
        ev,value=input_filter.check_pdb_dcd(pdb_file,'pdb')
        if(ev == 0):
            error.append('check input pdb file: '+pdb_file)
            return error
        if(value == 0):
            error.append( 'input pdb file, '+pdb_file+', is not a valid pdb file')
            return error
        try:
            m1 = sasmol.SasMol(0)
            m1.read_pdb(pdb_file)
            number_of_frames = m1.number_of_frames()
        except:
            error.append('could not open PDB file '+pdb_file+' to check number of frames')
            return error
        if(number_of_frames < 1):
            error.append('PDB file has no frames : '+pdb_file)
            return error
        #
        number_of_frames_list = []
        for trajectory_file in trajectory_name_list:
            error=input_filter.check_file_exists(trajectory_file)
            if(len(error) != 0):
                #error.append('input trajectory file, '+trajectory_file+', does not exist')
                return error
            if trajectory_file[-3:] == 'dcd':
                infile_type = 'dcd'
                ev,value=input_filter.check_pdb_dcd(trajectory_file,'dcd')
                if(ev == 0):
                    error.append('check input trajectory filename : '+trajectory_file)
                    return error
                elif(value == 0):
                    error.append( 'input trajectory file, '+trajectory_file+', is not a valid dcd file')
                    return error

                value=input_filter.certify_pdb_dcd(pdb_file,trajectory_file)
                if(value == 0):
                    error.append('input pdb file '+pdb_file+' and dcd file '+trajectory_file+', are not compatible')
                    return error

            elif trajectory_file[-3:] == 'pdb':
                infile_type = 'pdb'
                ev,value=input_filter.check_pdb_dcd(trajectory_file,'pdb')
                if(ev == 0):
                    error.append('check input trajectory filename : '+trajectory_file)
                    return error
                elif(value == 0):
                    error.append( 'input trajectory file, '+trajectory_file+', is not a valid pdb file')
                    return error

                value=input_filter.certify_pdb_pdb(pdb_file,trajectory_file)
                if(value == 0):
                    error.append('input pdb file '+pdb_file+' and pdb file '+trajectory_file+', are not compatible')
                    return error

            try:
                m1 = sasmol.SasMol(0)
                if infile_type == 'dcd':
                    dcdinputfile = m1.open_dcd_read(trajectory_file)
                    number_of_frames = dcdinputfile[2]
                elif infile_type == 'pdb':
                    m1.read_pdb(trajectory_file)
                    number_of_frames = m1.number_of_frames()
            except:
                error.append('could not open trajectory file '+trajectory_file+' to check number of frames')
                return error
            if(number_of_frames < 1):
                error.append('trajectory file has no frames : '+trajectory_file)
                return error
            number_of_frames_list.append(number_of_frames)

    ## check for sas files
    if merge_option in [0,1]:
        #
        sas_path_list = sas_paths.split(',') ## @NOTE to ZHL: !
        if len(sas_path_list) != number_of_runs:
            error.append('number of SAS folders "%d" does not match number of runs "%d"'%(len()))
            return error
        #
        if sas_type not in [1,2,3]:
            error.append("sas_type %d not supported!"%sas_type)
            return error
        #
        number_spec_files_list = []
        for i,sas_path in enumerate(sas_path_list):
            head, my_sas_type = os.path.split(sas_path)
            dict_sas_type = {1:'sascalc', 2:'cryson', 3:'crysol'}
            name_sas_type = dict_sas_type[sas_type]
            if my_sas_type != name_sas_type:
                error.append('the SAS type "'+name_sas_type+ '" you entered is not compatiable with the SAS type in the SAS data path you selected')
                return error

            ev,rv,wv=input_filter.check_permissions(sas_path)
            if(not ev or not rv):
                error.append('permission error in input file path '+sas_path+'  [code = '+str(ev)+str(rv)+str(wv)+']')
                if(ev==False):
                    error.append('sas path "'+sas_path+'" does not exist')
                elif(rv==False):
                    error.append('read permission not allowed for sas path "'+sas_path+'"')
                return error
            #
            if(sas_type==1):
                suffix=['*.iq','*.log']
                extra = ['*.inf','*.crd','*.ans','*.pr']
            elif(sas_type==2):
                suffix=['*.int','*.log']
                extra = ['*.sav','*.flm','*.alm','*.ans']
            elif(sas_type==3):
                suffix=['*.int','*.log']
                extra = ['*.sav','*.flm','*.alm','*.ans']
            spec_files = glob.glob(os.path.join(sas_path,suffix[0]))
            log_files = glob.glob(os.path.join(sas_path,suffix[1]))
            #spec_files.sort()
            #log_files.sort()
            extra_files_list = []
            for ex in extra:
                #extra_files_list.append(glob.glob(os.path.join(sas_path,suffix[0])).sort())
                extra_files_list.append(glob.glob(os.path.join(sas_path,suffix[0])))
            number_spec_files_list.append(len(spec_files))
            #
            if merge_option == 0:
                if number_spec_files_list[i] != number_of_frames_list[i]:
                    error.append('number of SAS files in folder "%s" does not match number of frames of the PDB/DCD file in "%s"'%(sas_path, trajectory_name_list[i]))

                if number_spec_files_list[i]==0:
                    error.append("there is no scattering files found for the selected sas-type: '"+name_sastype+"' in folder: '"+sas_path+"'")
                    return error

    ## some general naming
    if merge_option == 0:
        filetype = 'PDB/DCD and/or SAS'
        unit = 'frame and/or SAS file'
        units = 'frames and/or SAS files'
        data_place_list = trajectory_name_list
    elif merge_option == 1:
        filetype = 'SAS'
        unit = 'SAS file'
        units = 'SAS files'
        number_of_frames_list = number_spec_files_list
        data_place_list = sas_path_list
    elif merge_option == 2:
        filetype = 'PDB/DCD'
        unit = 'frame'
        units = 'frames'
        data_place_list = trajectory_name_list

    ## checking merge_type options
    if(merge_type_option == 1):
        weight_file_list = local_value.split(',')
        if len(weight_file_list) != number_of_runs:
            error.append('number of weight files "%s" does not match number of runs "%s"'%(len(weight_file_list),number_of_runs))
            return error
        for j,weight_file in enumerate(weight_file_list):
            number_of_frames = number_of_frames_list[j]
            try:
                error=input_filter.check_file_exists(weight_file)
                if(len(error) != 0):
                    #error.append('weight file does not exist :'+weight_file)
                    return error
            except:
                error.append('weight file does not exist :'+weight_file)
                return error
            #try:
            if True:
                infile = open(weight_file,'r').readlines()
                frames = [] ; weights = []
                for lin in infile:
                    words = lin.split()
                    if not len(words): continue
                    if words[0][0]=='#': continue
                    else:
                        if locale.atof(words[2]) not in [0.0,1.0]:
                            error.append('weight file column 3 can only have 0 or 1 : '+words[2]+' was found')
                            return error
                        weights.append(locale.atof(words[2]))
                        frames.append(locale.atoi(words[0]))
                frames.sort()
                local_count = 0

                for this_weight in weights:
                    if this_weight == 0.0:
                        local_count += 1

                if local_count == len(weights):
                    error.append('all weights in weight file are zero which means you will not extract any structure or SAS profiles')
                    return error


                for i in range(len(frames)):
                    if(frames[i] < 1):
                        error.append('weight file colume 1 can only have positive integers : "'+str(frames[i])+'" was found in the weight file')
                        return error
                    if(frames[i] > number_of_frames):
                        error.append('there are '+str(number_of_frames)+' '+units+' in "'+data_place_list[j]+'" : '+unit+' number '+str(frames[i])+' was found in the weight file')
                        return error
                    if(i>0 and frames[i]==frames[i-1]):
                        error.append('redundant '+unit+' number "'+str(frames[i])+'" found in the weight file')
                        return error
            #except:
            else:
                error.append('encountered an unknown error reading weight_file: '+weight_file)
                return error
    elif(merge_type_option == 2):
        try:
            sampling_frequency = locale.atoi(local_value)
            if sampling_frequency < 1:
                error.append('the sampling frequency needs to be a positive number : you entered "'+local_value+'"')
                return error
            for j,number_of_frames in enumerate(number_of_frames_list):
                if sampling_frequency > number_of_frames:
                    error.append('the sampling frequency need to be smaller than the number of '+units+' in "'+data_place_list[j]+'" ("'+str(number_of_frames)+'") : you entered "'+local_value+'"')
                    return error
        except:
            error.append('encountered an unknown error reading the sampling frequency : '+local_value)
            return error

    return error
