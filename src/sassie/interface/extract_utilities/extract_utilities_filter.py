import os
import sys
import locale
import glob
import re
import sasmol.system as system
import sassie.simulate.constraints.constraints as constraints
import sassie.interface.input_filter as input_filter

def check_extract_utilities(variables, **kwargs):
    run_name = variables['run_name'][0]
    pdb_filename = variables['pdb_filename'][0]
    trajectory_filename = variables['trajectory_filename'][0]
    option = variables['option'][0]
    local_value = variables['local_value'][0]
    output_filename = variables['output_filename'][0]
    extract_trajectory = variables['extract_trajectory'][0]
    extract_sas = variables['extract_sas'][0]
    path = variables['path'][0]
    sas_type = variables['sas_type'][0]
    sas_paths = variables['sas_paths'][0]

    error = []
    error = input_filter.check_name(run_name)
    if error:
        return error

    if 'no_file_check' not in kwargs:
        ev, rv, wv = input_filter.check_permissions(path)
        if not ev or not rv or not wv:
            error.append(f'permission error in input file path {path}  [code = {ev}{rv}{wv}]')
            if not ev:
                error.append('path does not exist')
            elif not rv:
                error.append('read permission not allowed')
            elif not wv:
                error.append('write permission not allowed')
            return error

    if not extract_trajectory and not extract_sas:
        error.append("at least one of the options ('extract trajectory' and 'extract SAS') needs to be checked")
        return error

    # checking trajectory data
    if extract_trajectory:
        try:
            if output_filename[0] == '.' or output_filename[-4:] not in ['.pdb', '.dcd'] or len(output_filename) < 5:
                error.append(f'output filename must be greater than four characters long and end with .pdb or .dcd : {output_filename}')
                return error
        except Exception as e:
            return error

        error = input_filter.check_file_exists(pdb_filename)
        if error:
            return error
        ev, value = input_filter.check_pdb_dcd(pdb_filename, 'pdb')
        if value == 0:
            error.append(f'input pdb file, {pdb_filename}, is not a valid pdb file')
            return error

        try:
            m1 = system.Molecule(0)
            m1.read_pdb(pdb_filename)
            number_of_frames = m1.number_of_frames()
            print(f'> found {number_of_frames} frames in PDB file')
        except Exception as e:
            error.append(f'could not read frames in pdb file {pdb_filename}')
            return error

        error = input_filter.check_file_exists(trajectory_filename)
        if error:
            return error
        ev, value = input_filter.check_pdb_dcd(trajectory_filename, 'dcd')
        if value == 0:
            ev, value = input_filter.check_pdb_dcd(trajectory_filename, 'pdb')
            if value == 0:
                error.append(f'input trajectory file, {trajectory_filename}, is not a valid pdb file')
                return error
            infile_type = 'pdb'
        else:
            infile_type = 'dcd'

        if infile_type == 'dcd':
            value = input_filter.certify_pdb_dcd(pdb_filename, trajectory_filename)
            if value == 0:
                error.append(f'input pdb file {pdb_filename} and dcd file {trajectory_filename}, are not compatible')
                return error
        elif infile_type == 'pdb':
            ev, value = input_filter.certify_pdb_pdb(pdb_filename, trajectory_filename)
            if value == 0:
                error.append(f'input pdb file {pdb_filename} and dcd file {trajectory_filename}, are not compatible')
                return error

        try:
            m1 = system.Molecule(0)
            if infile_type == 'dcd':
                dcdinputfile = m1.open_dcd_read(trajectory_filename)
                number_of_frames = dcdinputfile[2]
            elif infile_type == 'pdb':
                m1.read_pdb(trajectory_filename)
                number_of_frames = m1.number_of_frames()
        except Exception as e:
            error.append(f'could not read frames in trajectory file {trajectory_filename}')
            return error

    # checking sas data
    if extract_sas:
        if sas_type not in [0, 1, 2, 3]:
            error.append(f"sas_type {sas_type} not supported!")
            return error
        sas_paths = [x.strip() for x in sas_paths.split(',')]
        for sas_path in sas_paths:
            print(f'sas_path: {sas_path}')
            ev, rv, wv = input_filter.check_permissions(sas_path)
            if not ev or not rv:
                error.append(f'permission error in input file path {sas_path}  [code = {ev}{rv}{wv}]')
                if not ev:
                    error.append(f'sas path "{sas_path}" does not exist')
                    return error
                elif not rv:
                    error.append(f'read permission not allowed for sas path "{sas_path}"')
                    return error
            if sas_type == 0:
                suffix = ['*.iq', '*.log']
                extra = ['*.inf', '*.crd', '*.ans', '*.pr']
            elif sas_type == 1:
                suffix = ['*.iq', '*.log']
                extra = ['*.inf', '*.crd', '*.ans', '*.pr']
            elif sas_type == 2:
                suffix = ['*.int', '*.log']
                extra = ['*.sav', '*.flm', '*.alm', '*.ans']
            elif sas_type == 3:
                suffix = ['*.int', '*.log']
                extra = ['*.sav', '*.flm', '*.alm', '*.ans']
            spec_files = glob.glob(os.path.join(sas_path, suffix[0]))
            log_files = glob.glob(os.path.join(sas_path, suffix[1]))
            extra_files_list = []
            for ex in extra:
                extra_files_list.append(glob.glob(os.path.join(sas_path, suffix[0])))
            number_spec_files = len(spec_files)

            dict_sastype = {0: 'sascalc', 1: 'xtal2sas', 2: 'cryson', 3: 'crysol'}
            name_sastype = dict_sastype[sas_type]
            if number_spec_files == 0:
                error.append(f"there are no scattering files found for the selected sas-type: '{name_sastype}' in folder: '{sas_path}'")
                return error

            if extract_trajectory:
                if number_spec_files != number_of_frames:
                    error.append("number of SAS files does not match number of frames in the trajectory files")
                    return error
            else:
                number_of_frames = number_spec_files

    if extract_trajectory and extract_sas:
        filetype = 'PDB/DCD and/or SAS'
        unit = 'frame and/or SAS file'
        units = 'frames and/or SAS files'
    elif extract_trajectory:
        filetype = 'PDB/DCD'
        unit = 'frame'
        units = 'frames'
    elif extract_sas:
        filetype = 'SAS'
        unit = 'SAS file'
        units = 'SAS files'

    # checking extracting options
    if option == 'single_frame':
        try:
            this_value = locale.atoi(local_value)
            if this_value > number_of_frames:
                error.append(f'there are {number_of_frames} {units} in your data path : you requested {unit} number {this_value}')
                return error
            if this_value < 1:
                error.append(f'you entered: "{local_value}": the {unit} number to be extracted needs to be a positive integer')
                return error
        except Exception as e:
            error.append(f'you entered "{local_value}" : the number of {units} to be extracted needs to be an integer')
            return error

    elif option == 'range':
        try:
            regex_num_hyphen_num = re.compile(r'\s*[+-]?\d+\s*-\s*[+-]?\d+\s*')
            match_num_hyphen_num = regex_num_hyphen_num.findall(local_value)[0]
            if match_num_hyphen_num != local_value.strip():
                error.append(f'range needs to be two integers separated by a hyphen : you entered "{local_value}"')
                return error
            regex_num_hyphen = re.compile(r'\s*[+-]?\d+\s*-')
            match = regex_num_hyphen.match(match_num_hyphen_num)
            tr1 = locale.atoi(match_num_hyphen_num[match.start():match.end() - 1])
            tr2 = locale.atoi(match_num_hyphen_num[match.end():])
            if tr1 > tr2:
                error.append(f'range needs to be from low to higher integer : you entered "{local_value}"')
                return error
            if tr1 < 1:
                error.append(f'lower limit of the range needs to be greater than 0 : you entered "{tr1}"')
                return error
            if tr1 == tr2:
                error.append(f'lower and higher limits in the range should be different : you entered "{local_value}"')
                return error
            if tr1 > number_of_frames:
                error.append(f'lower limit of the range needs to be equal or smaller than the maximum number of {units} : you entered "{tr1}"')
                return error
            if tr2 > number_of_frames:
                error.append(f'higher limit of the range needs to be equal or smaller than the maximum number of {units} : you entered "{tr2}"')
                return error
        except Exception as e:
            error.append(f'range needs to be two integers separated by a hyphen! : you entered "{local_value}"')
            return error

    elif option == 'text_file':
        try:
            error = input_filter.check_file_exists(local_value)
            if error:
                return error
        except Exception as e:
            error.append(f'file : {local_value} does not exist')
            return error
        try:
            with open(local_value, 'r') as infile:
                data = []
                for line in infile:
                    lin = line.split()
                    if lin:
                        this_value = locale.atoi(lin[0])
                        data.append(this_value)
                data.sort()
                for i in range(len(data)):
                    if data[i] < 1:
                        error.append(f'text file can only have positive integers : "{data[i]}" was found in the text file')
                        return error
                    if data[i] > number_of_frames:
                        error.append(f'there are {number_of_frames} {units} in your data path : you requested {unit} number {data[i]} in the text file')
                        return error
                    if i > 0 and data[i] == data[i - 1]:
                        error.append(f'redundant {unit} number "{data[i]}" found in the text file')
                        return error
        except Exception as e:
            error.append(f'encountered an unknown error reading text_file: {local_value}')
            return error

    elif option == 'weight_file':
        try:
            error = input_filter.check_file_exists(local_value)
            if error:
                return error
        except Exception as e:
            error.append(f'file : {local_value} does not exist')
            return error
        try:
            with open(local_value, 'r') as infile:
                frames = []
                weights = []
                for line in infile:
                    words = line.split()
                    if not words:
                        continue
                    if words[0][0] == '#':
                        continue
                    else:
                        if locale.atof(words[2]) not in [0.0, 1.0]:
                            error.append(f'weight file column 3 can only have 0 or 1 : {words[2]} was found')
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
                    if frames[i] < 1:
                        error.append(f'weight file column 1 can only have positive integers : "{frames[i]}" was found in the weight file')
                        return error
                    if frames[i] > number_of_frames:
                        error.append(f'there are {number_of_frames} {units} in your data path : {unit} number {frames[i]} was found in the weight file')
                        return error
                    if i > 0 and frames[i] == frames[i - 1]:
                        error.append(f'redundant {unit} number "{frames[i]}" found in the weight file')
                        return error
        except Exception as e:
            error.append(f'encountered an unknown error reading weight_file: {local_value}')
            return error

    elif option == 'sampling_frequency':
        try:
            sampling_frequency = locale.atoi(local_value)
            if sampling_frequency < 1:
                error.append(f'the sampling frequency needs to be a positive number : you entered "{local_value}"')
                return error
            if sampling_frequency > number_of_frames:
                error.append(f'the sampling frequency needs to be smaller than the number of {units} in your data path : you entered "{local_value}"')
                return error
        except Exception as e:
            error.append(f'encountered an unknown error reading the sampling frequency : {local_value}')
            return error

    return error