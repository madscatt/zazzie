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
import os
import sys
import string
import locale
import glob
import sassie.interface.input_filter as input_filter


def parse_grammar(basis):

    import re

#    print 'basis in parse_grammar: ', basis

    error = []

    p = re.compile(
        '\s*(Rg|RG|X2)\s*(>|<|>=|<=)\s*(\'[\d\s\.]+\'|"[\d\s\.]+"|[\d\.]+)\s*')
    if p.search(basis):
        error.append('expression "' + basis + '" not understood!')
        error.append('you may need to change "Rg|RG|X2" to lower cases')
        return error

    p = re.compile(
        '\s*(rg|x2)\s*(>|<|>=|<=)\s*(\'[\d\s\.]+\'|"[\d\s\.]+"|[\d\.]+)\s*')
    list_indices = []
    for m in p.finditer(basis):
        list_indices.append([m.start(), m.end()])

#    print 'list_indices: ',list_indices

    if not len(list_indices):
        error.append('grammar wrong in expression "' + basis + '"')
        return error

    words = []
    words.append(basis[0:list_indices[0][0]].strip())
    for i in range(1, len(list_indices)):
        words.append(basis[list_indices[i - 1][1]: list_indices[i][0]].strip())
    words.append(basis[list_indices[-1][1]:].strip())

#    print 'words: ',words

    if len(words[0]) and words[0][0] != '(' or len(words[-1]) and words[-1][-1] != ')':
        error.append('grammar wrong in expression "' + basis + '"')
        return error

    for word in words[1:-1]:
        if word[0] == '(' or word[-1] == ')':
            # print word
            error.append('grammar is wrong in expression "' + basis + '"')
            return error

    r = re.compile('(\(\s*\)|\(\s*(and|or)|(and|or)\s*(and|or)|(and|or)\s*\))')
    for word in words:
        if r.search(word):
            error.append('grammar wrong in expression "' + basis + '"')
            return error

    complete_words = ''.join(words)
    if complete_words.replace('and', '').replace('or', '').replace('(', '').replace(')', '').strip():
        error.append('expression "' + basis + '" not understood!')
        return error
    if complete_words.count('(') != complete_words.count(')'):
        error.append('parenthesis in expression "' +
                     basis + '" does not match!')
        return error

    return error


def check_chi_square_filter(variables, **kwargs):

    runname = variables['runname'][0]
    saspaths = variables['saspaths'][0]
    sasintfiles = variables['sasintfiles'][0]
    io = variables['io'][0]
    basis_string = variables['basis_string'][0]
    sastype = variables['sastype'][0]
    plotflag = variables['plotflag'][0]
    reduced_x2 = variables['reduced_x2'][0]
    number_of_weight_files = variables['number_of_weight_files'][0]
    plotflag = variables['plotflag'][0]
    reduced_x2 = variables['reduced_x2'][0]
    weight_file_names = variables['weight_file_names'][0]

    error = []
    error = input_filter.check_name(runname)
    if(error != []):
        return error

#   saspaths existence and permissions is checked later
#   this generic path check is now consistent with other module filters, but...
#   path variable isn't used here or in chi_square_filter module
#    if "no_file_check" not in kwargs:
#
#        ev, rv, wv = input_filter.check_permissions(path)
#
#        if(not ev or not rv):
#            error.append('permission error in input file path ' +
#                        path + '  [code = ' + str(ev) + str(rv) + str(wv) + ']')
#           if(ev == False):
#                error.append('path does not exist')
#            elif(rv == False):
#                error.append('read permission not allowed')
#            return error

    sas_paths = [x.strip() for x in saspaths.split(',')]
    sasintfiles = [x.strip() for x in sasintfiles.split(',')]

    if(len(sas_paths) != len(sasintfiles)):
        error.append(
            'number of interpolated data files does not match number of SAS data paths')
        return error

    Q_exps = []  # to keep track of Q values for all sasintfiles in case they are different
    for sasintfile in sasintfiles:
        error = input_filter.check_file_exists(sasintfile)
        if(len(error) != 0):
            error.append('check interpolated data file path + filename : ' +
                         os.path.basename(sasintfile))
            return error

        else:
            try:
                lines = open(sasintfile, 'r').readlines()
            except:
                error.append("unable to open and read your interpolated data file : " +
                             os.path.basename(sasintfile))
                return error

#           check interpolated sasintfile for three columns (q,I(q),error)
            sum_I = 0.0
            Q_exp = []
            for line in lines:
                words = line.split()
                if len(words) == 0 or words[0][0] == '#':
                    continue
                elif len(words) > 3 and words[3][0] != '#':
                    error.append('the number of columns should be 3 in your interpolated data file : ' +
                                 os.path.basename(sasintfile) + '  \n' + line)
                    return error
                elif len(words) < 3:
                    error.append('the number of columns should be 3 in your interpolated data file : ' +
                                 os.path.basename(sasintfile) + '  \n' + line)
                    return error
                elif len(words) == 3 and words[2][0] == '#':
                    error.append('the number of columns should be 3 in your interpolated data file : ' +
                                 os.path.basename(sasintfile) + '  \n' + line)
                    return error
                else:
                    try:
                        q = locale.atof(words[0])
                        I = locale.atof(words[1])
                        Error = locale.atof(words[2])
                        Q_exp.append(q)
                    except:
                        #    error.append('the following line cannot be parsed in your SAS data file : '+os.path.basename(sasintfile)+'  \n'+line)
                        #    return error
                        pass
                    else:
                        if q < 0.0:
                            error.append('negative Q value found in your SAS data file : ' +
                                         os.path.basename(sasintfile) + '  \n' + line)
                            return error
                        # elif I<0.0:
                        #    error.append('negative I value found in your SAS data file : '+os.path.basename(sasintfile)+'  \n'+line)
                        #    return error
                        # elif Error<0.0:
                        #    error.append('negative error value found in your SAS data file : '+os.path.basename(sasintfile)+'  \n'+line)
                        #    return error
                        else:
                            sum_I += abs(I)
            if sum_I == 0.0:
                error.append(
                    "no lines or all-zero intensity entries found in your SAS data file : " + os.path.basename(sasintfile))
                return error
        Q_exps.append(Q_exp)

#    print 'Q_exps: ', Q_exps

    if(plotflag != 0 and plotflag != 1):
        error.append(
            'plot flag needs to be 1 or 0 ... you entered: ' + str(plotflag))
        return error

    # check basis_string
    list_basis = basis_string.split(',')
#    print 'list_basis in filter: ', list_basis
#    print 'length of list_basis: ',len(list_basis)
#    print 'number of weight_files: ', number_of_weight_files
    if((len(list_basis) != number_of_weight_files) and (number_of_weight_files > 0)):
        error.append(
            'number of weight files does not match length of basis string')
        return error
    if(number_of_weight_files > 0):
        try:
            for basis in list_basis:
                error = parse_grammar(basis)
                if len(error) > 0:
                    print 'error in parse_grammar: ', error
                    return error
        except:
            pass

    # check weight_file_names
    list_names = weight_file_names.split(',')
#    print 'list_names in filter: ', list_names
    if((len(list_names) != number_of_weight_files) and (number_of_weight_files > 0)):
        error.append(
            'number of weight files does not match number of weight file names')
        return error

    # check sastype
#    print 'sastype: ', sastype
    if sastype not in [0, 1, 2, 3]:
        error.append("sastype %d not supported!" % sastype)
        return error

    # check saspaths
#    print 'sas_paths in filter: ',sas_paths

    count = 0
    for sas_path in sas_paths:
#        print 'sas_path in filter: ', sas_path
#        print 'count: ', count
        ev, rv, wv = input_filter.check_permissions(sas_path)
        if(not ev or not rv):
            error.append('permission error in input file path ' +
                         sas_path + '  [code = ' + str(ev) + str(rv) + str(wv) + ']')
            if(ev == False):
                error.append('sas path "' + sas_path + '" does not exist')
                return error
            elif(rv == False):
                error.append(
                    'read permission not allowed for sas path "' + sas_path + '"')
                return error
        if(sastype == 0 or sastype == 1):
            suffix = ['*.iq', '*.log']
            extra = ['*.inf', '*.crd', '*.ans', '*.pr']
        elif(sastype == 2):
            suffix = ['*.int', '*.log']
            extra = ['*.sav', '*.flm', '*.alm', '*.ans']
        elif(sastype == 3):
            suffix = ['*.int', '*.log']
            extra = ['*.sav', '*.flm', '*.alm', '*.ans']
        spec_files = glob.glob(os.path.join(sas_path, suffix[0]))
        log_files = glob.glob(os.path.join(sas_path, suffix[1]))

        number_spec_files = len(spec_files)
        dict_sastype = {0: 'sascalc', 1: 'xtal2sas', 2: 'cryson', 3: 'crysol'}
        name_sastype = dict_sastype[sastype]
        if number_spec_files == 0:
            error.append("there are no scattering files found for the selected sas-type: '" +
                         name_sastype + "' in folder: '" + sas_path + "'")
            return error

        Q_spec = []
        for line in open(spec_files[0]).readlines():
            words = line.split()
            if len(words) == 0 or words[0][0] == '#':
                continue
            else:
                try:
                    q = locale.atof(words[0])
                    I = locale.atof(words[1])
                except:
                    #error.append('the following line cannot be parsed in your SAS file : '+os.path.basename(spec_files[0])+'  \n'+line)
                    # return error
                    pass
                else:
                    Q_spec.append(q)
#        print 'Q_spec: ',Q_spec
#        print 'Q_exps[count]: ',Q_exps[count]
        Q_exp = Q_exps[count]
#        print 'Q_exp: ',Q_exp
        if Q_spec != Q_exp:
            if len(Q_spec) != len(Q_exp):
                error.append(
                    'The number of Q values in supplied interpolated data file does not match Q value in theoretical SAS profiles\n')
                return error
            else:
                for i in xrange(min(len(Q_spec), len(Q_exp))):
                    if (Q_spec[i] != Q_exp[i]):
                        diff_start = i
                        break
                if diff_start:
                    output_Q_spec, output_Q_exp = '...', '...'
                else:
                    output_Q_spec, output_Q_exp = '', ''
                Q_spec = Q_spec[diff_start:]
                Q_exp = Q_exp[diff_start:]
                output_Q_spec += str(Q_spec[:min(10, len(Q_spec))])
                if len(Q_spec) > 10:
                    output_Q_spec += '...'
                output_Q_exp += str(Q_exp[:min(10, len(Q_exp))])
                if len(Q_exp) > 10:
                    output_Q_exp += '...'
                error.append(
                    'Q value or spacing in supplied interpolated data file does not match Q value or spacing in theoretical SAS profiles\n')
                error.append('Q value in your SAS file : ' + os.path.basename(spec_files[
                    0]) + ':\n' + output_Q_spec + '\n  does not match that in your experimental file:\n' + output_Q_exp)
                return error

        # check file size in the SAS folder
        file_size = os.path.getsize(spec_files[0])

        q_in_spec_file_0 = []
        spec_file_0 = open(spec_files[0]).readlines()

        for line in spec_file_0:
            lin = string.split(line)
            q_in_spec_file_0.append(lin[0])

        if len(spec_files) > 1:
            for f in spec_files[1:]:
                q_in_spec_file = []
                this_spec_file = open(f).readlines()
                for line in this_spec_file:
                    lin = string.split(line)
                    q_in_spec_file.append(lin[0])
                if os.path.getsize(f) != file_size:
                    if len(q_in_spec_file_0) != len(q_in_spec_file):
                        error.append('number of q-values of "' + f + '" does not match "' + os.path.basename(
                            spec_files[0]) + '".\nPlease verify the files in the SAS path are consistent!')
                        return error
                if q_in_spec_file_0 != q_in_spec_file:
                    error.append('q-values of "' + f + '" do not match "' + os.path.basename(
                        spec_files[0]) + '".\nPlease verify the files in the SAS path are consistent!')
                    return error

        number_log_files = len(log_files)

        if number_spec_files != number_log_files:
            error.append(
                "the number of log files does not match the number of theoretical SAS profiles in sas path: '" + sas_path + '"')

        count += 1

    # check io
#    print len(io)
    if(len(sas_paths) != len(io)):
        error.append(
            'number I(0) values does not match number of of SAS data paths')
        return error
    for i in range(len(io)):
        if io[i] <= 0.0:
            error.append('I(0): "' + str(io[i]) + '" should be greater than 0')
            return error

    # return
    return error
