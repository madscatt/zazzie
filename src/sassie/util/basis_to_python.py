# -*- coding: utf-8 -*-

#    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#       BASIS TO PYTHON
#
#                             initial coding         :   Joseph E. Curtis
#       8/2024        --      Python 3               :   Susan Krueger
#       8/2024        -- added support for = and !=  :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Basis to Python** is the method that converts a basis string to a python-readable string.

    **Inputs:**

    basis string with VMD-like syntax defining selection basis for each component

    **Outputs:**
    
    python-readable basis string

'''

import shlex

DEBUG = False
# DEBUG = True

# NOTE: <, >, and ! must be included individually in the ignore_list; otherwise, !=, >=, <=, and == will not be parsed correctly (SK)
ignore_list = ['>', '<', '!', '=', '!=', '==', '>=', '<=',
               'and', 'not', 'or', '(', ')', '.']
word_list = ['and', 'not', 'or']


def parse_basis(basis):
    '''
    **Parse Basis** parses the basis string with VMD-like syntax to create a python-readable basis string.

    '''

    kd = {}
    # the 'atom' keyword currently doesn't work (8/2024) (SK)
    kd['atom'] = ['atom[i] == ', 'atom[i] ']
    kd['index'] = ['index[i] == ', 'index[i] ']
    kd['name'] = ['name[i] == ', 'name[i] ']
    kd['resname'] = ['resname[i] == ', 'resname[i] ']
    kd['chain'] = ['chain[i] == ', 'chain[i] ']
    kd['resid'] = ['resid[i] == ', 'resid[i] ']
    kd['occupancy'] = ['occupancy[i] == ', 'occupancy[i] ']
    kd['beta'] = ['beta[i] == ', 'beta[i] ']
    kd['segname'] = ['segname[i] == ', 'segname[i] ']
    kd['element'] = ['element[i] == ', 'element[i] ']
    kd['charge'] = ['charge[i] == ', 'charge[i] ']
    kd['moltype'] = ['moltype[i] == ', 'moltype[i] ']
    kd['rg'] = ['rg[i] == ', 'rg[i] ']
    kd['x2'] = ['x2[i] == ', 'x2[i] ']

    if DEBUG:
        print('original_basis = ' + basis + '\n')

    lexer = shlex.shlex(basis)
    original_tokenlist = []
    for token in lexer:
        original_tokenlist.append(str(token))
    if DEBUG:
        print('original_tokenlist = ', original_tokenlist)

    tokenlist = []
    number_of_tokens = len(original_tokenlist)
    if DEBUG:
        print('original_tokenlist, len(original_tokenlist) = ',
              original_tokenlist, len(original_tokenlist))

    i = 0
    while (i < number_of_tokens):
        this_token = original_tokenlist[i]
        if DEBUG:
            print('i, this_token = ', i, this_token)
        if DEBUG:
            print('this_token = ', this_token),
        if i < number_of_tokens-1:
            next_token = original_tokenlist[i+1]
            if DEBUG:
                print('next_token = ', next_token)
            if (this_token in ignore_list):
                if DEBUG:
                    print('this_token is in ignore_list')
                if (next_token in ignore_list):
                    if DEBUG:
                        print('next_token is in ignore_list')
                    if (this_token in word_list or next_token in word_list):
                        if DEBUG:
                            print('this_token or next_token is in word_list')
                        tokenlist.append(this_token+' '+next_token)
                    else:
                        if DEBUG:
                            print('this_token and next_token are not in word_list')
                        tokenlist.append(this_token+next_token)
                    i += 1
                elif (this_token == '='):
                    if DEBUG:
                        print('this_token is =')
#                        print('this_token: ', this_token)
#                        print('tokenlist: ', tokenlist)
                    this_token = this_token.replace('=', '==')
                    if DEBUG:
                        print('this_token: ', this_token)
                    tokenlist.append(this_token)
                    if DEBUG:
                        print('tokenlist after replace: ', tokenlist)
# NOTE: #the line below had to be commented out in order for both "=" and "==" to be parsed correctly (SK)
#                    i += 1
                else:
                    if DEBUG:
                        print('this_token is not =')
#                        print('this_token: ', this_token)
                    tokenlist.append(this_token)

            # elif(this_token not in kd and not this_token.replace(".", "", 1).isdigit()):
            #    tokenlist.append('"'+this_token+'"')

            elif (this_token not in kd and not this_token.isdigit()):
                # added this if statement to handle the case where the token is in "" but is not the last token in the list
                # example: not name "CA" and resid > 2 (SK)
                if (this_token[0] != '"' and this_token not in ignore_list):
                    tokenlist.append('"'+this_token+'"')
                else:
                    tokenlist.append(this_token)
                # end of addition (SK)
            else:
                tokenlist.append(this_token)
            # print('this_token = ',this_token)
        else:
            if (this_token not in kd and not this_token.isdigit()):

                if (this_token[0] != '"' and this_token not in ignore_list):
                    tokenlist.append('"'+this_token+'"')
                else:
                    tokenlist.append(this_token)
            else:
                tokenlist.append(this_token)
        i += 1

    if DEBUG:
        print('tokenlist = ', tokenlist)

    number_of_tokens = len(tokenlist)
    new_basis = ''

    for i in range(number_of_tokens):
        this_word = tokenlist[i]
        try:
            next_word = tokenlist[i+1]
            nwe = True
        except:
            nwe = False

        if this_word in kd:
            if DEBUG:
                print('this_word,next_word = ', this_word, next_word)
            if nwe:
                if next_word not in ignore_list:
                    new_basis += kd[this_word][0]
                else:
                    new_basis += kd[this_word][1]

        else:
            new_basis += ' '+this_word+' '

    if DEBUG:
        print('new_basis = ', new_basis, '\n')

    return new_basis


def concatenate_decimal(string):
    '''
    **Concatenate Decimal** concatenates the decimal numbers in the basis string that are separated by a space and a decimal point.

    '''

    import re
    if DEBUG:
        print('in concatenate_decimal')
        print('input string = ', string)
    find_decimals = re.findall(r"\d+\s+.\s+\d+", string)
    if DEBUG:
        print('find decimals = ', find_decimals)
    for m in find_decimals:
        new_string = m.replace(' ', '')
        if DEBUG:
            print('previous string, new string = ', m, new_string)
        string = string.replace(m, new_string)
        if DEBUG:
            print('returned string = ', string)
    return string


def clean_up_weight_basis(basis_string):
    '''
    **Clean Up Weight Basis** cleans up the basis strings that define the weight files in chi-square filter.
    The syntax of this basis string is different from that used for Sascalc, TAMC, Build Utilities, etc.

    '''

#    import sasconfig as sasconfig
#    if sasconfig.__cluster__ == True:
#        new_basis_string = string.split(basis_string.encode('ascii','ignore'),',')
#        new_basis_string = string.split(str(basis_string),',')
#    else:
#        new_basis_string = string.split(basis_string,',')

    new_basis_string = str(basis_string).split(',')

    python_basis = []

    for basis in new_basis_string:
        # the following line is no longer needed due to changes in parse_basis to handle ==, !=, >=, <= (SK)
        #        basis = basis.replace('=', '==')
        this_python_basis = parse_basis(basis)
        this_python_basis = this_python_basis.replace('"', '')
        this_python_basis = concatenate_decimal(this_python_basis)

        python_basis.append(this_python_basis)

    return python_basis


if __name__ == '__main__':
    import sys

    basis = []
    basis.append('(name CA and name NH) or resid > 43')
    basis.append('(name CA and name NH) or resid > 43 and resid < 57')
    basis.append('segname HC1 and (resid >= 210 and resid <=214)')
    basis.append('segname TEST and resid < 210')
    basis.append('resname TIP3')
    basis.append('resname GLY and index != 1')
    basis.append('index = 523')
    # "1+": the quotes are needed for the entry to be recognized as a string
    basis.append('charge "1+"')
    basis.append('name "CA" and resid > 43')
    basis.append('not name "CA" and resid > 43')
    basis.append('resid = 10')
    basis.append('resid == 10')
    basis.append('resid != 10')
    basis.append('moltype protein')
    basis.append('not moltype nucleic')
    basis.append('occupancy < "1.00"')
    basis.append('occupancy <= "1.00"')
    basis.append('occupancy = "1.00"')
    # the following results in occupancy == "1.00" as well, but the basis_to_python_filter requires the = sign, so there won't be mix ups when using !=, >=, <=, ==, etc. (SK)
    basis.append('occupancy "1.00"')

    for i in range(len(basis)):
        print('#####')
        print('basis ' + str(i) + ' = ' + str(basis[i]))
        python_basis = parse_basis(basis[i])
        print('python_basis = ', python_basis)
        print('\n', '\n')


#####
    # The following is for testing whether the right basis is obtained from a PDB file and written to a new PDB file

    import sasmol.system as system
    m = system.Molecule(0)
    m.read_pdb('min3.pdb')

#    basis = 'not name "CA" and (resid >= 23 and resid <= 68) and resid != 43'
# manually input test basis string
    basis = input('test basis: ')
    print('test basis = ', basis)
    python_basis = parse_basis(basis)
    print('python_basis = ', python_basis)

    sub_mol = system.Molecule(0)

    frame = 0
    error, mask = m.get_subset_mask(python_basis)

    if (len(error) > 0):
        print('error getting mask = ', error)
    else:
        print('mask = ', mask)

    import numpy
    print(numpy.sum(mask))

    if (numpy.sum(mask) > 0):
        error = m.copy_molecule_using_mask(sub_mol, mask, frame)

    if (len(error) > 0):
        print('error copying molecule using mask = ', error)
    else:
        sub_mol.write_pdb('test.pdb', frame, 'w')
        print('new PDB file written')

    sys.exit()
#####

    # The following is for testing the weight basis in chi_square filter
    # The syntax of this basis is different from that used for Sascalc, TAMC and Build Utilities

    rg = [float(x) for x in range(20)]
    x2 = [0.1*float(x) for x in range(20)]

#    basis_string = 'rg != 0.5'
#    basis_string = 'rg = 0.5 or rg = 1.2, rg < 2.5'
    basis_string = 'rg >= 0.5, rg <= 2.5'

    print('\n\nbasis_string = ', basis_string)
    python_basis = clean_up_weight_basis(basis_string)
    print('python_basis = ', python_basis)

    for basis in python_basis:
        print('basis = ', basis)
        mask_array = []

        try:
            for i in range(len(rg)):
                if (eval(basis)):
                    mask_array.append(1)
                else:
                    mask_array.append(0)
        except:
            print('ERROR: failed to parse basis')
            sys.exit()

        print('mask_arrary = ', mask_array)
