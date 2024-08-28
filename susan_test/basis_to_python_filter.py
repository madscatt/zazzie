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
#       BASIS TO PYTHON FILTER
#
#       7/2023        --      Initial coding         :   Susan Krueger
#       8/2024        --      Python 3               :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Basis to Python Filter** is the method that checks the syntax of a basis
    string to make it more likely to work after it is converted to a
    python-readable string.

    **Inputs:**

        Basis string

    **Outputs:**

        Error string

    Called from modules that use **Basis to Python** to convert a basis string
    with VMD-like syntax to a python-readable string.

'''
import re

DEBUG = False


def check_basis_syntax(basis):
    """
    Method to check that a basis string has the correct syntax for conversion to a python-readable string.

    Note
    ----

    Allowed keywords (must be in lower case):

    **atom:** accepts ATOM or HETATM; not currently implemented (8/2024)

    **segname:** accepts alphnumeric (case-sensitive, no more than 4 characters); double quotes or no quotes are OK

    **name, resname:**  accept A-Z (uc) and 0-9; no more than 4 characters, double quotes or no quotes are OK (name "CA" or name CA or name HT1; resname "GLY" or resname TIP3)

    **moltype:** accepts protein or nucleic (lc, double quotes or no quotes OK); nucleic will find both dna and rna, since there are ambiguities with dna and rna base names; DO NOT use to distinguish between dna and rna; instead use segname, chain, etc.

    **element:** accepts A-Z (uc) and optional a-z (lc)

    **charge:** accepts 0-9 and optional + or -; note that the charge parameter MUST have double quotes to be recognized as a string (charge "1+", charge '1-' or charge "1")

    **chain:** accepts lc, uc and numbers 0-9 (single letter or digit); double quotes or no quotes are OK

    **index, resid:** accept integer; must use >, <, >=, <=, ==, =, !=; no quotes

    **beta, occupancy:** accept string; must use >, <, >=, <=, ==, =, !=; MUST have double quotes and two decimal places (beta = "0.00", occupancy = "1.00")

    Parameters
    ----------

        basis:  string
            The basis string with VMD-like syntax used to select a portion of a PDB file

    Returns
    -------

        error: string
            The error message generated when a check fails. If there are no failures, the error is blank.

    """

    if DEBUG:
        # added the hash marks to see if there are any leading or trailing spaces
        print('basis in check_basis_syntax: #' + basis + '#')
    error = []

# check to see if the keyword is correct, but is not in all lower case
    caps_check = re.compile(
        r'(NAME|Name|RESNAME|Resname)\s*(\"[A-Z\d]{1,4}"|[A-Z\d]{1,4})\s*|(Element|ELEMENT)\s*(\"[A-Z][a-z]?"|[A-Z][a-z]?)\s*|(Chain|CHAIN)\s*(\"[\w]"|[\w])\s*|(Segname|SEGNAME)\s*(\"[\w]{1,4}"|[\w]{1,4})\s*|(Charge|CHARGE)\s*(\"(d+\-?|\d+\+?)")\s*|(Moltype|MOLTYPE)\s*(\"?(protein|nucleic)\"?)\s*|(Index|INDEX|Resid|RESID)\s*(>|<|>=|<=|==|=|!=)\s*([\d]+)\s*|(Beta|BETA|Occupancy|OCCUPANCY)\s*(>|<|>=|<=|==|=|!=)\s*(\"[\d]+\.[\d]{2}\")\s*')

    if caps_check.search(basis):
        error.append('Expression "' + basis + '" not understood!')
        error.append('you may need to change "Index|INDEX|Resid|RESID|Beta|BETA|Occupancy|OCCUPANCY|Segname|SEGNAME|Chain|CHAIN|Name|NAME|Resname|RESNAME|Element|ELEMENT|Charge|CHARGE|Moltype|MOLTYPE" to all lower case')
        print('error: ', error)
        return error

# check to see if the basis contains at least one of the allowed (lc) keywords
    syntax_check = re.compile(
        r'(name|resname)\s*(\"[A-Z\d]{1,4}"|[A-Z\d]{1,4})\s*|element\s*(\"[A-Z][a-z]?"|[A-Z][a-z]?)\s*|chain\s*(\"[\w]"|[\w])\s*|segname\s*(\"[\w]{1,4}"|[\w]{1,4})\s*|charge\s*(\"(d+\-?|\d+\+?)")\s*|moltype\s*(\"?(protein|nucleic)\"?)\s*|(index|resid)\s*(>|<|>=|<=|==|=|!=)\s*([\d]+)\s*|(beta|occupancy)\s*(>|<|>=|<=|==|=|!=)\s*(\"[\d]+\.[\d]{2}\")\s*')


# if a keyword is found, then check to see if the syntax of the entire expression is correct
    list_indices = []
    for n in syntax_check.finditer(basis):
        list_indices.append([n.start(), n.end()])

    if DEBUG:
        print('list_indices, len(list_indices): ',
              list_indices, len(list_indices))
# example segment TEST
# example resid 10
    if not len(list_indices):
        error.append('(0) syntax wrong in expression "' + basis +
                     '". No keywords found or incorrect use of keyword arguments.')
#        print('error: ', error)
        return error

# split the expression into individual words to look for 'and', 'or', 'not'
    words = []
    words.append(basis[0:list_indices[0][0]].strip())
    for i in range(1, len(list_indices)):
        words.append(basis[list_indices[i - 1][1]: list_indices[i][0]].strip())
    words.append(basis[list_indices[-1][1]:].strip())

    if DEBUG:
        print('words: ', words)
        print('len(words[0]): ', len(words[0]))

# check to make sure that "and or "or" is used if there are multiple keywords
# example segname SEG1 resname GLY
# example (segname SEG1) resname GLY
# example segname SEG1 not resname GLY
    patterns = [
        r'\band\b',
        r'\bor\b',
        r'\bnot\b',
        r'\band\s+not\b',
        r'\bor\s+not\b',
        r'\band\s*\(not\b',
        r'\band\s+not\s*\(',
        r'\bor\s*\(not\b',
        r'\bor\s+not\s*\(',
        r'\bor\s\('
    ]

# Compile the patterns into a single regex
    combined_pattern = re.compile('|'.join(patterns))
    if len(list_indices) > 1:
        found_connectors = [
            word for word in words if combined_pattern.search(word)]
        if DEBUG:
            print('found, leng: ', found_connectors, len(found_connectors))
# example resname A (segname test)
# example resname A name B
        if (len(found_connectors) == 0):
            error.append('(1) syntax wrong in expression "' + basis +
                         '". You must use "and" or "or" to separate expressions.')
            return error
# check if 'not' is used by itself in cases where there is more than one keyword
# example resname A not resname B
        if ('not' in found_connectors):
            error.append('(1) syntax wrong in expression "' + basis +
                         '". You must use "and" or "or" to separate expressions.')
#            print('error: ', error)
            return error

# make sure 'and', 'or' is not used if there is only one keyword and that 'not' is used correctly
# example but resname A
# example resname A and
    if len(words[0]) and words[0][0] != '(' and len(words[0]) and words[0][0] != 'n' or len(words[-1]) and words[-1][-1] != ')':
        if DEBUG:
            if len(words[0]):
                print('words[0][0]: ', words[0][0])
            if len(words[-1]):
                print('words[-1][-1]: ', words[-1][-1])
        error.append('(2) syntax wrong in expression "' +
                     basis + '". Possible incorrect use of "and", "or" or "not".')
#        print('error: ', error)
        return error

# If syntax is OK, check the entire expression to make sure that parentheses, and, or are used correctly
    connector_check = re.compile(
        '(\(\s*\)|\(\s*(and|or)|(and|or)\s*(and|or)|(and|or)\s*\))')
# example: resname A ( and not ) resname B
# example segname test (and name CA)
# example (segname test and) name CA
# example resname A and or name CA
# example resname A ()
    for word in words:
        if DEBUG:
            print('word1: ', word)
        if connector_check.search(word):
            error.append('(3) syntax wrong in expression "' + basis + '"')
#            print('error: ', error)
            return error

# test for wrong words in an expression containing allowed keywords
    complete_words = ''.join(words)
    if DEBUG:
        print('complete words: ', complete_words)
        test = complete_words.replace('and', '').replace('or', '').replace(
            'not', '').replace('(', '').replace(')', '').strip()
        print('test: ', test)
# example n0t name CA
    if complete_words.replace('and', '').replace('or', '').replace('not', '').replace('(', '').replace(')', '').strip():
        error.append('(4) expression "' + basis + '" not understood!')
#        print('error: ', error)
        return error

# test to make sure parentheses match
# example (resname A
    if complete_words.count('(') != complete_words.count(')'):
        error.append('(5) parenthesis in expression "' +
                     basis + '" do not match!')
#        print('error: ', error)
        return error

    return error


if __name__ == '__main__':

    import sys
    import basis_to_python as basis_to_python
    basis = input('test expression: ')
    print('basis: ', basis)
    error = check_basis_syntax(basis)
#    print('error, len(error): ', error, len(error))
    if len(error):
        print('error: ', error)
        sys.exit()
    else:
        try:
            python_basis = basis_to_python.parse_basis(basis)
            print('python_basis: ', python_basis)
        except:
            error.append(
                'unable to convert input string ' + basis + ' to python readable string')
            print('error: ', error)
            sys.exit()

    # The following is for testing whether the right basis is obtained from a PDB file and written to a new PDB file

    import sasmol.system as system
    m = system.Molecule(0)
    m.read_pdb('min3.pdb')
#    m.read_pdb('pai_vn_start.pdb')
#    m.read_pdb('pfvu5_fit.pdb')

    sub_mol = system.Molecule(0)

    frame = 0
    error, mask = m.get_subset_mask(python_basis)

    if (len(error) > 0):
        print('error = ', error)

    import numpy
    print(numpy.sum(mask))

    if (numpy.sum(mask) > 0):
        error = m.copy_molecule_using_mask(sub_mol, mask, frame)

    if (len(error) > 0):
        print('error = ', error)
    else:
        sub_mol.write_pdb('test_basis_syntax.pdb', frame, 'w')
        print('new PDB file written')
