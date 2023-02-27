# -*- coding: utf-8 -*-
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

from __future__ import print_function

import os
import glob
import readline
import itertools
from operator import itemgetter
import subprocess

import sassie.util.sasconfig as sasconfig

residue_dictionary = {
    'ALA': 'A',
    'ARG': 'R',
    'ASP': 'D',
    'ASN': 'N',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HSD': 'H',
    'HIS': 'H',
    'HSE': 'H',
    'HSP': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'GUA': 'G',
    'CYT': 'C',
    'ADE': 'A',
    'THY': 'T',
    'URA': 'U',
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U',
    'DA': 'A',
    'DC': 'C',
    'DG': 'G',
    'DT': 'T'
}

nucleic_res_charmm_dict = {
    'A': 'ADE',
    'C': 'CYT',
    'G': 'GUA',
    'T': 'THY',
    'U': 'URA',
    'DA': 'ADE',
    'DC': 'CYT',
    'DG': 'GUA',
    'DT': 'THY'
}

#    'PCA': 'E', # Pyroglutamic acid
#    'MLY': 'K', # Methyllysine
#    'MSE': 'M', # Selenomethionine
#    'UNK': 'X', # Unknown


def conv_aa3to1(aa3):
    """
    Convert three letter amino acid codes to single letter ones.
    Unknown codes set to 'X'

    @type  aa3  :  string
    @param aa3  :  Three letter amino acid code
    @rtype      :  string
    @return     :  Single letter amino acid code

    """

    try:
        aa1 = residue_dictionary[aa3]
    except:
        aa1 = 'X'

    return aa1


def conv_aa1to3(aa1):

    map = {v: k for k, v in residue_dictionary.iteritems()}

    try:
        aa3 = map[aa1]
    except:
        aa3 = 'UNK'

    return aa3


def uniquify_list(lst):
    """
    Return only unique members of input list. Retains ordering unlike simple 
    use of set(lst).

    @type  lst:   list
    @param lst:   List to be filtered
    @rtype:       list
    @return:      List of unique elements from the input list

    """

    seen = set()
    seen_add = seen.add

    u_lst = [x for x in lst if not (x in seen or seen_add(x))]

    return u_lst


def get_grouped_ranges(input_data):
    """
    Convert list of digits to range strings so that [1,2,3,6,8] becomes
    ['1-3','6-8'].

    @type  input_data:  list
    @param input_data:  List of parsed PDB lines (dictionaries)
    @rtype:             list
    @return:            List with grouped numeric ranges

    """

    str_list = []

    groups = group_ranges(input_data)

    for group in groups:
        if group[0] != group[-1]:
            str_list.append('{0:d}-{1:d}'.format(group[0], group[-1]))
        else:
            str_list.append(str(group[0]))

    return str_list


def group_ranges(input_data):
    """
    Group list of digits to give first and last residues of the range as a 
    tuple, so [1,2,3,6,8] becomes [(1,3), (6,8)].

    @type  input_data:  list
    @param input_data:  List of parsed PDB lines (dictionaries)
    @rtype:             list
    @return:            List of tuples: (start_no, end_no)

    """

    data = uniquify_list(sorted(input_data))

    end_points = []

    # Use difference of index in list to value to
    for ndx, group in itertools.groupby(
            enumerate(data), lambda i_x: i_x[0] - i_x[1]):
        lst = map(itemgetter(1), group)
        end_points.append((lst[0], lst[-1]))

    return end_points


def group_sequential_nums(input_data):
    """
    Group list of digits to give first and last residues of the range as a
    tuple, so [1,2,3,6,8] becomes [[1,2,3], [6,8]].

    @type  input_data:  list
    @param input_data:  List of parsed PDB lines (dictionaries)
    @rtype:             list
    @return:            List of tuples: (start_no, end_no)

    """

    data = uniquify_list(sorted(input_data))

    ranges = []

    # Use difference of index in list to value to
    for ndx, group in itertools.groupby(
            enumerate(data), lambda i_x: i_x[0] - i_x[1]):
        lst = map(itemgetter(1), group)
        ranges.append(lst)

    return ranges


def is_sequential_numbers(numbers):
    """
    Check if input list of numbers is sequential.

    @type numbers  :  list
    @param  numbers:  List of numbers
    @rtype         :  boolean
    @return:       :  Is the input numerical list sequential
    """

    return numbers == range(numbers[0], numbers[-1] + 1)


def convert_mkd_html(mkd_filepath):
    """
    Run Pandoc to convert input markdown file to HTML

    @type mkd_filepath  :  string
    @param  mkd_filepath:  Path to the input markdown format text file
    @rtype              :  string
    @return:            :  Path to the HTML format output file
    """

    pandoc_exe = os.path.join(sasconfig.__bin_path__, 'pandoc')

    html_base = os.path.splitext(mkd_filepath)[0]
    html_filepath = html_base + '.html'

    options = [pandoc_exe, '-o', html_filepath, mkd_filepath]
    subprocess.check_call(options)

    return html_filepath


def cmdline_path_complete(text, state):
    """
    Return the state-th completion of text from the filepaths listed
    by glob.glob.

    @type text  :  string
    @param text :  Text of path to be completed
    @type state :  integer
    @param state:  Option number of completion to return
    @rtype      :  string
    @return     :  Completed path
    """
    return (glob.glob(text + '*') + [None])[state]


def get_command_line_filepath(prompt):
    """
    Get filepath from user on commadline - with tab completion.

    @type prompt :  string
    @param prompt:  Text of path to be completed
    @rtype       :  string
    @return      :  Filepath
    """

    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")
    readline.set_completer(cmdline_path_complete)

    filepath = raw_input(prompt)

    return filepath
