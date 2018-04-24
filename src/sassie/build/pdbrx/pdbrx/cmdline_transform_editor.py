#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Command line editor of BIOMT records

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

"""

from __future__ import print_function

import sys
import json
from copy import deepcopy

import numpy as np
import yaml

# make Python 2.x input behave as in Python 3
try:
    input = raw_input
except NameError:
    pass

def check_array_input(input_txt, dimensions):
    '''
    Check that input text can be converted into an array.

    @type input_txt :  str
    @param input_txt:  Text supposedly containing an array
    @type dimensions :  tuple
    @param dimensions:  Expected dimensions of array
    @rtype :  bool, np.array
    @return:  Input validity flag
              Text converted to array
    '''

    flag = True

    try:
        matrix = np.array(input_txt)
        matrix = 1.0 * matrix
        if matrix.shape != dimensions:
            flag = False

    except TypeError:
        matrix = None

    return flag, matrix


def check_rotation(input_rot):
    '''
    Check that input rotation matrix is valid (i.e. text can be converted to
    a 3*3 array)

    @type input_rot :  str
    @param input_rot:  Text format rotation matrix
    @rtype :  bool, np.array
    @return:  Input validity flag
              Text converted to array
    '''

    flag, matrix = check_array_input(input_rot, (3, 3))

    return flag, matrix


def check_translation(input_trans):
    '''
    Check that input translation vector is valid (i.e. text can be converted to
    a 3*1 array)

    @type input_trans :  str
    @param input_trans:  Text format translation vector
    @rtype :  bool, np.array
    @return:  Input validity flag
              Text converted to array
    '''

    flag, matrix = check_array_input(input_trans, (3,))

    return flag, matrix


def get_user_transform():
    '''
    Get BIOMT style rotation matrix and translation vector from user.

    @return:
    '''

    flag = True
    rot = None
    trans = None

    print("Enter rotation matrix (3 x 3, using Python list notation "
          "i.e.: [[1,0,0],[0,1,0],[0,0,1]]):")
    rot_user = eval(input())
    check_rot, rot = check_rotation(rot_user)

    if not check_rot:
        flag = False
        print("Rotation must be a 3 x 3 array of numeric values")

    else:

        print("Enter translation vector (3 x 1, using Python list notation "
              "i.e.: [1,0,0]):")
        trans_user = eval(input())
        check_trans, trans = check_translation(trans_user)

        if not check_trans:
            flag = False
            print("Error: Translation must be a 3 x 1 array of numeric values")

    if ((rot == np.identity(3)).all()) and ((trans == np.array([0.0, 0.0, 0.0])).all()):
        print("A second identity transform will not be added")
        flag = False

    return flag, rot, trans


def init_biomt(subdivs):
    '''
    Create a blank BIOMT record to be applied to the chosen subdivisions
    (chains or segments).

    @type subdivs :  list
    @param subdivs:  Segment or chain identifiers
    @rtype :  dict
    @return:  Description of blank BIOMT
    '''

    biomt_rec = {
        'subdivs': subdivs,
        'auth_bio_unit': '',
        'soft_bio_unit': '',
        'rot': [],
        'trans': []
    }

    biomt_rec['rot'].append(np.identity(3))
    biomt_rec['trans'].append(np.array([0.0, 0.0, 0.0]))

    return biomt_rec


def add_transform(biomt_rec):
    '''
    Get BIOMT transformations from user and add to record

    @type biomt_rec : dict
    @param biomt_rec: Description of unit transforms to create biological unit
    @return:
    '''

    input_txt = ''
    valid_transform = False

    while not valid_transform and input_txt not in ['x', 'X']:

        valid_transform, rot, trans = get_user_transform()

        if valid_transform:
            biomt_rec['rot'].append(rot)
            biomt_rec['trans'].append(trans)
        else:
            print('No transform added, to return to other '
                  'options press "X", any other key to try again.')
            input_txt = sys.stdin.read(1)

    return


def edit_transform(biomt_rec):
    '''
    Allow user editing of input BIOMT data

    @type biomt_rec : dict
    @param biomt_rec: Description of unit transforms to create biological unit
    @return:
    '''

    print_biomt(biomt_rec)

    rec_no = 0
    no_entries = len(biomt_rec.rot)

    if no_entries > 1:

        while rec_no == 0:

            print("\nSelect record to edit")
            txt = input()
            try:
                rec_no = int(txt)
                if rec_no <= 0 or rec_no >= no_entries:
                    print('Invalid record selected')
                    rec_no = 0
            except:
                print('Invalid record selected')
                rec_no = 0

        input_txt = ''
        valid_transform = False

        while not valid_transform and input_txt not in ['x', 'X']:

            valid_transform, rot, trans = get_user_transform()

            if valid_transform:
                biomt_rec['rot'][rec_no] = rot
                biomt_rec['trans'][rec_no] = trans
            else:
                print('Transform not edited, to return to other options '
                      'press "X", any other key to try again.')
                input_txt = sys.stdin.read(1)

    else:

        print('Only the identity tranformation exists, add a new transform')

    return


def print_biomt(biomt_rec):
    '''
    Print the existing data in a BIOMT record

    @type biomt_rec : dict
    @param biomt_rec: Description of unit transforms to create biological unit
    @return:
    '''

    no_transforms = len(biomt_rec.rot)

    print('0. Identity (cannot be edited)')

    for i in range(1, no_transforms):
        print('#' + str(i))
        print('Rotation:')
        print(biomt_rec.rot[i])
        print('Translation:')
        print(biomt_rec.trans[i])

    return


def select_segnames(segnames):
    '''
    Prompt user for a selection fo segment names to apply BIOMT to.

    @type segnames :  list
    @param segnames:  Valid segnames contained in structure
    @return:
    '''

    seg_list = ', '.join(segnames)

    print("Select list of segments for the transformation to be applied to "
          "using Python list notation e.g. ['A','B'])")
    print("Available segments are: " + seg_list)

    selected = ''

    while not selected:

        try:
            selected = eval(input())
            if not set(selected).issubset(set(segnames)):
                print('All values in list must be available segments')
                selected = ''
        except:
            print('Invalid input, try again')
            selected = ''

    return list(selected)


def prepare_biomt_json(biomt_rec):
    '''
    Convert BIOMT information into JSOn for output

    @type biomt_rec : dict
    @param biomt_rec: Description of unit transforms to create biological unit
    @return:
    '''

    json_rec = deepcopy(biomt_rec)

    rot = json_rec['rot']
    trans = json_rec['trans']

    for i in range(len(rot)):
        rot[i] = rot[i].tolist()
        trans[i] = trans[i].tolist()

    return json_rec


def user_biomt(segnames_json):
    '''
    Get user to edit BIOMT information from the command line

    @type  segnames_json:  file
    @param segnames_json:  Valid segment names in JSON format
    @return:
    '''

    segnames = yaml.safe_load(segnames_json)

    selected_segnames = select_segnames(segnames)

    biomt_rec = init_biomt(selected_segnames)

    print('An identity transform already exists, enter a new transform:')
    add_transform(biomt_rec)

    input_txt = ''

    while input_txt not in ['q', 'Q']:

        #print('Press "E" to edit an existing transform, '
        print('Press "A" to add a new transform or '
              '"Q" to accept current BIOMT and quit')
        input_txt = sys.stdin.read(1)

        #if input_txt in ['e', 'E']:
        #    edit_transform(biomt_rec)
        if input_txt in ['a', 'A']:
            add_transform(biomt_rec)

    return json.dumps(prepare_biomt_json(biomt_rec))
