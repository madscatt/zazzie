#! /usr/bin/env python

import sys
import numpy as np
import json
import yaml
from copy import deepcopy

def check_array_input(input_txt, dimensions):

    flag = True

    try:
        matrix = np.array(input_txt)
        matrix = 1.0 * matrix
        if matrix.shape != dimensions:
            flag = False

    except TypeError:
        matrix = None

    return flag, matrix

def check_rotation(input_txt):

    flag, matrix = check_array_input(input_txt,(3,3))

    return flag, matrix

def check_translation(input_txt):

    flag, matrix = check_array_input(input_txt,(3,))

    return flag, matrix

def get_user_transform():

    flag = True
    rot = None
    trans = None

    print "Enter rotation matrix (3 x 3, using Python list notation i.e.: [[1,0,0],[0,1,0],[0,0,1]]):"
    rot_txt = input()
    check_rot, rot = check_rotation(rot_txt)

    if not check_rot:
        flag = False
        print "Rotation must be a 3 x 3 array of numeric values"

    else:

        print "Enter translation vector (3 x 1, using Python list notation i.e.: [1,0,0]):"
        trans_txt = input()
        check_trans, trans = check_translation(trans_txt)

        if not check_trans:
            flag = False
            print "Error: Translation must be a 3 x 1 array of numeric values"

    if ((rot == np.identity(3)).all()) and ((trans == np.array([0.0,0.0,0.0])).all()):
        print "A second identity transform will not be added"
        flag = False

    return flag, rot, trans

def init_biomt(subdivs):

    biomt_rec =  {
                    'subdivs': subdivs,
                    'auth_bio_unit': '',
                    'soft_bio_unit': '',
                    'rot': [],
                    'trans': []
                 }

    biomt_rec['rot'].append(np.identity(3))
    biomt_rec['trans'].append(np.array([0.0,0.0,0.0]))

    return biomt_rec

def add_transform(biomt_rec):

    input_txt = ''
    valid_transform = False

    while not valid_transform and not input_txt in ['x','X']:

        valid_transform, rot, trans = get_user_transform()

        if valid_transform:
            biomt_rec['rot'].append(rot)
            biomt_rec['trans'].append(trans)
        else:
            print 'No transform added, to return to other options press "X", any other key to try again.'
            input_txt = sys.stdin.read(1)

    return

def edit_transform(biomt_rec):

    print_biomt(biomt_rec)

    rec_no = 0
    no_entries = len(biomt_rec.rot)

    if no_entries > 1:

        while rec_no == 0:

            print "\nSelect record to edit"
            txt = raw_input()
            try:
                rec_no = int(txt)
                if rec_no <= 0 or rec_no >= no_entries:
                    print 'Invalid record selected'
                    rec_no = 0
            except:
                print 'Invalid record selected'
                rec_no = 0

        input_txt = ''
        valid_transform = False

        while not valid_transform and not input_txt in ['x','X']:

            valid_transform, rot, trans = get_user_transform()

            if valid_transform:
                biomt_rec['rot'][rec_no] = rot
                biomt_rec['trans'][rec_no] = trans
            else:
                print 'Transform not edited, to return to other options press "X", any other key to try again.'
                input_txt = sys.stdin.read(1)

    else:

        print 'Only the identity tranformation exists, add a new transform'

    return

def print_biomt(biomt_rec):

    no_transforms = len(biomt_rec.rot)

    print '0. Identity (cannot be edited)'

    for i in range(1,no_transforms):

        print '#' + str(i)
        print 'Rotation:'
        print biomt_rec.rot[i]
        print 'Translation:'
        print biomt_rec.trans[i]

    return

def select_segnames(segnames):

    seg_list = ', '.join(segnames)

    print "Select list of segments for the transformation to be applied to (using Python list notation e.g. ['A','B'])"
    print "Available segments are: " + seg_list

    selected = ''

    while not selected:

        try:
            selected = input()
            if not set(selected).issubset(set(segnames)):
                print 'All values in list must be available segments'
                selected = ''
        except:
            print 'Invalid input, try again'
            selected = ''

    return list(selected)

def prepare_biomt_json(biomt_rec):

    json_rec = deepcopy(biomt_rec)

    rot = json_rec['rot']
    trans = json_rec['trans']

    for i in range(len(rot)):
        rot[i] = rot[i].tolist()
        trans[i] = trans[i].tolist()

    return json_rec

def user_biomt(segnames_json):

    segnames = yaml.safe_load(segnames_json)

    selected_segnames = select_segnames(segnames)

    biomt_rec = init_biomt(selected_segnames)

    print 'An identity transform already exists, enter a new transform:'
    add_transform(biomt_rec)

    input_txt = ''

    while input_txt not in ['q','Q']:

        print 'Press "E" to edit an existing transform, "A" to add a new transform or "Q" to accept current BIOMT and quit'
        input_txt = sys.stdin.read(1)

        if input_txt in ['e','E']:
            edit_transform(biomt_rec)
        elif input_txt in ['a','A']:
            add_transform(biomt_rec)

    return json.dumps(prepare_biomt_json(biomt_rec))