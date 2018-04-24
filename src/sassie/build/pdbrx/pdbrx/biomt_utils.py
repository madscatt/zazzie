# -*- coding: utf-8 -*-
"""
Methods to process information for biomt data

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

import yaml
import numpy


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

    biomt_rec['rot'].append(numpy.identity(3))
    biomt_rec['trans'].append(numpy.array([0.0, 0.0, 0.0]))

    return biomt_rec

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
        matrix = numpy.array(input_txt)
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



def check_biological_unit(biomt_unit):
    """
    Check biological unit transform from user is valid

    @return:
    """

    valid = True

    recs = set(['subdivs', 'auth_bio_unit',
                    'soft_bio_unit', 'rot', 'trans'])

    if recs == set(biomt_unit.keys()):
        pass

    else:

        valid = False

    return valid

def biomt_json2data(json_biomt_rec):
    """
    Convert JSON format BIOMT records into numpy arrays for use in
    coordinate transforms

    @type json_biomt_rec :  str
    @param json_biomt_rec:  BIOMT records in JSON format
    @rtype :
    @return:
    """

    biomt_rec = yaml.safe_load(json_biomt_rec)
    valid = check_biological_unit(biomt_rec)

    if valid:

        for i in range(len(biomt_rec['rot'])):
            biomt_rec['rot'][i] = numpy.array(biomt_rec['rot'][i])
            biomt_rec['trans'][i] = numpy.array(biomt_rec['trans'][i])

    else:

        biomt_rec = None

    return biomt_rec

