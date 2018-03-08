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

