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

import string

'''
    PDB_RECORD_READER

    This provides schemas and line parsing functions for reading PDBs

'''


def split(input_str):
    """
    Wraps string.split so it can be used as a function in the parser

    @type  input_str:  string
    @param input_str:  String to be split
    @rtype:            List
    @return:           List of strings split at whitespace from input
    """

    return string.split(input_str)

# Schemas for interprating different PDB records and creating dictionary
# Schema format:
#
#    1) record name
#    2) start index
#    3) stop as for Python range (i.e. no. characters counting from 0) or None for single character
#    4) processing function or None for `string.strip`
#
#    so a truncated (but valid) schema for a PDB ATOM record cold be:
#
#    ATOM_RECORD_SCHEMA = (
#    ...    ('record', 0,  6,   None),   # 0,  record name
#    ...    ('name', 6, 11,   int),      # 1,  serial
#    )

rec_schemas = {

    # Commented out schemas are for coordinate record types which are handled by
    # the standard sasmol PDB reader

    #    'ATOM': (
    #        ('record', 0, 6, None),       # 0,  record type name
    #        ('index', 6, 11, int),        # 1,  atom serial/index no
    #        ('atomname', 12, 16, None),   # 2,  atom name
    #        ('alt', 16, None, None),      # 3,  altLoc
    #        ('resname', 17, 20, None),    # 4,  residue name
    #        ('chain', 21, None, None),    # 5,  chain ID
    #        ('resid', 22, 26, int),       # 6,  residue sequence number
    #        ('insert', 26, None, None),   # 7,  insertion code
    #        ('x', 30, 38, float),         # 8,  x
    #        ('y', 38, 46, float),         # 9,  y
    #        ('z', 46, 54, float),         # 10, z
    #        ('occ', 54, 60, float),       # 11, occupancy
    #        ('b', 60, 66, float),         # 12, temperature factor
    #        ('segid', 72, 76, None),      # 13, segment ID
    #        ('element', 76, 78, None),    # 14, element
    #        ('charge', 78, 80, None),     # 15, charge
    #    ),
    #
    #    'HETATM': (
    #        ('record', 0, 6, None),       # 0,  record type name
    #        ('index', 6, 11, int),        # 1,  atom serial/index no
    #        ('atomname', 12, 16, None),   # 2,  atom name
    #        ('alt', 16, None, None),      # 3,  altLoc
    #        ('resname', 17, 20, None),    # 4,  residue name
    #        ('chain', 21, None, None),    # 5,  chain ID
    #        ('resid', 22, 26, int),       # 6,  residue sequence number
    #        ('insert', 26, None, None),   # 7,  insertion code
    #        ('x', 30, 38, float),         # 8,  x
    #        ('y', 38, 46, float),         # 9,  y
    #        ('z', 46, 54, float),         # 10, z
    #        ('occ', 54, 60, float),       # 11, occupancy
    #        ('b', 60, 66, float),         # 12, temperature factor
    #        ('segid', 72, 76, None),      # 13, segment ID
    #        ('element', 76, 78, None),    # 14, element
    #        ('charge', 78, 80, None),     # 15, charge
    #    ),
    #
    #    'TER': (
    #        ('record', 0, 6, None),  # 0,  record type name
    #        ('text', 6, 80, str),    # 1,  text
    #    ),
    #
    #    'END': (
    #        ('record', 0, 6, None),  # 0,  record type name
    #        ('text', 6, 80, str),    # 1,  text
    #    ),
    #
    #    'ENDMDL': (
    #        ('record', 0, 6, None),  # 0,  record type name
    #        ('text', 6, 80, str),    # 1,  text
    #    ),

    'MASTER': (
        ('record', 0, 6, None),  # 0,  record type name
        ('text', 6, 80, str),    # 1,  text
    ),

    #    'ANISOU': (
    #        ('record', 0, 6, None),       # 0,  record type name
    #        ('index', 6, 11, int),        # 1,  atom serial/index no
    #        ('atomname', 12, 16, None),   # 2,  atom name
    #        ('alt', 16, None, None),      # 3,  altLoc
    #        ('resname', 17, 20, None),    # 4,  residue name
    #        ('chain', 21, None, None),    # 5,  chain ID
    #        ('resid', 22, 26, int),       # 6,  residue sequence number
    #        ('insert', 26, None, None),   # 7,  insertion code
    #        ('u11', 28, 35, int),         # 8,  U(1,1)
    #        ('u22', 35, 42, int),         # 9,  U(2,2)
    #        ('u33', 42, 49, int),         # 10, U(3,3)
    #        ('u12', 50, 56, int),         # 8,  U(1,2)
    #        ('u13', 56, 63, int),         # 9,  U(1,3)
    #        ('u23', 64, 70, int),         # 10, U(2,3)
    #        ('element', 76, 78, None),   # 11, element
    #        ('charge', 78, 80, None),     # 12, charge
    #
    #    ),

    'SITE': (
        ('record', 0, 6, None),  # 0,  record type name
        ('text', 6, 80, str),    # 1,  text
    ),

    'ORIGXn': (
        ('record', 0, 6, None),  # 0,  record type name
        ('on1', 10, 20, float),  # 1.  On1
        ('on2', 20, 30, float),  # 2.  On2
        ('on3', 30, 40, float),  # 3.  On3
        #('tn', 45, 55, float),   # 4,  Tn - sometimes missing
        ('tn', 45, 55, None),    # 4,  Tn
    ),

    'SCALEn': (
        ('record', 0, 6, None),  # 0,  record type name
        ('sn1', 10, 20, float),  # 1.  Sn1
        ('sn2', 20, 30, float),  # 2.  Sn2
        ('sn3', 30, 40, float),  # 3.  Sn3
        ('un', 45, 55, float),   # 4,  Un
    ),

    'HET': (
        ('record', 0, 6, None),       # 0,  record type name
        ('resname', 7, 10, None),     # 1,  residue name
        ('chain', 12, None, None),    # 2,  chain ID
        ('resid', 13, 17, int),       # 3,  residue sequence number
        ('insert', 17, None, None),   # 4,  insertion code
        ('natoms', 21, 25, int),      # 5,  number of atoms
        ('text', 30, 70, None),       # 6,  description
    ),

    'HETNAM': (
        ('record', 0, 6, None),    # 0,  record type name
        ('cont', 8, None, None),   # 1,  continuation line no.
        ('resname', 11, 14, None),  # 2,  hetatm id
        ('text', 15, 70, None),    # 3,  chemical name
    ),

    'CRYST1': (
        # unit cell parameters, space group, and Z value
        ('record', 0, 6, None),         # 0,  record type name
        ('a', 6, 15, float),            # 1,  a (Angstrom)
        ('b', 15, 24, float),           # 2,  b (Angstrom)
        ('c', 24, 33, float),           # 3,  c (Angstrom)
        ('alpha', 33, 40, float),       # 4,  alpha
        ('beta', 40, 47, float),        # 5,  beta
        ('gamma', 47, 54, float),       # 6,  gamma
        ('space_group', 55, 66, None),  # 7.  space group
        ('z', 66, 70, int),             # 8.  Z value
    ),

    'REMARK': (
        ('record', 0, 6, None),  # 0,  record type name
        ('num', 7, 11, int),     # 1,  remark number (type)
        ('text', 11, 79, str),   # 2,  text
    ),

    'HEADER': (
        ('record', 0, 6, None),    # 0,  record type name
        ('class', 10, 50, None),   # 1,  classification
        ('date', 50, 59, None),    # 2,  date
        ('pdb_id', 62, 66, None),  # 3,  PDB ID of this file
    ),

    'TITLE': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 10, 80, None),  # 2,  text
    ),

    'JRNL': (
        ('record', 0, 6, None),  # 0,  record type name
        ('text', 12, 80, None),  # 1,  text
    ),

    'OBSLTE': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 11, 75, None),  # 2,  text
    ),

    'SPLIT': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 11, 80, None),  # 2,  text
    ),

    'CAVEAT': (
        ('record', 0, 6, None),    # 0,  record type name
        ('cont', 8, 10, None),     # 1,  continuation line no.
        ('pdb_id', 11, 15, None),  # 2,  PDB ID of this file
        ('text', 15, 79, None),    # 3,  text
    ),

    'COMPND': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 10, 80, None),  # 2,  text
    ),

    'SOURCE': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 10, 79, None),  # 2,  text
    ),

    'KEYWDS': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 10, 79, None),  # 2,  text
    ),

    'EXPDTA': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 10, 79, None),  # 2,  text
    ),

    'NUMMDL': (
        ('record', 0, 6, None),       # 0,  record type name
        ('no_models', 10, 14, None),  # 1,  number of models
    ),

    'MDLTYP': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 10, 80, None),  # 2,  text
    ),

    'AUTHOR': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 10, 79, None),  # 2,  text
    ),

    'REVDAT': (
        ('record', 0, 6, None),    # 0,  record type name
        ('mod_num', 7, 10, None),  # 1,  modification no.
        ('cont', 10, 12, None),    # 2,  continuation line no.
        ('text', 13, 66, None),    # 3,  text
    ),

    'SPRSDE': (
        ('record', 0, 6, None),  # 0,  record type name
        ('cont', 8, 10, None),   # 1,  continuation line no.
        ('text', 11, 75, None),  # 2,  text
    ),

    'DBREF': (
        ('record', 0, 6, None),  # 0,  record type name
        ('text', 7, 80, None),   # 1,  text
    ),

    'DBREF1': (
        ('record', 0, 6, None),  # 0,  record type name
        ('text', 7, 80, None),   # 1,  text
    ),

    'DBREF2': (
        ('record', 0, 6, None),  # 0,  record type name
        ('text', 7, 80, None),   # 1,  text
    ),

    'SEQADV': (
        ('record', 0, 6, None),        # 0,  record type name
        ('pdb_id', 7, 11, None),       # 1,  PDB ID
        ('resname', 12, 15, None),     # 2,  residue name
        ('chain', 16, None, None),     # 3,  chain ID
        #('resid', 18, 22, int),        # 4,  residue sequence number - missing in some records
        ('resid', 18, 22, None),       # 4,  residue sequence number
        ('insert', 23, None, None),    # 5,  insertion code
        ('db', 24, 28, None),          # 6,  database
        ('dbid', 29, 38, None),        # 7,  database accession number
        ('db_resname', 39, 42, None),  # 8,  database residue name
        ('db_resid', 43, 48, None),    # 9,  database sequence number
        ('text', 49, 70, None),        # 10, conflict comment

    ),

    'SEQRES': (
        ('record', 0, 6, None),        # 0,  record type name
        ('cont', 7, 10, None),         # 1,  serial number
        ('chain', 11, None, None),     # 2,  chain ID
        ('num_res', 13, 17, int),      # 3,  number of residues in chain
        ('resnames', 19, 70, split),   # 4,  list of residue names
    ),


    'MODRES': (
        ('record', 0, 6, None),      # 0,  record type name
        ('pdb_id', 7, 11, None),     # 1,  PDB ID
        ('resname', 12, 15, None),   # 2,  residue name
        ('chain', 16, None, None),   # 3,  chain ID
        ('resid', 18, 22, int),      # 4,  residue sequence number
        ('insert', 22, None, None),  # 5,  insertion code
        ('std_res', 24, 27, None),   # 6,  name of standard residue
        ('text', 29, 70, None),      # 7,  description
    ),

    'HETSYN': (
        ('record', 0, 6, None),     # 0,  record type name
        ('cont', 8, 10, None),      # 1,  continuation line no.
        ('resname', 11, 14, None),  # 2,  residue name for HETATM
        ('text', 15, 70, None),     # 3,  text
    ),

    'FORMUL': (
        ('record', 0, 6, None),     # 0,  record type name
        ('comp_no', 8, 10, None),   # 1,  component no.
        ('resname', 12, 15, None),  # 2,  residue name for HETATM
        ('cont', 16, 18, None),     # 3,  continuation no.
        ('water', 18, None, None),  # 4,  '*' for water
        ('text', 19, 70, None),     # 5,  text
    ),


    'SSBOND': (
        ('record', 0, 6, None),       # 0,  record type name
        ('serial', 7, 10, int),       # 1,  serial number
        ('resname1', 11, 14, None),   # 2,  residue name 1
        ('chain1', 15, None, None),   # 3,  chain ID 1
        ('resid1', 17, 21, int),      # 4.  residue sequence number 1
        ('insert1', 21, None, None),  # 5,  insertion code1
        ('resname2', 25, 28, None),   # 6,  residue name 2
        ('chain2', 29, None, None),   # 7,  chain ID 2
        ('resid2', 31, 35, int),      # 8.  residue sequence number 2
        ('insert2', 35, None, None),  # 9,  insertion code 2
        ('sym1', 59, 65, int),        # 10,  symmetry operator for 1st residue
        ('sym2', 66, 72, int),        # 11,  symmetry operator for 2nd residue
    ),

    'HELIX': (
        ('record', 0, 6, None),       # 0,  record type name
        ('serial', 7, 10, int),       # 1,  serial number
        ('id', 11, 14, None),         # 2,  helix ID
        ('resname1', 15, 18, None),   # 3,  residue name start of helix (1)
        ('chain1', 19, None, None),   # 4,  chain ID 1
        ('resid1', 21, 25, int),      # 5.  residue sequence number 1
        ('insert1', 25, None, None),  # 6,  insertion code1
        ('resname2', 27, 30, None),   # 7,  residue name end of helix (2)
        ('chain2', 31, None, None),   # 8,  chain ID 2
        ('resid2', 33, 37, int),      # 9.  residue sequence number 2
        ('insert2', 37, None, None),  # 10,  insertion code 2
        ('class', 38, 40, int),       # 11,  helix class
        ('text', 40, 70, None),       # 12,  comment
        ('length', 71, 76, int),      # 13, helix length
    ),

    'SHEET': (
        ('record', 0, 6, None),          # 0,  record type name
        ('strand', 7, 10, int),          # 1,  strand number
        ('id', 11, 14, None),            # 2,  sheet ID
        ('no_strands', 14, 16, None),    # 3,  number of strands
                                         # no conversion as can be blank
        ('resname1', 17, 20, None),      # 4,  residue name start of strand (1)
        ('chain1', 21, None, None),      # 5,  chain ID 1
        ('resid1', 22, 26, int),         # 6.  residue sequence number 1
        ('insert1', 26, None, None),     # 7,  insertion code1
        ('resname2', 28, 31, None),      # 8,  residue name end of strand (2)
        ('chain2', 32, None, None),      # 9,  chain ID 2
        ('resid2', 33, 37, int),         # 10.  residue sequence number 2
        ('insert2', 37, None, None),     # 11,  insertion code 2
        #('sense', 38, 40, int),
        ('sense', 38, 40, None),         # 12,  sense of strand (wrt last strand)
                                         # No conversion as can be blank
        ('bonding', 41, 70, None)        # 13,  bonding information
        #
        #        ('cur_atomname', 41, 45, None),  # 13,  atom in current strand
        #        ('cur_resname', 45, 48, None),   # 14,  residue name in current strand
        #        ('cur_chain', 49, None, None),   # 15,  chain in current strand
        #        ('cur_resid', 50, 54, int),      # 16,  sequence number in current strand
        #        ('cur_insert', 54, None, None),  # 17,  insertion code in current strand
        #        ('pre_atomname', 56, 60, None),  # 18,  atom in previous strand
        #        ('pre_resname', 60, 63, None),   # 19,  residue name in previous strand
        #        ('pre_chain', 64, None, None),   # 20,  chain in previous strand
        #        ('pre_resid', 65, 69, int),      # 21,  sequence number in previous strand
        #        ('pre_insert', 69, None, None),    # 22,  insertion code in previous strand
    ),

    'TURN': (
        ('record', 0, 6, None),       # 0,  record type name
        ('serial', 7, 10, int),       # 1,  serial number
        ('id', 11, 14, None),         # 2,  helix ID
        ('resname1', 15, 18, None),   # 3,  residue name start of helix (1)
        ('chain1', 19, None, None),   # 4,  chain ID 1
        ('resid1', 20, 24, int),      # 5.  residue sequence number 1
        ('insert1', 24, None, None),  # 6,  insertion code1
        ('resname2', 26, 29, None),   # 7,  residue name end of helix (2)
        ('chain2', 30, None, None),   # 8,  chain ID 2
        ('resid2', 31, 35, int),      # 9.  residue sequence number 2
        ('insert2', 35, None, None),  # 10,  insertion code 2
        ('text', 40, 70, None),       # 11,  comment
    ),

    'LINK': (
        ('record', 0, 6, None),       # 0,  record type name
        ('atomname1', 12, 16, None),  # 1,  atom name 1
        ('altloc1', 16, None, None),  # 2,  alternate location indicator 1
        ('resname1', 17, 20, None),   # 3,  residue name 1
        ('chain1', 21, None, None),   # 4,  chain ID 1
        ('resid1', 22, 26, int),      # 5.  residue sequence number 1
        ('insert1', 26, None, None),  # 6,  insertion code1
        ('atomname2', 42, 46, None),  # 7,  atom name 2
        ('altloc2', 46, None, None),  # 8,  alternate location indicator 2
        ('resname2', 47, 50, None),   # 9,  residue name 2
        ('chain2', 51, None, None),   # 10,  chain ID 2
        ('resid2', 52, 56, int),      # 11.  residue sequence number 2
        ('insert2', 56, None, None),  # 12,  insertion code 2
        ('sym1', 59, 65, int),        # 13,  symmetry operator for 1st residue
        ('sym2', 66, 72, int),        # 14,  symmetry operator for 2nd residue
    ),

    'CISPEP': (
        ('record', 0, 6, None),       # 0,  record type name
        ('serial', 7, 10, int),       # 1,  serial number
        ('resname1', 11, 14, None),   # 2,  residue name 1
        ('chain1', 15, None, None),   # 3,  chain ID 1
        ('resid1', 17, 21, int),      # 4.  residue sequence number 1
        ('insert1', 21, None, None),  # 5,  insertion code1
        ('resname2', 25, 28, None),   # 6,  residue name end of helix (2)
        ('chain2', 29, None, None),   # 7,  chain ID 2
        ('resid2', 31, 35, int),      # 8.  residue sequence number 2
        ('insert2', 35, None, None),  # 9,  insertion code 2
        ('model', 43, 46, int),       # 10,  specific model
        ('angle', 53, 59, None),      # 11,  angle measurement
    ),

    #    'CONECT': (
    #        ('record', 0, 6, None),     # 0,  record type name
    #        ('resnames', 6, 31, None),  # 1,  list of atom indices
    #    ),

    'NONSTD': (
        ('record', 0, 6, None),  # 0,  record type name
        ('text', 6, 80, None),   # 1,  text
    ),

}


def parse_line(line, schema):
    """
    Parse an record line from a PDB via schema.
    Based on the PDB ATOM parsing code written by Alisue (lambdalisue@hashnote.net):
    https://gist.github.com/lambdalisue/7209305

    Schema format::

        1) record name
        2) start index
        3) stop as for Python range (i.e. no. characters counting from 0) or
           None for single character
        4) processing function or None for `string.strip`

        so a truncated (but valid) schema for a PDB ATOM record could be:

        ATOM_RECORD_SCHEMA = (
        ...    ('record', 0,  6,   None),   # 0,  record name
        ...    ('name', 6, 11,   int),      # 1,  serial
        )

    @type  line:   string
    @param line:   Line from a text file.
    @type  schema: list
    @param schema: List of 4-tuples. The format of the tuples is given above.
    @rtype:        dictionary, list
    @return:
                   0. Dictionary with keys specific to the record name in the
                      schema.
                   1. List of errors encountered (if any)
    """

    vals = {}
    err = []

    for record, start, end, fn in schema:
        if end is None:
            end = start + 1
        if fn is None:
            fn = string.strip
        try:
            vals[record] = fn(line[start:end])
        except ValueError as e:
            err.append("Python error message: {0:s}".format(e))

    return vals, err


def process_runon_line(rec_lines, element):
    """
    Combine strings where information in a PDB record is split between multiple
    lines.

    @type  rec_lines:  list
    @param rec_lines:  List of strings to be combine to form single line
    @type  element:    string
    @param element:    Key to record to be combined
                       (record name as key)
    @rtype:            string
    @return:           All element entries of rec_lines joing by a space

    """

    return ' '.join(list(x[element] for x in rec_lines))
