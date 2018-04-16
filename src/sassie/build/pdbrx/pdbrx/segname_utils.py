# -*- coding: utf-8 -*-
"""
Methods to process segment based information

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

import json
import numpy
import sassie.build.pdbscan.pdbscan.data_struct as data_struct

def split_segnames(other_self, mol, ndx, new_segname):
    """
    Split an existing segment and name the newly created segment.

    @type ndx :  int
    @param ndx:  Index of the residue selected for segment break (in list
                     of residue descritions)
    @type new_segname :  str
    @param new_segname:  Name to be applied to the newly created segment
    @return:
    """

    resid_desc = other_self.resid_descriptions

    last_ndx = len(other_self.resid_descriptions) - 1

    current_segname = other_self.resid_descriptions[ndx][0]


    if ndx != 0:
        previous_segname = other_self.resid_descriptions[ndx - 1][0]
    else:
        previous_segname = ''

    if previous_segname == current_segname:

        updated_data = []

        for i in range(len(other_self.resid_descriptions)):

            line = other_self.resid_descriptions[i]

            if i >= ndx and other_self.resid_descriptions[i][0] == current_segname:
                line[0] = new_segname

            updated_data.append(line)

        other_self.resid_descriptions = numpy.array(updated_data)

        mol.segnames().append(new_segname)

    return

def rename_segment(other_self, mol, ndx, new_segname):
    """
    Change the name of selected segment (the one including the selected
        residue).

    @type ndx :  int
    @param ndx:  Index of the user selected residue
    @type new_segname :  str
    @param new_segname:  New name for segment
    @return:
    """

    target_segname = other_self.resid_descriptions[ndx][0]

    updated_data = []

    for line in other_self.resid_descriptions:
        if line[0] == target_segname:
            line[0] = new_segname
        updated_data.append(line)

    other_self.resid_descriptions = numpy.array(updated_data)

    #dumfile = open('dum.txt', 'w')
    #dumfile.write('#ORIGINAL : \n')
    #dumfile.write('type(mol.segnames()) \n' + str(type(mol.segnames()))+ '\n')
    #for value in mol.segnames():
    #    dumfile.write(value + "\n")

    temp_segnames = [x if (x != target_segname)
                         else new_segname for x in mol.segnames()]

    mol.setSegnames(temp_segnames)

    #dumfile.write('#FINAL : \n')
    #dumfile.write('type(mol.segnames()) \n' + str(type(mol.segnames())) + '\n')
    #dumfile.write('#mol.segnames() : \n')
    #for value in mol.segnames():
    #    dumfile.write(value + "\n")

    #dumfile.close()

    ## assign new segname for each residue if changed (DEBUGGING)
    ## assign new segname for each residue if changed (DEBUGGING)
    ## assign new segname for each residue if changed (DEBUGGING)

    temp_segname = []

    for this_segname in mol.segname():

        if this_segname not in mol.segnames():
            temp_segname.append(new_segname)      
        else:
            temp_segname.append(this_segname)      

    mol.setSegname(temp_segname)

    #dumfile.write('#mol.segname() : \n')
    #for value in mol.segname():
    #    dumfile.write(value + "\n")
#
#/    dumfile.close()

    return

def join_segnames(other_self, mol, ndx):
    """
    Join segment containing ndx-th residue to the previous segment.

    @type ndx:          integer
    @param ndx:         Index of the residue that starts segment to join
                            previous segment
    @return:            Updated array of lines containing: segnames, indices,
                            resids, resnames, chains, moltypes
    """

    last_ndx = len(other_self.resid_descriptions) - 1

    current_segname = other_self.resid_descriptions[ndx][0]
    current_moltype = other_self.resid_descriptions[ndx][-1]

    if ndx != 0:

        previous_segname = other_self.resid_descriptions[ndx - 1][0]
        previous_moltype = other_self.resid_descriptions[ndx - 1][-1]

        moltype_match = (previous_moltype == current_moltype)
        resid_match = (other_self.resid_descriptions[ndx - 1][2] < other_self.resid_descriptions[ndx][2])

    else:
        previous_segname = ''
        # No previous segment, so joining makes no sense
        moltype_match = False
        resid_match = False

    segname_mismatch = (previous_segname != current_segname)

    acceptable_join = moltype_match and resid_match and segname_mismatch

    error = ''

    if acceptable_join:

        updated_data = []

        for i in range(len(other_self.resid_descriptions)):

            line = other_self.resid_descriptions[i]

            if i >= ndx and other_self.resid_descriptions[i][0] == current_segname:
                line[0] = previous_segname

            updated_data.append(line)

        other_self.resid_descriptions = numpy.array(updated_data)

        mol.segnames.remove(current_segname)

    else:

        if not segname_mismatch:
            error = 'Segments with the same name cannot be joined'
        elif not resid_match:
            error = 'Joined segment must start with higher resid'
        else:
            error = 'Joined segments must have same moltype'

    return error

def get_segment_starts(other_self):
    """
    Get indicies where the resid descriptions change segment name.

    @return:
    """

    new_breaks = numpy.where(other_self.resid_descriptions[:-1, 0] != other_self.resid_descriptions[1:, 0])[0]

    if (new_breaks != other_self.starting_breaks).any():

        new_breaks += 1
        new_breaks = numpy.append([0], new_breaks)

        start_segnames = {}

        # Note residue descriptions give last atom in that residue
        # account for that here
        for start_ndx in new_breaks:
            if start_ndx == 0:
                sasmol_index = 0
            else:
                sasmol_index = int(other_self.resid_descriptions[start_ndx - 1][1]) + 1

            start_segnames[sasmol_index] = other_self.resid_descriptions[start_ndx][0]

    else:
        start_segnames = {}

    return json.dumps(start_segnames)

def valid_segname(segname, segnames):
    """
    Check that the input segment name is valid to use for a new segment,
    i.e is 4 characters long or less and not an existing segment name.

    @type segname :  str
    @param segname:  Proposed segment name
    @rtype :  bool
    @return:  Is the input segname valid
    """

    valid = False

    if len(segname) <= 4 and segname not in segnames:
        valid = True

    return valid

def convert_segname_start(data):
    """
    Convert JSON segname_start dictionary to have integer keys and ascii
    strings (as opposed to unicode)

    @rtype :  dictionary
    @return:  Keys = atom index of start of segment as integer,
    Value = segname
    """

    return dict((int(x[0]), x[1].encode('ascii')) for x in data.items())


def redefine_segments(mol, segname_starts):
    """
    Alter segmentation in mol to follow user input then scan to obtain
    information for further processing.

    @type  segname_starts:  dictionary
    @param segname_starts:  Keys = atom index of start of segment,
                                Value = segname
    """

    update_segments(mol, segname_starts)

    # Blank the previous segname_info
    mol.segname_info = data_struct.Info(scan_type='segname')

    # Get updated sequence/missing residue etc. information
    mol.segment_scan(initialize=False)

    return

def update_segments(mol, segname_starts):
    """
    Alter segmentation in mol to follow user input.

    @type  segname_starts:  dictionary
    @param segname_starts:  Keys = atom index of start of segment,
                                Value = segname
    """

    natoms = mol.natoms()

    segnames = []
    indices = mol.index()

    new_segname = segname_starts[0]

    for ndx in range(natoms):

        atom_index = indices[ndx]

        if indices[ndx] in segname_starts:
            new_segname = segname_starts[atom_index]

        segnames.append(new_segname)

    mol.update_segnames(segnames)

    return

def get_segment_moltype(mol, segname):
    """
    Get the Moltype of the chosen segment

    @type segname :  str
    @param segname:  Name of segment for which moltype is required
    @rtype : str
    @return: Moltype of segmanent
    """

    mapping = set(zip(mol.segname(), mol.moltype()))
    moltypes = [x[1] for x in mapping if x[0] == segname]

    moltype = '/'.join(moltypes)

    return moltype


def create_residue_descriptions(mol):
    """
    Filter information from self.mol to provide residue descriptions for
    used in displaying structure contents for segment editing.

    @rtype : list
    @return: List of tuples describing segname, first atomic index,
        resid, resname, chain and moltype for each residue
    """

    segnames = mol.segname()
    indices = mol.index()
    resids = mol.resid()
    resnames = mol.resname()
    chains = mol.chain()
    moltypes = mol.moltype()

    data = numpy.array(zip(segnames, indices, resids,
                            resnames, chains, moltypes))

    mask = numpy.where(data[:-1, 2] != data[1:, 2])[0]
    mask = numpy.append(mask, [len(data) - 1])

    residue_descriptions = data[mask]

    return residue_descriptions

def create_sequence_report(mol, seq_segnames): 

    sequence_report = 'Current residue sequences for each segment (uppercase letters have coordinates while lowercase letters do not have coordinates): \n\n'

    for segname in seq_segnames:
        seq = mol.segname_info.sequence_to_fasta(
                    segname, missing_lower=True)
        sequence_report += 'segname ' + segname + ':\n'
        sequence_report += seq + '\n\n'

    return sequence_report 

