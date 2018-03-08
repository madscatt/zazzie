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

import numpy
import sassie.build.pdbscan.pdbscan.data_struct as data_struct

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

