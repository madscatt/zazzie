# -*- coding: utf-8 -*-
"""
Methods to process information in fasta format

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

import re
from textwrap import TextWrapper
import sassie.build.pdbscan.pdbscan.pdbscan_utils as utils

def create_sequence_report(mol, seq_segnames):

        sequence_report = 'Current residue sequences for each segment (uppercase letters have coordinates while lowercase letters do not have coordinates): \n\n'

        for segname in seq_segnames:
            seq = mol.segname_info.sequence_to_fasta(
                    segname, missing_lower=True)
            sequence_report += 'segname ' + segname + ':\n'
            sequence_report += seq + '\n\n'

        return sequence_report


def parse_fasta_file(filename):
    """
    Parse FASTA files. Adapted from:
    http://www.petercollingridge.co.uk/python-bioinformatics-tools/fasta-parser

    @type  filename: string
    @param filename: Path to FASTA file to be read.
    @rtype:   list, dictionary
    @return:  1. Sequence identifiers found in order they appear in the file

              2. Dictionary where keys are sequence identifiers and values are
                sequences as strings (single letter amino acid codes used).
    """

    # We are going to store ordered sequence names and the sequences
    # themselves
    order = []
    sequences = {}

    with open(filename) as fasta_file:
        for line in fasta_file:
            # Title line for sequences in FASTA files start with >
            if line.startswith('>'):
                name = line[1:].rstrip('\n')
                name = name.replace('_', ' ')
                order.append(name)
                sequences[name] = ''
            elif len(order) > 0:
                # If we have seen a title but not in this line
                # add the contents of the line to the currently named sequence
                # Note: * = chain ending character in FASTA so is removed
                sequences[name] += line.rstrip('\n').rstrip('*')

    return order, sequences

def list_fasta_sequences(sequences, ordered_names):
    """
    Create text to present input FASTA sequences as options for user
    selection and list of numbered options.

    @type sequences : dict
    @param sequences: FASTA sequences, keys = sequence name
    @type ordered_names : list
    @param ordered_names: Names of FASTA sequences ordered as they were
                              should be presented as options
    @rtype : list, list
    @return: List of strings containing FASTA sequence formatted for output
                 List of tuples of option number and sequence name
    """

    options = enumerate(ordered_names)

    rep = []

    for num, name in options:

        rep.append('(' + str(num) + ') ' + name)

        wrapper = TextWrapper(width=50)

        seq = sequences[name]

        for line in wrapper.wrap(seq):
            rep.append(line)

    return rep, options


def reformat_fasta(fasta):
    """
    Reformat input FASTA string to a single line and upper case characters
    @type fasta :  string
    @param fasta:  FASTA sequence (possibly multiple lines)
    @rtype :       string
    @return:       Single line, all caps FASTA sequence
    """

    edited_fasta = fasta.upper()

    return edited_fasta.replace('\n', '')

def validate_fasta(fasta, moltype):
    """
    Check fasta sequence is valid (i.e. contains only letters corresponding
    to monomer units appropriate to the moltype).

    @type fasta :  string
    @param fasta:  FASTA sequence to be validated
    @type moltype : string
    @param moltype: Type of polymer being edited (protein/dna/rna)
    @return:
    """

    valid = False

    if moltype == 'protein':
        valid_res = 'GALMFWKQESPVICYHRNDTX'
    else:
        valid_res = 'ACGTUX'

    allowed = [x for x in fasta if x in valid_res]

    if (len(allowed) == len(fasta)) and len(fasta) > 0:
        valid = True

    return valid

def match_fasta_model(mol, segname, new_fasta):
    """
    Check FASTA sequence from user to ensure it makes sense with the known
    coordinate sequence

    @type  segname  :  string
    @param segname  :  Name of the segment to which the FASTA sequence
                           is to be matched.
    @type  new_fasta:  string
    @param new_fasta:  FASTA format protein/nucleic acid sequence
    @rtype :  MatchObject
    @return:  Details of the matching region of the input fasta and
                  existing sequence of the selected segment
    """

    segname_info = mol.segname_info

    # The current sequence here contains '.' for all missing residues
    # (i.e. with no coordinates - these may be internal gaps or
    # terminal residues from header
    # Note: at present at least the same length of terminal residues must
    # be provided
    current_fasta = segname_info.sequence_to_fasta(
            segname, for_matching=True)

    match = re.search(current_fasta, new_fasta)

    return match

def complete_sequence_fasta(mol, segname, new_fasta):
    """
    Combine the input FASTA sequence with that existing in the coordinates
    of the selected segment.

    @type  segname  :  string
    @param segname  :  Name of the segment to which the FASTA sequence
                           is to be matched.
    @type  new_fasta:  string
    @param new_fasta:  FASTA format protein/nucleic acid sequence
    @rtype :           boolean
    @return:           Was the sequence correctly inserted into the


    @todo: Check what to do in the case of sequence gaps
    """

    segname_info = mol.segname_info
    model_no = mol.model_no

        # Sequence as a list of (resid, resname) tuples
    seq = segname_info.sequence[segname]

    first_coor_resid, first_coor_resname = segname_info.get_first_coor_resid(
            segname, model_no=model_no)

    if segname not in segname_info.missing_resids[model_no]:

        segname_info.missing_resids[model_no][segname] = {}

    missing_resids = segname_info.missing_resids[model_no][segname]

    # Match object will contain location of the coordinate sequence
    # relative to the start of the input FASTA sequence
    match = match_fasta_model(mol, segname, new_fasta)

    inserted = False

    if match:

        match_start = match.start()

        # new_fasta has residues before the existing sequence
        if match.start() > 0:

            position = 0

            first_resid = first_coor_resid - match_start

            # Resid 0 cannot exist
            if first_resid <= 0:
                first_resid = first_resid - 1

            while position < match_start:

                resid = first_resid + position
                resname = utils.conv_aa1to3(new_fasta[position])

                segname_info.prepend_residue_to_sequence(
                        segname, resid, resname)
                missing_resids[resid] = resname

                position += 1

                # Skip resid = 0
                if position == 0:
                    position += 1

        else:

            position = 0
            #### NOT TESTED
            first_resid = first_coor_resid

        # Deal with residues between the start and finish of the
        # input coordinates

        coor_start_ndx = seq.index((first_coor_resid, first_coor_resname))

        for ndx in range(coor_start_ndx, len(seq)):

            res = seq[ndx]

            resid = first_resid + position

            if res[1] == '':

                new_resname = utils.conv_aa1to3(new_fasta[position])

                seq[ndx] = (resid, new_resname)

                missing_resids[resid] = new_resname

            position += 1

        # new_fasta sequence has residues after those in existing sequence

        if position < len(new_fasta):

            for i in range(position, len(new_fasta)):

                resid = first_resid + i

                resname = utils.conv_aa1to3(new_fasta[i])

                segname_info.add_residue_to_sequence(
                        segname, resid, resname)
                segname_info.add_missing_resid(
                        segname, resid, resname, model_no=model_no)
                missing_resids[resid] = resname

    inserted = True

    return inserted

