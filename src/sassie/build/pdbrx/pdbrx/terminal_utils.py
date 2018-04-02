# -*- coding: utf-8 -*-
"""
Methods to manage user input from terminal

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

import json

import sassie.build.pdbscan.pdbscan.pdbscan_utils as pdbscan_utils
import sassie.build.pdbscan.pdbscan.report as report

#import sassie.build.pdbrx.pdbrx.segname_utils as segname_utils
#import sassie.build.pdbrx.pdbrx.fasta_utils as fasta_utils
#import sassie.build.pdbrx.pdbrx.biomt_utils as biomt_utils

from . import segname_utils as segname_utils
from . import fasta_utils as fasta_utils
from . import biomt_utils as biomt_utils

from . import cmdline_segname_editor as cmd_segname_edit
from . import cmdline_transform_editor

# make Python 2.x input behave as in Python 3
try:
    input = raw_input
except NameError:
    pass


def handle_terminal_user_input(other_self, mol):
    """
    Present user with options to edit segmentation, sequence and BIOMT from
    the command line.

    @return:
    """

    process_segment_input(other_self, mol)
    
    process_sequence_input(other_self, mol)
    
    process_biomt_input(other_self, mol)

    return

### segment terminal input processing methods

def process_segment_input(other_self, mol):

    accepted_segmentation = False

    print(
                "Do you wish to edit the system segmentation? (answer [y]es/[n]o)")

    while not accepted_segmentation:

        choice = input().lower()

        if choice in ['y', 'yes']:

            #ui_output = cmd_segname_edit.SegnameEditor(
            #    mol.segnames(), other_self.resid_descriptions, max_row=20).get_segment_starts()
            #ui_output = cmd_segname_edit.SegnameEditor(
            #    other_self, mol, max_row=20).get_segment_starts()
            dum = cmd_segname_edit.SegnameEditor(
                other_self, mol, max_row=20)

            ui_output = segname_utils.get_segment_starts(other_self)

            segname_starts = json.loads(
                ui_output, object_hook = segname_utils.convert_segname_start)

            if segname_starts:
                segname_utils.redefine_segments(mol, segname_starts)

            dumfile = open('dum.txt', 'a')
            dumfile.write('# BEFORE PREPAREDNESS\n')
            dumfile.write('type(mol.segnames()) \n' + str(type(mol.segnames())) + '\n')
            for value in mol.segnames():
                dumfile.write(value + '\n')
            dumfile.close()

            mol.check_segname_simulation_preparedness()

            accepted_segmentation = True

        elif choice in ['n', 'no']:

            accepted_segmentation = True

    return

### fasta terminal input processing methods

def get_user_fasta_sequence(other_self, segname, moltype):
    """
    Get FASTA sequence from user to complete a segment with missing residues

    @type segname : string
    @param segname: Segment name selected to have sequence altered
    @type moltype : string
    @param moltype: Type of polymer being edited (protein/dna/rna)
    @rtype :  string
    @return:  FASTA sequence from user
    """

    valid_fasta = False


    prompt = 'Enter filename containing FASTA sequence for segment {0:s}: '.format(
                segname)

    fasta_file = pdbscan_utils.get_command_line_filepath(prompt)

    ordered_names, sequences = fasta_utils.parse_fasta_file(fasta_file)

    if len(ordered_names) == 1:
        chosen = ordered_names[0]
    else:

        print('Choose sequence to use: \n')

        rep, options = fasta_utils.list_fasta_sequences(
                    sequences, ordered_names)

        for line in rep:
            print(line)

        user_input = -1

        while input not in options:
            user_input = input('')
            try:
                user_input = int(user_input)
            except ValueError:
                user_input = -1

        chosen = ordered_names[user_input]

    fasta_sequence = fasta_utils.reformat_fasta(sequences[chosen])

    valid_fasta = fasta_utils.validate_fasta(fasta_sequence, moltype)

    if not valid_fasta:
        fasta_sequence = ''

    return fasta_sequence

def process_sequence_input(other_self, mol):

    accepted_sequences = False
     
    seq_segnames = mol.segname_info.sequence.keys()

    while not accepted_sequences:

        print("Current sequences (lowercase indicates residues not in coordinates): ")

        for segname in seq_segnames:
            seq = mol.segname_info.sequence_to_fasta(
                    segname, missing_lower=True)
            print(segname + ':')
            print(seq)

        print("Do you want to edit any sequences? (answer [y]es/[n]o)")
        choice = input().lower()

        if choice in ['y', 'yes']:

            if len(seq_segnames) > 1:

                print("Which segment do you wish to provide a sequence for?")
                segname = input().strip()

            else:
                segname = seq_segnames[0]

            if segname in seq_segnames:

                moltype = segname_utils.get_segment_moltype(mol, segname)
                fasta_sequence = self.get_user_fasta_sequence(other_self,
                        segname, moltype)

                if fasta_sequence:

                    success = fasta_utils.complete_sequence_fasta(mol, 
                            segname, fasta_sequence)

                    if not success:

                        print(
                                "FASTA did not match existing sequence description")

                else:

                    print("Invalid FASTA sequence")

            else:

                print("Invalid segname selected")

        elif choice in ['n', 'no']:

            accepted_sequences = True

    return

### biomt terminal input processing methods

def get_biological_unit_transform(other_self, mol):
    """
    Check matrix for biological unit transform from user

    @return:
    """

    segnames_json = json.dumps(mol.segnames())

    ui_output = cmdline_transform_editor.user_biomt(segnames_json)

    user_biomt = biomt_utils.biomt_json2data(ui_output)

    if user_biomt:

        biomol_nos = mol.segname_info.biomt.keys()
        biomol_no = sorted(biomol_nos)[-1] + 1

        mol.segname_info.biomt[biomol_no] = user_biomt

    else:
        pass

    return

def process_biomt_input(other_self, mol):

    if mol.segname_info.biomt:

        print("Current biological unit transforms: ")

        for line in report.create_biomt_summary(mol.segname_info.biomt):

            print(line)

    else:
        print("There are no existing biological unit transforms")

    choice_made = False

    while not choice_made:
        print(
                "Do you want to add a new biological unit transform? (answer [y]es/[n]o)")
        choice = input().lower()

        if choice in ['y', 'yes']:

            self.get_biological_unit_transform()
            choice_made = True

            if mol.segname_info.biomt:

                print("Updated biological unit transforms: ")

                for line in report.create_biomt_summary(mol.segname_info.biomt):
                    print(line)

        elif choice in ['n', 'no']:
            choice_made = True

    return
