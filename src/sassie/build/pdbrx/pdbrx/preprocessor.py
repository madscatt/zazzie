# -*- coding: utf-8 -*-
"""
Preprocessor to finalize system description after PDB Scan, this is the first
step in PDB Rx

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

import logging
import re
import json
from textwrap import TextWrapper

import yaml
import numpy as np
import sassie.build.pdbscan.pdbscan.data_struct as data_struct
import sassie.build.pdbscan.pdbscan.pdbscan_utils as utils

from . import cmdline_segname_editor as cmd_segname_edit
from . import cmdline_transform_editor

from . import sassie_web_segname_editor as sassie_web_segname_editor

# make Python 2.x input behave as in Python 3
try:
    input = raw_input
except NameError:
    pass


class PreProcessor(object):
    """
    Preprocessor checks and edits an input SasMol ready for model building
    """

    def __init__(self, *args, **kwargs):
        """
        Setup logging and atrributes for preprocessing

        @type mol   :  SasMol
        @keyword mol:  Molecule data
        @type ui_type    :  str
        @keyword ui_type :  Choice of UI type
        """
       
        if 'logger' in kwargs:
            self.logger = kwargs['logger']
        else:
            #TODO: the following has not been verified to work
            self.logger = logging.getLogger(__name__)

        if 'mol' in kwargs:
            self.mol = kwargs['mol']
        else:
            # TODO: Make a better decision on what to do here
            self.mol = None

        if 'ui' in kwargs:
            self.ui_type = kwargs['ui']
        else:
            # TODO: Not strictly needed as it is now passed from gui_mimic
            self.ui_type = 'terminal'

        if 'json' in kwargs:
            self.json = kwargs['json']
        else:
            self.json = False

        if 'default_subs' in kwargs:
            self.default_subs = kwargs['default_subs']
        else:
            self.default_subs = False

        self.logger.info('ui_type = ' + self.ui_type)
        self.logger.info('ui_type = ' + self.ui_type)
        self.logger.info('ui_type = ' + self.ui_type)

        return

    def convert_segname_start(self, data):
        """
        Convert JSON segname_start dictionary to have integer keys and ascii
        strings (as opposed to unicode)

        @rtype :  dictionary
        @return:  Keys = atom index of start of segment as integer,
                  Value = segname
        """

        return dict((int(x[0]), x[1].encode('ascii')) for x in data.items())

    def get_user_segmentation(self):
        """
        Get new segmentation from user
        @rtype :  dictionary
        @return:  Keys = atom index of start of segment,
                  Value = segname
        """

        mol = self.mol

        self.resid_descriptions = self.create_residue_descriptions()

        if self.ui_type == 'terminal':

            ui_output = cmd_segname_edit.SegnameEditor(
                mol.segnames(), self.resid_descriptions, max_row=20).get_segment_starts()

        else:
            ui_output = sassie_web_segname_editor.SegnameEditor(
                mol.segnames(), self.resid_descriptions, self.json).get_segment_starts()

            self.logger.info('SHOULD NOT GET HERE')
            self.logger.info('SHOULD NOT GET HERE')
            self.logger.info('SHOULD NOT GET HERE')
            self.logger.info('SHOULD NOT GET HERE')
            # TODO: something for a real GUI
            ui_output = []

        if ui_output:
            segname_starts = json.loads(
                ui_output, object_hook=self.convert_segname_start)
        else:
            segname_starts = {}

        return segname_starts

    def redefine_segments(self, segname_starts):
        """
        Alter segmentation in mol to follow user input then scan to obtain
        information for further processing.

        @type  segname_starts:  dictionary
        @param segname_starts:  Keys = atom index of start of segment,
                                Value = segname
        """

        self.update_segments(segname_starts)

        # Blank the previous segname_info
        self.mol.segname_info = data_struct.Info(scan_type='segname')

        # Get updated sequence/missing residue etc. information
        self.mol.segment_scan(initialize=False)

        return

    def update_segments(self, segname_starts):
        """
        Alter segmentation in mol to follow user input.

        @type  segname_starts:  dictionary
        @param segname_starts:  Keys = atom index of start of segment,
                                Value = segname
        """

        mol = self.mol

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

    def parse_fasta_file(self, filename):
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

    def list_fasta_sequences(self, sequences, ordered_names):
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

    def get_segment_moltype(self, segname):
        """
        Get the Moltype of the chosen segment

        @type segname :  str
        @param segname:  Name of segment for which moltype is required
        @rtype : str
        @return: Moltype of segmanent
        """

        mol = self.mol

        mapping = set(zip(mol.segname(), mol.moltype()))
        moltypes = [x[1] for x in mapping if x[0] == segname]

        moltype = '/'.join(moltypes)

        return moltype

    def get_user_fasta_sequence(self, segname, moltype):
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

        if self.ui_type == 'terminal':

            prompt = 'Enter filename containing FASTA sequence for segment {0:s}: '.format(
                segname)

            fasta_file = utils.get_command_line_filepath(prompt)

            ordered_names, sequences = self.parse_fasta_file(fasta_file)

            if len(ordered_names) == 1:
                chosen = ordered_names[0]
            else:

                print('Choose sequence to use: \n')

                rep, options = self.list_fasta_sequences(
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

            fasta_sequence = self.reformat_fasta(sequences[chosen])

            valid_fasta = self.validate_fasta(fasta_sequence, moltype)

        else:

            pass

        if not valid_fasta:
            fasta_sequence = ''

        return fasta_sequence

    def reformat_fasta(self, fasta):
        """
        Reformat input FASTA string to a single line and upper case characters
        @type fasta :  string
        @param fasta:  FASTA sequence (possibly multiple lines)
        @rtype :       string
        @return:       Single line, all caps FASTA sequence
        """

        edited_fasta = fasta.upper()

        return edited_fasta.replace('\n', '')

    def validate_fasta(self, fasta, moltype):
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

    def match_fasta_model(self, segname, new_fasta):
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

        segname_info = self.mol.segname_info

        # The current sequence here contains '.' for all missing residues
        # (i.e. with no coordinates - these may be internal gaps or
        # terminal residues from header
        # Note: at present at least the same length of terminal residues must
        # be provided
        current_fasta = segname_info.sequence_to_fasta(
            segname, for_matching=True)

        match = re.search(current_fasta, new_fasta)

        return match

    def complete_sequence_fasta(self, segname, new_fasta):
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

        segname_info = self.mol.segname_info
        model_no = self.mol.model_no

        # Sequence as a list of (resid, resname) tuples
        seq = segname_info.sequence[segname]

        first_coor_resid, first_coor_resname = segname_info.get_first_coor_resid(
            segname, model_no=model_no)

        if segname not in segname_info.missing_resids[model_no]:

            segname_info.missing_resids[model_no][segname] = {}

        missing_resids = segname_info.missing_resids[model_no][segname]

        # Match object will contain location of the coordinate sequence
        # relative to the start of the input FASTA sequence
        match = self.match_fasta_model(segname, new_fasta)

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

    def get_biological_unit_transform(self):
        """
        Check matrix for biological unit transform from user

        @return:
        """

        segnames_json = json.dumps(self.mol.segnames())

        if self.ui_type == 'terminal':

            ui_output = cmdline_transform_editor.user_biomt(segnames_json)

        else:

            ui_output = []

        user_biomt = self.biomt_json2data(ui_output)

        if user_biomt:

            biomol_nos = self.mol.segname_info.biomt.keys()
            biomol_no = sorted(biomol_nos)[-1] + 1

            self.mol.segname_info.biomt[biomol_no] = user_biomt

        else:
            pass

        return

    def check_biological_unit(self, biomt_unit):
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

    def biomt_json2data(self, json_biomt_rec):
        """
        Convert JSON format BIOMT records into numpy arrays for use in
        coordinate transforms

        @type json_biomt_rec :  str
        @param json_biomt_rec:  BIOMT records in JSON format
        @rtype :
        @return:
        """

        biomt_rec = yaml.safe_load(json_biomt_rec)
        valid = self.check_biological_unit(biomt_rec)

        if valid:

            for i in range(len(biomt_rec['rot'])):
                biomt_rec['rot'][i] = np.array(biomt_rec['rot'][i])
                biomt_rec['trans'][i] = np.array(biomt_rec['trans'][i])

        else:

            biomt_rec = None

        return biomt_rec

    def create_residue_descriptions(self):
        """
        Filter information from self.mol to provide residue descriptions for
        used in displaying structure contents for segment editing.

        @rtype : list
        @return: List of tuples describing segname, first atomic index,
                 resid, resname, chain and moltype for each residue
        """

        segnames = self.mol.segname()
        indices = self.mol.index()
        resids = self.mol.resid()
        resnames = self.mol.resname()
        chains = self.mol.chain()
        moltypes = self.mol.moltype()

        data = np.array(zip(segnames, indices, resids,
                            resnames, chains, moltypes))

        mask = np.where(data[:-1, 2] != data[1:, 2])[0]
        mask = np.append(mask, [len(data) - 1])

        residue_descriptions = data[mask]

        return residue_descriptions



    def handle_sassie_web_edit_options(self):

        """
        Present user with options to edit segmentation, sequence and BIOMT from
        the commandline.

        @return:
        """

        mol = self.mol
        self.resid_descriptions = self.create_residue_descriptions()

        accepted_segmentation = False

        self.logger.info('UI_TYPE = ' + self.ui_type)
        self.logger.info('UI_TYPE = ' + self.ui_type)
        self.logger.info('UI_TYPE = ' + self.ui_type)

        self.logger.info('GETTING READY TO CHAT WITH SASSIE WEB') 
        self.logger.info('GETTING READY TO CHAT WITH SASSIE WEB') 
        self.logger.info('GETTING READY TO CHAT WITH SASSIE WEB') 
            
        se = sassie_web_segname_editor.SegnameEditor(\
                mol.segnames(), self.resid_descriptions, self.json)
            
         
        self.logger.info('TYPE SE.ANSWER = ' + str(type(se.answer)))
        self.logger.info('TYPE SE.ANSWER = ' + str(type(se.answer)))
        self.logger.info('TYPE SE.ANSWER = ' + str(type(se.answer)))
           
        choice = json.dumps(se.answer['_response'])
            
        self.logger.info('SASSIE_WEB CHOICE = ' + choice) 
        self.logger.info('SASSIE_WEB CHOICE = ' + choice) 
        self.logger.info('SASSIE_WEB CHOICE = ' + choice) 

        if choice == 'yes':
            choice = 'no'

        self.logger.info('SASSIE_WEB CHOICE = ' + choice) 
        self.logger.info('SASSIE_WEB CHOICE = ' + choice) 
        self.logger.info('SASSIE_WEB CHOICE = ' + choice) 

        accepted_segmentation = True

        while not accepted_segmentation:

            if choice in ['y', 'yes']:

                segname_starts = self.get_user_segmentation()

                if segname_starts:
                    self.redefine_segments(segname_starts)
                    mol.check_segname_simulation_preparedness()

                accepted_segmentation = True

            elif choice in ['n', 'no']:

                accepted_segmentation = True

        accepted_sequences = False
        accepted_sequences = True
     
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

                    moltype = self.get_segment_moltype(segname)
                    fasta_sequence = self.get_user_fasta_sequence(
                        segname, moltype)

                    if fasta_sequence:

                        success = self.complete_sequence_fasta(
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

        if mol.segname_info.biomt:

            print("Current biological unit transforms: ")

            import sassie.build.pdbscan.pdbscan.report as report
            for line in report.create_biomt_summary(mol.segname_info.biomt):

                print(line)
        else:
            print("There are no existing biological unit transforms")

        choice_made = False
        choice_made = True

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

    def handle_terminal_edit_options(self):
        """
        Present user with options to edit segmentation, sequence and BIOMT from
        the commandline.

        @return:
        """

        mol = self.mol
        self.resid_descriptions = self.create_residue_descriptions()

        accepted_segmentation = False

        print(
                "Do you wish to edit the system segmentation? (answer [y]es/[n]o)")

        while not accepted_segmentation:

            choice = input().lower()

            if choice in ['y', 'yes']:

                segname_starts = self.get_user_segmentation()

                if segname_starts:
                    self.redefine_segments(segname_starts)
                    mol.check_segname_simulation_preparedness()

                accepted_segmentation = True

            elif choice in ['n', 'no']:

                accepted_segmentation = True

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

                    moltype = self.get_segment_moltype(segname)
                    fasta_sequence = self.get_user_fasta_sequence(
                        segname, moltype)

                    if fasta_sequence:

                        success = self.complete_sequence_fasta(
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

        if mol.segname_info.biomt:

            print("Current biological unit transforms: ")

            import sassie.build.pdbscan.pdbscan.report as report
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

    def user_edit_options(self):
        """
        Get user input from terminal or other source. ONLY TERMINAL CURRENTLY

        @return:
        """

        if self.ui_type == 'terminal':

            self.handle_terminal_edit_options()

        elif self.ui_type == 'sassie_web':

            self.handle_sassie_web_edit_options()

        return
