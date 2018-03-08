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

import json

import sassie.build.pdbscan.pdbscan.pdbscan_utils as utils

#import sassie.build.pdbrx.pdbrx.segname_utils as segname_utils
#import sassie.build.pdbrx.pdbrx.fasta_utils as fasta_utils
#import sassie.build.pdbrx.pdbrx.biomt_utils as biomt_utils

from . import segname_utils as segname_utils
from . import fasta_utils as fasta_utils
from . import biomt_utils as biomt_utils

from . import cmdline_segname_editor as cmd_segname_edit
from . import cmdline_transform_editor

from . import sassie_web_editor as sassie_web_editor

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

    def get_user_segmentation(self):
        """
        Get new segmentation from user
        @rtype :  dictionary
        @return:  Keys = atom index of start of segment,
                  Value = segname
        """

        mol = self.mol

        self.resid_descriptions = segname_utils.create_residue_descriptions(mol)

        if self.ui_type == 'terminal':

            ui_output = cmd_segname_edit.SegnameEditor(
                mol.segnames(), self.resid_descriptions, max_row=20).get_segment_starts()

        else:
            ui_output = sassie_web_editor.SegnameEditor(
                mol.segnames(), self.resid_descriptions, self.json, self.logger).get_segment_starts()

            self.logger.info('SHOULD NOT GET HERE')
            self.logger.info('SHOULD NOT GET HERE')
            self.logger.info('SHOULD NOT GET HERE')
            self.logger.info('SHOULD NOT GET HERE')
            # TODO: something for a real GUI
            ui_output = []

        if ui_output:
            segname_starts = json.loads(
                ui_output, object_hook=segname_utils.convert_segname_start)
        else:
            segname_starts = {}

        return segname_starts


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

        else:

            pass

        if not valid_fasta:
            fasta_sequence = ''

        return fasta_sequence

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

        user_biomt = biomt_utils.biomt_json2data(ui_output)

        if user_biomt:

            biomol_nos = self.mol.segname_info.biomt.keys()
            biomol_no = sorted(biomol_nos)[-1] + 1

            self.mol.segname_info.biomt[biomol_no] = user_biomt

        else:
            pass

        return

    def handle_sassie_web_edit_options(self, pdbscan_report):

        """
        Present user with options to edit segmentation, sequence and BIOMT from
        the commandline.

        @return:
        """

        mol = self.mol
        self.resid_descriptions = segname_utils.create_residue_descriptions(mol)

        accepted_segmentation = False

        self.logger.info('UI_TYPE = ' + self.ui_type)
            
        sassie_query_object  = sassie_web_editor.SegnameEditor(\
                mol.segnames(), self.resid_descriptions, self.json, pdbscan_report, self.logger)
        
        choice = sassie_query_object.answer["_response"]["button"]
            
        if choice == 'yes':
            choice = 'no'

        if choice == 'no':
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
       
        seq_segnames = mol.segname_info.sequence.keys()

        sequence_report = self.create_sequence_report(mol, seq_segnames) 
        sassie_query_object = sassie_web_editor.FastaEditor(self.json, sequence_report, self.logger) 

        choice = sassie_query_object.answer["_response"]["button"]

        self.logger.info('SEQUENCE EDIT CHOICE = ' + choice)
        
        if choice == "no":
            accepted_sequences = True

        while not accepted_sequences:

            sequence_report = self.create_sequence_report(mol, seq_segnames) 
            sassie_query_object = sassie_web_editor.FastaEditor(self.json, sequence_report, self.logger) 

            choice = sassie_query_object.answer["_response"]["button"]

            self.logger.info('SEQUENCE CHOICE = ' + choice)

#            print("Do you want to edit any sequences? (answer [y]es/[n]o)")

            choice = "no"

            if choice in ['y', 'yes']:

                if len(seq_segnames) > 1:

                    print("Which segment do you wish to provide a sequence for?")
            #        segname = input().strip()

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
        #    choice = input().lower()

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
        self.resid_descriptions = segname_utils.create_residue_descriptions(mol)

        accepted_segmentation = False

        print(
                "Do you wish to edit the system segmentation? (answer [y]es/[n]o)")

        while not accepted_segmentation:

            choice = input().lower()

            if choice in ['y', 'yes']:

                segname_starts = self.get_user_segmentation()

                if segname_starts:
                    segname_utils.redefine_segments(mol, segname_starts)
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

                    moltype = segname_utils.get_segment_moltype(mol, segname)
                    fasta_sequence = self.get_user_fasta_sequence(
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

    def user_edit_options(self, pdbscan_report):
        """
        Get user input from terminal or other source. ONLY TERMINAL CURRENTLY

        @return:
        """

        if self.ui_type == 'terminal':

            self.handle_terminal_edit_options()

        elif self.ui_type == 'sassie_web':

            self.handle_sassie_web_edit_options(pdbscan_report)
            self.logger.info('I SHOULD QUIT HERE')
            self.logger.info('I SHOULD QUIT HERE')
            self.logger.info('I SHOULD QUIT HERE')
            import sys ; sys.exit()

        return
