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

import json

import sassie.build.pdbscan.pdbscan.pdbscan_utils as pdbscan_utils

#import sassie.build.pdbrx.pdbrx.segname_utils as segname_utils
#import sassie.build.pdbrx.pdbrx.fasta_utils as fasta_utils
#import sassie.build.pdbrx.pdbrx.biomt_utils as biomt_utils

from . import segname_utils as segname_utils
from . import fasta_utils as fasta_utils
from . import biomt_utils as biomt_utils

from . import sassie_web_editor as sassie_web_editor


def handle_sassie_web_user_input(other_self, mol, pdbscan_report):

    """
    Present user with options to edit segmentation, sequence and BIOMT from
    the commandline.

    @return:
    """

    mvars = other_self.mvars
    log = other_self.log

    process_segment_input(other_self, mol, pdbscan_report)

    process_sequence_input(other_self, mol)

    process_biomt_input(other_self, mol)

    return


def process_segment_input(other_self, mol, pdbscan_report):

    mvars = other_self.mvars
    log = other_self.log

    log.info('processing_segment_input_sassie-web')
            
    sassie_query_object  = sassie_web_editor.SegnameEditor(\
                mol.segnames(), other_self.resid_descriptions, other_self.json_variables, pdbscan_report, log)
        
    choice = sassie_query_object.answer["_response"]["button"]

    if choice == "yes":

        in_loop = True

        while in_loop:

            sassie_query_object.display_and_query_segments_loop(mol, pdbscan_report)

            segment_choice = sassie_query_object.answer["_response"]["button"]

            if segment_choice == "accept":
                in_loop = False

            elif segment_choice == "split":
                pass

            elif segment_choice == "join":
                pass

            elif segment_choice == "rename":
                pass

        #segname_starts = get_user_segmentation()

        #if segname_starts:
        #    segname_utils.redefine_segments(mol, segname_starts)
        #    mol.check_segname_simulation_preparedness()


def process_sequence_input(other_self, mol):

    mvars = other_self.mvars
    log = other_self.log

    log.info('processing_sequence_input_sassie-web')
       
    seq_segnames = mol.segname_info.sequence.keys()
    sequence_report = fasta_utils.create_sequence_report(mol, seq_segnames)
            
    sassie_query_object  = sassie_web_editor.FastaEditor(\
                other_self.json_variables, sequence_report, log)
        
    choice = sassie_query_object.answer["_response"]["button"]

    if choice == "yes":

        in_loop = True

        while in_loop:

            #seq_segnames = mol.segname_info.sequence.keys()
            #sequence_report = fasta_utils.create_sequence_report(mol, seq_segnames)

            sassie_query_object.display_and_query_sequence_loop(mol, sequence_report)

            sequence_choice = sassie_query_object.answer["_response"]["button"]

            if sequence_choice == "done":
                in_loop = False

            elif sequence_choice == "submit":
                pass

#    seq_segnames = mol.segname_info.sequence.keys()
#
#    sequence_report = segname_utils.create_sequence_report(mol, seq_segnames) 

    return

def process_biomt_input(other_self, mol):

    mvars = other_self.mvars
    log = other_self.log

    log.info('processing_biomt_input_sassie-web')


    if mol.segname_info.biomt:
        print("Current biological unit transforms: ")

        import sassie.build.pdbscan.pdbscan.report as report
        for line in report.create_biomt_summary(mol.segname_info.biomt):

            print(line)
    else:
        print("There are no existing biological unit transforms")



#    sequence_report = fasta_utils.create_sequence_report(mol, seq_segnames)
            
    #sassie_query_object  = sassie_web_editor.FastaEditor(\
#                other_self.json_variables, sequence_report, log)
        
    #choice = sassie_query_object.answer["_response"]["button"]

    #if choice == "yes":
#
#        in_loop = True
#
#        while in_loop:
#

#            sassie_query_object.display_and_query_sequence_loop(mol, sequence_report)

#            sequence_choice = sassie_query_object.answer["_response"]["button"]

#            if sequence_choice == "done":
#                in_loop = False

#            elif sequence_choice == "submit":
#                pass
#

    return
