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
import yaml

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

def final_segment_selection(other_self, mol, prep_report):

    mvars = other_self.mvars
    log = other_self.log

    log.info('processing_final_segment_selection_sassie-web')

    segname_list = mol.segnames()

    choice_made = False

    while not choice_made:

        if len(segname_list) > 1:

            sassie_query_object  = sassie_web_editor.SegnameChoice(\
                other_self, segname_list, prep_report, log)            

            selected_segnames_ndx = map(int, sassie_query_object.answer["_response"]["segname_listbox"])

            log.info("selected_segnames_ndx = " + ''.join(str(x) for x in selected_segnames_ndx))
            selected_segnames = []
            for i in xrange(len(segname_list)):
                if i in selected_segnames_ndx:
                    selected_segnames.append(segname_list[i])

            log.info("selected_segnames = " + ''.join(selected_segnames))

        else:

            selected_segnames = [segname_list[0]]

        if selected_segnames:

            if not segname_utils.check_segment_status(mol, selected_segnames):
                pass
                # TODO: This is just a stub - needs filling out
                #self.process_non_ff(selected_segnames)

                #if segname_check_segment_status(mol, selected_segnames):
                #    choice_made = True

            else:
                choice_made = True

        else:
            ###TODO: evaluate whether to tell user to select at least one segment??
            print("No segnames selected")

    return selected_segnames

def process_segment_input(other_self, mol, pdbscan_report):

    mvars = other_self.mvars
    log = other_self.log

    log.info('processing_segment_input_sassie-web')
            
    sassie_query_object  = sassie_web_editor.SegnameEditor(\
                other_self,  mol.segnames(), pdbscan_report, log)
        
    choice = sassie_query_object.answer["_response"]["button"]

    if choice == "yes":

        in_loop = True

        while in_loop:

            sassie_query_object.display_and_query_segments_loop(other_self, mol, pdbscan_report)

            segment_choice = sassie_query_object.answer["_response"]["button"]

            if segment_choice == "accept":
                log.info("in segment choice = accept\n")
                in_loop = False

            elif segment_choice == "split":
                log.info("in segment choice = split\n")

                new_segname_query = sassie_query_object.query_new_segname()

                new_segname = new_segname_query["_response"]["new_segname"]

                log.info("new segname = " + new_segname)

                ndx = int(sassie_query_object.answer["_response"]["segment_list_box"])

                if segname_utils.valid_segname(new_segname, mol.segnames()):
                    segname_utils.split_segnames(other_self, mol, ndx, new_segname)

            elif segment_choice == "join":
                log.info("in segment choice = join\n")

                ndx = int(sassie_query_object.answer["_response"]["segment_list_box"])

                error = segname_utils.join_segnames(other_self, mol, ndx)

                #for k,v in new_segname_query.iteritems():
                #    if isinstance(v, dict):
                #        try:
                #            for k2,v2 in v.iteritems():
                #                 log.info("key2, value2 = " + k2 + "\t" + v2 + "\n")
                #                 log.info("type(key2), type(value2) = " + str(type(k2)) + "\t" + str(type(v2)) + "\n")
                #        except:
                #            log.info("type(key2), type(value2) = " + str(type(k2)) + "\t" + str(type(v2)) + "\n")
                #            if isinstance(v2,list):
                #                st = ''.join(v2)
                #                log.info("list = " + st)
#
#                    else: 
#                        log.info("key, value = " + k + "\t" + v + "\n")

            elif segment_choice == "rename":
                log.info("in segment choice = rename\n")

                new_segname_query = sassie_query_object.query_new_segname()

                new_segname = new_segname_query["_response"]["new_segname"]

                log.info("new segname = " + new_segname)

                ndx = int(sassie_query_object.answer["_response"]["segment_list_box"])

                if segname_utils.valid_segname(new_segname, mol.segnames()):
                    segname_utils.rename_segment(other_self, mol, ndx, new_segname)
                else:
                    log.info("valid segname = False\n")
                    log.info("new_segname = " + new_segname + "\n")
                    log.info("ndx = " + str(ndx) + "\n")

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

            seq_segnames = mol.segname_info.sequence.keys()
            sequence_report = fasta_utils.create_sequence_report(mol, seq_segnames)

            sassie_query_object.display_and_query_sequence_loop(mol, sequence_report)

            sequence_choice = sassie_query_object.answer["_response"]["button"]

            if sequence_choice == "done":
                in_loop = False

            elif sequence_choice == "submit":
           
                ###TODO: catch if submit is selected WITHOUT selecting a file to upload or use!!!
 
                sequence_segname_index = int(sassie_query_object.answer["_response"]["sequence_listbox"])

                log.info("sequence_segname_index = " + str(sequence_segname_index))

                sequence_filename = sassie_query_object.answer["sequence_lrfile"]
    
                log.info("sequence_filename = " + sequence_filename[0])

                ordered_names, sequences = fasta_utils.parse_fasta_file(sequence_filename[0])

                if len(ordered_names) == 1:
                    chosen = ordered_names[0]
                else:
                    pass
                    ###TODO: add user input to pick a single sequence from list

                fasta_sequence = fasta_utils.reformat_fasta(sequences[chosen])
                moltype = segname_utils.get_segment_moltype(mol, seq_segnames[sequence_segname_index])
                log.info("moltype = " + str(moltype))
                valid_fasta = fasta_utils.validate_fasta(fasta_sequence, moltype)

                if not valid_fasta:
                    ###TODO: deal with failure case
                    pass
                else:
                    success = fasta_utils.complete_sequence_fasta(mol, seq_segnames[sequence_segname_index],\
                                    fasta_sequence)

                    if not success:
                        log.info("Fasta seuqnce did not match sequence from chosen segment")

                #for k,v in sassie_query_object.answer.iteritems():
                #    if isinstance(v, dict):
                #        try:
                #            for k2,v2 in v.iteritems():
                #                 log.info("key2, value2 = " + k2 + "\t" + v2 + "\n")
                #                 log.info("type(key2), type(value2) = " + str(type(k2)) + "\t" + str(type(v2)) + "\n")
                #        except:
                #            log.info("type(key2), type(value2) = " + str(type(k2)) + "\t" + str(type(v2)) + "\n")
                #            if isinstance(v2,list):
                #                st = ''.join(v2)
                #                log.info("list = " + st)
#
#                    else: 
#                        log.info("key, value = " + k + "\t" + v + "\n")
#
#

#    seq_segnames = mol.segname_info.sequence.keys()
#
#    sequence_report = segname_utils.create_sequence_report(mol, seq_segnames) 

    return

def process_biomt_input(other_self, mol):

    mvars = other_self.mvars
    log = other_self.log

    log.info('processing_biomt_input_sassie-web')

    sassie_query_object  = sassie_web_editor.BiomtEditor(\
                other_self.json_variables, mol, log)
        
    choice = sassie_query_object.answer["_response"]["button"]

    if choice == "yes":

        in_loop = True

        while in_loop:

                segnames_json = json.dumps(mol.segnames())

                #ui_output = cmdline_transform_editor.user_biomt(segnames_json)

                segnames = yaml.safe_load(segnames_json)

                sassie_query_object.display_and_query_biomt_loop(mol)

                #selected_segnames = sassie_query_object.answer["_response"]["q1"]
                selected_segnames = sassie_query_object.answer["_response"]["biomt_listbox"]

                log.info("selected_segnames = " + ''.join(selected_segnames))

                sassie_query_object.process_biomt_matrix(other_self, mol)

                

#                log.info("type(sassie_query_object) = " + str(type(sassie_query_object)))
#                log.info("type(sassie_query_object.answer) = " + str(type(sassie_query_object.answer)))

#                for k,v in sassie_query_object.answer.iteritems():
#                    if isinstance(v, dict):
#                        try:
#                            for k2,v2 in v.iteritems():
#                                log.info("key2, value2 = " + k2 + "\t" + v2 + "\n")
#                                log.info("type(key2), type(value2) = " + str(type(k2)) + "\t" + str(type(v2)) + "\n")
#                        except:
##                            log.info("type(key2), type(value2) = " + str(type(k2)) + "\t" + str(type(v2)) + "\n")
#                            if isinstance(v2,list):
#                                st = ''.join(v2)
#                                log.info("list = " + st)
#
#                            
#                    else: 
#                        log.info("key, value = " + k + "\t" + v + "\n")

                choice = sassie_query_object.answer["_response"]["button"]
                if choice == "skip":
                    in_loop = False
                if choice == "submit":
                    in_loop = False

    #            biomt_rec = init_biomt(selected_segnames)

#   #             print('An identity transform already exists, enter a new transform:')
    #            sassie_query_object.add_transform(biomt_rec)
#
#                sassie_query_object.display_and_query_biomt_loop(mol, sequence_report)
#
#                biomt_choice = sassie_query_object.answer["_response"]["button"]
#
#                if biomt_choice == "done":
#                        in_loop = False
#
#                elif biomt_choice == "edit":
#                        sassie_query_object.edit_transform(biomt_rec)
#
#                elif biomt_choice == "add":
#                        sassie_query_object.add_transform(biomt_rec)

#TODO: HERE

       # return json.dumps(prepare_biomt_json(biomt_rec))
       
         #user_biomt = self.biomt_json2data(ui_output)
#
#                if user_biomt:
#
#                        biomol_nos = mol.segname_info.biomt.keys()
#                        biomol_no = sorted(biomol_nos)[-1] + 1
#
#                mol.segname_info.biomt[biomol_no] = user_biomt

    return
