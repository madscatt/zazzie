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
import numpy

import sassie.build.pdbscan.pdbscan.pdbscan_utils as pdbscan_utils
import sassie.build.pdbscan.pdbscan.data_struct as data_struct

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

    #process_segment_input(other_self, mol, pdbscan_report)

    #process_sequence_input(other_self, mol)

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
                new_segname = new_segname.encode('ascii','ignore')

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
                new_segname = new_segname.encode('ascii','ignore')

                log.info("new segname = " + new_segname)
                
                ndx = int(sassie_query_object.answer["_response"]["segment_list_box"])

                if segname_utils.valid_segname(new_segname, mol.segnames()):
                    segname_utils.rename_segment(other_self, mol, ndx, new_segname)
                   
                else:
                    log.info("valid segname = False\n")
                    log.info("valid segname = False\n")
                    log.info("new_segname = " + new_segname + "\n")
                    log.info("new_segname = " + new_segname + "\n")
                    log.info("ndx = " + str(ndx) + "\n")
                    log.info("ndx = " + str(ndx) + "\n")

        ui_output = segname_utils.get_segment_starts(other_self)

        segname_starts = json.loads(
                ui_output, object_hook = segname_utils.convert_segname_start) 

        if segname_starts:
            log.info("segname_starts has changed: updating segment information")
            segname_utils.redefine_segments(mol, segname_starts)

        mol.check_segname_simulation_preparedness()

    return

def process_sequence_input(other_self, mol):

    mvars = other_self.mvars
    log = other_self.log

    log.info('processing_sequence_input_sassie-web')

    mol.extract_sequence_info(subdiv_type='segname')       
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
                ###TODO: this should be handled by genapp

                try: 
                    sequence_segname_index = int(sassie_query_object.answer["_response"]["sequence_listbox"])
                    log.info("sequence_segname_index = " + str(sequence_segname_index))
                except:
                    log.info("user did not select a segment")
                    error = "please choose a segment for each file you wish to upload"
                    sassie_query_object.display_error(error)


                sequence_filename = sassie_query_object.answer["sequence_lrfile"]
    
                log.info("sequence_filename = " + sequence_filename[0])

                ordered_names, sequences = fasta_utils.parse_fasta_file(sequence_filename[0])

                try:

                    if len(ordered_names) == 1:
                        chosen = ordered_names[0]
                    else:
                        chosen = ordered_names[0]
                        ###TODO: add user input to pick a single sequence from list
                    fasta_sequence = fasta_utils.reformat_fasta(sequences[chosen])
                    moltype = segname_utils.get_segment_moltype(mol, seq_segnames[sequence_segname_index])
                    log.info("moltype = " + str(moltype))
                    valid_fasta = fasta_utils.validate_fasta(fasta_sequence, moltype)

                except:
                    log.info("unable to parse and validate input fasta file")
                    error = "unable to parse and validate input fasta file"
                    sassie_query_object.display_error(error)


                else:
                    success = fasta_utils.complete_sequence_fasta(mol, seq_segnames[sequence_segname_index],\
                                    fasta_sequence)

                    if not success:
                        log.info("Fasta sequence did not match sequence from chosen segment")
                        error = "Fasta sequence did not match sequence from chosen segment"
                        sassie_query_object.display_error(error)
    
#    seq_segnames = mol.segname_info.sequence.keys()
#
#    sequence_report = segname_utils.create_sequence_report(mol, seq_segnames) 

    return


def is_float(value):

    try:
        float(value)
        return True
    except:
        return False

def check_value(sassie_query_object, value):

    my_value = sassie_query_object.answer["_response"][value]

    if is_float(my_value):
        my_value = float(my_value)
        return True, my_value
    else:
        error = 'value entered for ' + value + ' needs to be a float or int'
        sassie_query_object.display_error(error)
        return False, my_value


def build_biomt_data(sassie_query_object, log):

    flag = True
    rot = None
    trans = None

    flag, a11 = check_value(sassie_query_object, 'a11')
    if not flag: return False, rot, trans
    flag, a12 = check_value(sassie_query_object, 'a12')
    if not flag: return False, rot, trans
    flag, a13 = check_value(sassie_query_object, 'a13')
    if not flag: return False, rot, trans
    flag, a21 = check_value(sassie_query_object, 'a21')
    if not flag: return False, rot, trans
    flag, a22 = check_value(sassie_query_object, 'a22')
    if not flag: return False, rot, trans
    flag, a23 = check_value(sassie_query_object, 'a23')
    if not flag: return False, rot, trans
    flag, a31 = check_value(sassie_query_object, 'a31')
    if not flag: return False, rot, trans
    flag, a32 = check_value(sassie_query_object, 'a32')
    if not flag: return False, rot, trans
    flag, a33 = check_value(sassie_query_object, 'a33')
    if not flag: return False, rot, trans
 
    flag, t1 = check_value(sassie_query_object, 't1')
    if not flag: return False, rot, trans
    flag, t2 = check_value(sassie_query_object, 't2')
    if not flag: return False, rot, trans
    flag, t3 = check_value(sassie_query_object, 't3')
    if not flag: return False, rot, trans

    log.info("biomt: a11 = " +  str(a11))
    log.info("biomt: a12 = " +  str(a12))
    log.info("biomt: a13 = " +  str(a13))
    log.info("biomt: a21 = " +  str(a21))
    log.info("biomt: a22 = " +  str(a22))
    log.info("biomt: a23 = " +  str(a23))
    log.info("biomt: a31 = " +  str(a31))
    log.info("biomt: a32 = " +  str(a32))
    log.info("biomt: a33 = " +  str(a33))

    log.info("biomt: t1 = " +  str(t1))
    log.info("biomt: t2 = " +  str(t2))
    log.info("biomt: t3 = " +  str(t3))

    rot_user = [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]]
    trans_user = [t1, t2, t3]

    check_rot, rot = biomt_utils.check_rotation(rot_user)
    if not check_rot:
        flag = False
        error = 'Rotation Matrix must be a 3 x 3 array of numeric values'
        sassie_query_object.display_error(error)    

    check_trans, trans = biomt_utils.check_translation(trans_user)
    if not check_trans:
        flag = False
        error = 'Error: Translation Vector must be a 3 x 1 array of numeric values'
        sassie_query_object.display_error(error)    

    if ((rot == numpy.identity(3)).all()) and ((trans == numpy.array([0.0, 0.0, 0.0])).all()):
        flag = False
        error = "A second identity transform will not be added"
        sassie_query_object.display_error(error)    

    return flag, rot, trans


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

                selected_segnames_ndx = map(int, sassie_query_object.answer["_response"]["biomt_listbox"])

                log.info("biomt: selected_segnames_ndx = " + ','.join(str(x) for x in selected_segnames_ndx))
                selected_segnames = []
                for i in xrange(len(segnames)):
                    if i in selected_segnames_ndx:
                        selected_segnames.append(segnames[i])

                log.info("biomt: selected_segnames = " + ','.join(selected_segnames))

                biomt_rec = biomt_utils.init_biomt(selected_segnames)

                flag = False
                while not flag:
                    sassie_query_object.process_biomt_matrix(other_self, mol)
                    flag, rot, trans = build_biomt_data(sassie_query_object, log)

                choice = sassie_query_object.answer["_response"]["button"]
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
