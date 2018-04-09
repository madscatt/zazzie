#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Command line editor of segment divides in SasMol objects

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

from __future__ import division  

import json
import sys
from StringIO import StringIO
import numpy 

import sassie.util.communication as communication
import sassie.build.pdbscan.pdbscan.report as pdbscan_report

class BiomtEditor():
        """
        Interface to for input/ouput of questions related to biomt entries
        """
    
        def __init__(self, json_variables, mol, log):
                """
                Setup the environment and display to show residue infromation and
                editing instructions to the user

                """

                self.json_variables = json_variables
                self.log = log

                self.answer = self.ask_question_edit_biomt(mol)
   
                return 
    
        def ask_question_edit_biomt(self, mol):

                timeout = 3600

                if mol.segname_info.biomt:
                        text_string = "Current biological unit transforms: \n"
                        for line in pdbscan_report.create_biomt_summary(mol.segname_info.biomt):
                                text_string += line + '\n'
                else:
                        text_string = "There are no existing biological unit transforms\n"

                pdbscan_dict = {}
                pdbscan_dict["id"] =  "text_2"
                pdbscan_dict["type"] = "textarea"
                pdbscan_dict["default"] = text_string
                pdbscan_dict["rows"] = text_string.count('\n') + 4
                pdbscan_dict["cols"] = 80
                pdbscan_dict["fontfamily"] = "monospace"
        
                label_dict = {}
                label_dict["id"] = "label_2"
                label_dict["type"] = "label"
                label_dict["label"] = "<p><hr><br>Do you wish to add or edit biomt records?<br><br>"
        
                my_question = {}
                my_question["id"] = "q2"
                my_question["title"] = "PDB Scan Biomt Report"
                my_question["text"] = "<p>Review system biomt report below</p><br><hr><br>"
                my_question["buttons"] = ["yes", "no"]
                my_question["fields"] = [ label_dict, pdbscan_dict ]
        
                self.log.info(json.dumps(pdbscan_dict))
        
                answer = communication.tcpquestion(self.json_variables, my_question, timeout)

                return answer        
    
        def display_and_query_biomt_loop(self, mol):

            timeout = 3600

            listbox_dict = {}
            listbox_dict["id"] = "biomt_listbox"
            listbox_dict["type"] = "listbox"
            listbox_dict["fontfamily"] = "monospace"

#        my_values, my_returns = self.create_sequence_data()

            my_values = mol.segnames()
            #my_returns = [i for i in range(len(my_values))]

            listbox_dict["values"] = my_values
            #listbox_dict["returns"] = my_returns
            listbox_dict["size"] = 10 
            listbox_dict["help"] = "select row(s) and choose an option below: command-click for non-adjacent rows (Mac) or control-click (Windows)"
            listbox_dict["header"] = "choose segment(s) you wish to assign a new biomt record\n click 'done' if you are finished\n\n"
            listbox_dict["fontsize"] = "0.93em"
            listbox_dict["multiple"] = "true"

            my_question = {}
            my_question["id"] = "q1"
            my_question["title"] = "PDB Rx Biomt Editor"
            my_question["text"] = "<p>Choose Segment(s) for New Biomt Record</p><br><hr><br>"
            my_question["buttons"] = ["done"]
            my_question["fields"] = [listbox_dict]

            self.answer = communication.tcpquestion(self.json_variables, my_question, timeout);

            return 

        def select_segnames(self, segnames, log):
            '''
            Prompt user for a selection fo segment names to apply BIOMT to.
        
            @type segnames :  list
            @param segnames:  Valid segnames contained in structure
            @return:
            '''
       
            log.info("IN SELECT_SEGNAMES\n") 

            seg_list = ', '.join(segnames)
        
            print("Select list of segments for the transformation to be applied to "
                "using Python list notation e.g. ['A','B'])")
            print("Available segments are: " + seg_list)
        
            selected = ''
        
            while not selected:
        
                try:
                   # selected = eval(input())
                    if not set(selected).issubset(set(segnames)):
                        print('All values in list must be available segments')
                        selected = ''
                except:
                   # print('Invalid input, try again')
                    selected = ''
        
            return list(selected)

        def edit_transform(self, biomt_rec):
            '''
            Allow user editing of input BIOMT data

            @type biomt_rec : dict
            @param biomt_rec: Description of unit transforms to create biological unit
            @return:
            '''

            #print_biomt(biomt_rec)

            rec_no = 0
            no_entries = len(biomt_rec.rot)
        
            if no_entries > 1:
                pass        
            #    while rec_no == 0:
#
            #        print("\nSelect record to edit")
            #        txt = input()
#                    try:
#                        rec_no = int(txt)
#                        if rec_no <= 0 or rec_no >= no_entries:
#                            print('Invalid record selected')
#                            rec_no = 0
#                    except:
#                        print('Invalid record selected')
#                        rec_no = 0
#
#                input_txt = ''
#                valid_transform = False
#        
#                while not valid_transform and input_txt not in ['x', 'X']:
#        
#                    valid_transform, rot, trans = get_user_transform()
#        
#                    if valid_transform:
#                        biomt_rec['rot'][rec_no] = rot
#                        biomt_rec['trans'][rec_no] = trans
#                    else:
#                        print('Transform not edited, to return to other options '
#                            'press "X", any other key to try again.')
#                        input_txt = sys.stdin.read(1)
#        
#            else:
#        
#                print('Only the identity tranformation exists, add a new transform')
        
            return

        def add_transform(self, biomt_rec):
            '''
            Get BIOMT transformations from user and add to record

            @type biomt_rec : dict
            @param biomt_rec: Description of unit transforms to create biological unit
            @return:
            '''

            input_txt = ''
            valid_transform = False

            #while not valid_transform and input_txt not in ['x', 'X']:
##
#                valid_transform, rot, trans = get_user_transform()
#
#                if valid_transform:
#                    biomt_rec['rot'].append(rot)
#                    biomt_rec['trans'].append(trans)
#                else:
#                    print('No transform added, to return to other '
#                        'options press "X", any other key to try again.')
#                    input_txt = sys.stdin.read(1)

            return


        def get_user_transform(self):
                '''
                Get BIOMT style rotation matrix and translation vector from user.

                @return:
                '''

                flag = True
                rot = None
                trans = None

                print("Enter rotation matrix (3 x 3, using Python list notation "
                        "i.e.: [[1,0,0],[0,1,0],[0,0,1]]):")
#                rot_user = eval(input())
                check_rot, rot = biomt_utils.check_rotation(rot_user)

                if not check_rot:
                        flag = False
                        print("Rotation must be a 3 x 3 array of numeric values")

                else:
                
                        print("Enter translation vector (3 x 1, using Python list notation "
                        "i.e.: [1,0,0]):")
#                        trans_user = eval(input())
#                        check_trans, trans = biomt_utils.check_translation(trans_user)

                        if not check_trans:
                            flag = False
                            print("Error: Translation must be a 3 x 1 array of numeric values")

                if ((rot == np.identity(3)).all()) and ((trans == np.array([0.0, 0.0, 0.0])).all()):
                        print("A second identity transform will not be added")
                        flag = False

                return flag, rot, trans


        def user_biomt(self, segnames_json):
                '''
                Get user to edit BIOMT information from sassie-web 

                @type  segnames_json:  file
                @param segnames_json:  Valid segment names in JSON format
                @return:
                '''

                #segnames = yaml.safe_load(segnames_json)

                #selected_segnames = select_segnames(segnames)

                #biomt_rec = init_biomt(selected_segnames)

                print('An identity transform already exists, enter a new transform:')
                self.add_transform(biomt_rec)

                input_txt = ''

                while input_txt not in ['q', 'Q']:

                #print('Press "E" to edit an existing transform, '
                #        '"A" to add a new transform or '
                #        '"Q" to accept current BIOMT and quit')
                #input_txt = sys.stdin.read(1)

                    if input_txt in ['e', 'E']:
                        self.edit_transform(biomt_rec)
                    elif input_txt in ['a', 'A']:
                        self.add_transform(biomt_rec)

                return json.dumps(prepare_biomt_json(biomt_rec))


class FastaEditor():
    """
    Interface to for input/ouput of questions related to sequences
    """
    
    def __init__(self, json_variables, sequence_report, log):
        """
        Setup the environment and display to show residue infromation and
        editing instructions to the user

        @type segnames :  list
        @param segnames:  List of segment names input
        @type resid_descriptions : list
        @param resid_descriptions: List of tuples describing segname,
                                   first atomic index, resid, resname, chain
                                   and moltype for each residue.
        @type max_row :  int
        @param max_row:  Maximum number of rows to be displayed in terminal
        """

        self.json_variables = json_variables
        self.log = log

        #seq_segnames = mol.segname_info.sequence.keys()
        #sequence_report = self.create_sequence_report(mol, seq_segnames)

	self.answer = self.ask_question_edit_sequence(sequence_report)
   
        return 
    
    def ask_question_edit_sequence(self, sequence_report):

        timeout = 3600

        pdbscan_dict = {} 
        pdbscan_dict["id"] =  "text_2"
        pdbscan_dict["type"] = "textarea"
        pdbscan_dict["default"] = sequence_report
        pdbscan_dict["rows"] = len(sequence_report) + 6
        pdbscan_dict["rows"] = sequence_report.count('\n') + 4
        pdbscan_dict["cols"] = 180
        pdbscan_dict["fontfamily"] = "monospace"

        label_dict = {}
        label_dict["id"] = "label_2"
        label_dict["type"] = "label"
        label_dict["label"] = "<p><hr><br>Do you wish to edit sequences?<br><br>"

        my_question = {}
        my_question["id"] = "q2"
        my_question["title"] = "PDB Scan Sequence Report"
        my_question["text"] = "<p>Review system sequence report below</p><br><hr><br>"
        my_question["buttons"] = ["yes", "no"]
        my_question["fields"] = [ pdbscan_dict, label_dict ]

        self.log.info(json.dumps(pdbscan_dict))

        answer = communication.tcpquestion(self.json_variables, my_question, timeout)
        
        return answer

    def display_and_query_sequence_loop(self, mol, sequence_report):

        timeout = 3600

        listbox_dict = {}
        listbox_dict["id"] = "sequence_listbox"
        listbox_dict["type"] = "listbox"
        listbox_dict["fontfamily"] = "monospace"

#        my_values, my_returns = self.create_sequence_data()

        my_values = mol.segnames()
        my_returns = [i for i in range(len(my_values))]

        listbox_dict["values"] = my_values
        listbox_dict["returns"] = my_returns
        listbox_dict["size"] = 10 
        listbox_dict["help"] = "select a row and choose an option below"
        listbox_dict["header"] = "choose a segment you wish to replace sequence with fasta file data\n and click 'submit' to upload file or 'done' if you are finished\n\n"
        listbox_dict["fontsize"] = "0.93em"

        lrfile_dict = {}
        lrfile_dict["id"] = "sequence_lrfile"
        lrfile_dict["type"] = "lrfile"
        lrfile_dict["label"] = "select a file to upload"
        lrfile_dict["help"] = "select a fasta file to use"

#        my_values, my_returns = self.create_segment_data()

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Sequence Editor"
        my_question["text"] = "<p>Choose Fasta Sequence File Upload Optons Below</p><br><hr><br>"
        my_question["buttons"] = ["submit", "done"]
        my_question["fields"] = [listbox_dict, lrfile_dict]

        self.log.info(json.dumps(lrfile_dict))

        self.answer = communication.tcpquestion(self.json_variables, my_question, timeout);

        return 

class SegnameEditor():
    """
    Interface to allow users to split, join and rename segments in
    SasMol objects
    """

    def __init__(self, segnames, resid_descriptions, json_variables, pdbscan_report, log):
        """
        Setup the environment and display to show residue infromation and
        editing instructions to the user

        @type segnames :  list
        @param segnames:  List of segment names input
        @type resid_descriptions : list
        @param resid_descriptions: List of tuples describing segname,
                                   first atomic index, resid, resname, chain
                                   and moltype for each residue.
        @type max_row :  int
        @param max_row:  Maximum number of rows to be displayed in terminal
        """

        self.segnames = segnames
        self.resid_descriptions = resid_descriptions
        self.json_variables = json_variables
        self.log = log

        # Get initial locations of segment name changes in description list
        self.starting_breaks = numpy.where(
            self.resid_descriptions[:-1, 0] != self.resid_descriptions[1:, 0])[0]

	self.answer = self.ask_question_edit_segmentation(pdbscan_report, log)

        return 
    
    def ask_question_edit_segmentation(self, pdbscan_report, log):

        timeout = 3600

        pdbscan_dict = {} 
        pdbscan_dict["id"] =  "text_1"
        pdbscan_dict["type"] = "textarea"
        pdbscan_dict["default"] = '\n'.join(pdbscan_report)
        pdbscan_dict["rows"] = len(pdbscan_report) + 6
        pdbscan_dict["cols"] = 180
        pdbscan_dict["fontfamily"] = "monospace"

        label_dict = {}
        label_dict["id"] = "label_1"
        label_dict["type"] = "label"
        label_dict["label"] = "<p><hr><br>Do you wish to edit the segment definitions?<br><br>"

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Scan Seqment Report"
        my_question["text"] = "<p>Review system segment report below</p><br><hr><br>"
        my_question["buttons"] = ["yes", "no"]
        my_question["fields"] = [ pdbscan_dict, label_dict ]

        self.log.info(json.dumps(pdbscan_dict))


        log.info("SENDING QUESTION VIA TCP\n")
        answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        log.info("DONE SENDING QUESTION VIA TCP\n")
        
        return answer

    def create_segment_data(self):

        my_values = []
        my_returns = []

        i = 0

        for row in self.resid_descriptions:
            my_values.append('{0:7s} {1:>6} {2:>10s} {3:>6s} {4:>12s}'.format(
                            row[0], row[2], row[3], row[4], row[5]))
            my_returns.append(i)
            i += 1

        return my_values, my_returns

    def display_and_query_segments_loop(self, mol, pdbscan_report):

        timeout = 3600

        listbox_dict = {}
        listbox_dict["id"] = "segment_list_box"
        listbox_dict["type"] = "listbox"
        listbox_dict["fontfamily"] = "monospace"

        my_values, my_returns = self.create_segment_data()

        h_list = ['Segname', 'Resid', 'Resname', 'Chain', 'Moltype']
        header = '{0:<7s} {1:>7s} {2:>10s}   {3:<8s} {4:<8s}'.format(h_list[0], h_list[1], h_list[2], h_list[3], h_list[4])

        listbox_dict["values"] = my_values
        listbox_dict["returns"] = my_returns
        listbox_dict["size"] = 10 
        listbox_dict["help"] = "select a row and choose an option below"
        listbox_dict["header"] = header                          
        listbox_dict["fontsize"] = "0.93em"

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Seqment Editor"
        my_question["text"] = "<p>Edit Segments Definitions Below</p><br><hr><br>"
        my_question["buttons"] = ["split", "join", "rename", "accept"]
        my_question["fields"] = [listbox_dict ]

        self.log.info(json.dumps(listbox_dict))

        self.answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        
        return 
   
    
