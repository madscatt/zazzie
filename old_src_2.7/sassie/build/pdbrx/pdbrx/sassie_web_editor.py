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

class BiomtChoice():
    """
    Interface to allow users to choose final biomt options
    """

    def __init__(self, other_self, biomt_list, biomol_report):
        """
        Setup the environment and display to show final biomt 
        selection instructions to the user

        """
        #self.log = log
        self.json_variables = other_self.json_variables

	self.answer = self.select_biomt(biomt_list, biomol_report)

        return 

    def select_biomt(self, available_biomt, biomol_report):

        timeout = 3600

        biomt_list = [str(x) for x in available_biomt]


        text_string = "Current biological unit transforms: \n"
        for line in biomol_report:
            text_string += line + '\n'

        report_dict = {}
        report_dict["id"] =  "text_1"
        report_dict["type"] = "textarea"
        report_dict["default"] = text_string
        report_dict["rows"] = text_string.count('\n') + 4
        report_dict["cols"] = 80
        report_dict["fontfamily"] = "monospace"

        listbox_dict = {}
        listbox_dict["id"] = "biomt_listbox"
        listbox_dict["type"] = "listbox"
        listbox_dict["fontfamily"] = "monospace"

        listbox_dict["values"] = biomt_list
        listbox_dict["size"] = 10
        listbox_dict["help"] = "select row(s) option below: command-click for non-adjacent rows (Mac) or control-click (Windows)"
        listbox_dict["header"] = "choose biomts(s) you wish to use and click 'submit'"
        listbox_dict["fontsize"] = "0.93em"
        listbox_dict["multiple"] = "true"

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Biomt Selection"
        my_question["text"] = "<p>Choose Biomts to Apply to Model</p><br><hr><br>"
        my_question["buttons"] = ["submit"]
        my_question["fields"] = [report_dict, listbox_dict]

        answer = communication.tcpquestion(self.json_variables, my_question, timeout)

        return answer
 
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
    
        def display_error(self, error):

            timeout = 3600

            error_dict = {} 
            error_dict["id"] =  "text_2"
            error_dict["type"] = "textarea"
            error_dict["default"] = error
            error_dict["rows"] = len(error) + 6
            error_dict["rows"] = error.count('\n') + 4
            error_dict["cols"] = 180
            error_dict["fontfamily"] = "monospace"

            my_question = {}
            my_question["id"] = "q1"
            my_question["title"] = "PDB Rx Error"
            my_question["text"] = "<p>Error Encountered: </p><br><hr><br>"
            my_question["buttons"] = ["continue"]
            my_question["icon"] = ["warning.png"]
            my_question["fields"] = [error_dict]

            answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        
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
            listbox_dict["header"] = "choose segment(s) you wish to assign a new biomt record\n click 'submit' if you are finished or 'skip' to exit without altering existing biomt record(s)\n\n"
            listbox_dict["fontsize"] = "0.93em"
            listbox_dict["multiple"] = "true"

            my_question = {}
            my_question["id"] = "q1"
            my_question["title"] = "PDB Rx Biomt Editor"
            my_question["text"] = "<p>Choose Segment(s) for New Biomt Record</p><br><hr><br>"
            my_question["buttons"] = ["submit", "skip"]
            my_question["fields"] = [listbox_dict]

            self.answer = communication.tcpquestion(self.json_variables, my_question, timeout);

            return 


        def process_biomt_matrix(self, other_self, mol):

            timeout = 3600

            mvars = other_self.mvars
            log = other_self.log

            log.info('processing_biomt_matrix_input_sassie-web')

            label_1_dict = {}
            label_1_dict["id"] = "label_1"
            label_1_dict["type"] = "label"
            label_1_dict["label"] = "enter elements of rotation matrix"

            a11_dict = {}
            a11_dict["id"] = "a11"
            a11_dict["type"] = "text"
            a11_dict["help"] = "enter value"
            a11_dict["norow"] = "true"
            a11_dict["label"] = "a11"
            a11_dict["default"] = "1"

            a12_dict = {}
            a12_dict["id"] = "a12"
            a12_dict["type"] = "text"
            a12_dict["help"] = "enter value"
            a12_dict["norow"] = "true"
            a12_dict["label"] = "a12"
            a12_dict["default"] = "0"

            a13_dict = {}
            a13_dict["id"] = "a13"
            a13_dict["type"] = "text"
            a13_dict["help"] = "enter value"
            a13_dict["norow"] = "false"
            a13_dict["label"] = "a13"
            a13_dict["default"] = "0"

            a21_dict = {}
            a21_dict["id"] = "a21"
            a21_dict["type"] = "text"
            a21_dict["help"] = "enter value"
            a21_dict["norow"] = "true"
            a21_dict["label"] = "a21"
            a21_dict["default"] = "0"

            a22_dict = {}
            a22_dict["id"] = "a22"
            a22_dict["type"] = "text"
            a22_dict["help"] = "enter value"
            a22_dict["norow"] = "true"
            a22_dict["label"] = "a22"
            a22_dict["default"] = "1"

            a23_dict = {}
            a23_dict["id"] = "a23"
            a23_dict["type"] = "text"
            a23_dict["help"] = "enter value"
            a23_dict["norow"] = "false"
            a23_dict["label"] = "a23"
            a23_dict["default"] = "0"

            a31_dict = {}
            a31_dict["id"] = "a31"
            a31_dict["type"] = "text"
            a31_dict["help"] = "enter value"
            a31_dict["norow"] = "true"
            a31_dict["label"] = "a31"
            a31_dict["default"] = "0"

            a32_dict = {}
            a32_dict["id"] = "a32"
            a32_dict["type"] = "text"
            a32_dict["help"] = "enter value"
            a32_dict["norow"] = "true"
            a32_dict["label"] = "a32"
            a32_dict["default"] = "0"

            a33_dict = {}
            a33_dict["id"] = "a33"
            a33_dict["type"] = "text"
            a33_dict["help"] = "enter value"
            a33_dict["norow"] = "false"
            a33_dict["label"] = "a33"
            a33_dict["default"] = "1"

            label_2_dict = {}
            label_2_dict["id"] = "label_2"
            label_2_dict["type"] = "label"
            label_2_dict["label"] = "enter elements of translation vector"

            t1_dict = {}
            t1_dict["id"] = "t1"
            t1_dict["type"] = "text"
            t1_dict["help"] = "enter value"
            t1_dict["norow"] = "true"
            t1_dict["label"] = "t1"
            t1_dict["default"] = "0"

            t2_dict = {}
            t2_dict["id"] = "t2"
            t2_dict["type"] = "text"
            t2_dict["help"] = "enter value"
            t2_dict["norow"] = "true"
            t2_dict["label"] = "t2"
            t2_dict["default"] = "0"

            t3_dict = {}
            t3_dict["id"] = "t3"
            t3_dict["type"] = "text"
            t3_dict["help"] = "enter value"
            t3_dict["norow"] = "false"
            t3_dict["label"] = "t3"
            t3_dict["default"] = "0"

            my_question = {}
            my_question["id"] = "q2"
            my_question["title"] = "PDB Rx Biomt Editor"
            my_question["text"] = "<p>Enter values for Biomt Matrix</p><br><hr><br>"
            my_question["buttons"] = ["submit"]
            my_question["fields"] = [label_1_dict, a11_dict, a12_dict, a13_dict, 
                                     a21_dict, a22_dict, a23_dict,
                                     a31_dict, a32_dict, a33_dict,
                                     label_2_dict, t1_dict, t2_dict, t3_dict]

            self.answer = communication.tcpquestion(self.json_variables, my_question, timeout);

            return

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

    def display_error(self, error):

        timeout = 3600

        error_dict = {} 
        error_dict["id"] =  "text_2"
        error_dict["type"] = "textarea"
        error_dict["default"] = error
        error_dict["rows"] = len(error) + 6
        error_dict["rows"] = error.count('\n') + 4
        error_dict["cols"] = 180
        error_dict["fontfamily"] = "monospace"

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Error"
        my_question["text"] = "<p>Error Encountered: </p><br><hr><br>"
        my_question["icon"] = ["warning.png"]
        my_question["buttons"] = ["continue"]
        my_question["fields"] = [error_dict]

        answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        
        return 

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
        listbox_dict["required"] = "true"

        lrfile_dict = {}
        lrfile_dict["id"] = "sequence_lrfile"
        lrfile_dict["type"] = "file"
        lrfile_dict["label"] = "select a file to upload"
        lrfile_dict["help"] = "select a fasta file to use"
        lrfile_dict["required"] = "true"

#        my_values, my_returns = self.create_segment_data()

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Sequence Editor"
        my_question["text"] = "<p>Choose Fasta Sequence File Upload Optons Below</p><br><hr><br>"
        button_1 = {}
        button_1["id"] = "submit"
        button_1["label"] = "submit"
       
        button_2 = {}
        button_2["id"] = "done"
        button_2["label"] = "done"
        button_2["skiprequired"] = "true"

      #  my_question["buttons"] = ["submit", "done"]
        my_question["buttons"] = [button_1, button_2]
        my_question["fields"] = [listbox_dict, lrfile_dict]

        self.log.info(json.dumps(listbox_dict))
        self.log.info(json.dumps(lrfile_dict))

        self.answer = communication.tcpquestion(self.json_variables, my_question, timeout);

        return 

class SegnameChoice():
    """
    Interface to allow users to split, join and rename segments in
    SasMol objects
    """

    def __init__(self, other_self, segname_list, prep_report, log):
        """
        Setup the environment and display to show final segment 
        selection instructions to the user

        """
        self.log = log
        self.json_variables = other_self.json_variables

	self.answer = self.select_segnames(segname_list, prep_report)

        return 
   
    def select_segnames(self, segname_list, prep_report):
  
        timeout = 3600

        prep_dict = {} 
        prep_dict["id"] =  "text_1"
        prep_dict["type"] = "textarea"
        prep_dict["default"] = '\n'.join(prep_report)
        prep_dict["rows"] = len(prep_report) + 6
        prep_dict["cols"] = 180
        prep_dict["fontfamily"] = "monospace"

        label_dict = {}
        label_dict["id"] = "label_1"
        label_dict["type"] = "label"
        label_dict["label"] = "<p><hr><br>Choose Segments<br><br>"

        listbox_dict = {}
        listbox_dict["id"] = "segname_listbox"
        listbox_dict["type"] = "listbox"
        listbox_dict["fontfamily"] = "monospace"

        listbox_dict["values"] = segname_list
        listbox_dict["size"] = 10
        listbox_dict["help"] = "select row(s) and choose an option below: command-click for non-adjacent rows (Mac) or control-click (Windows)"
        listbox_dict["header"] = "choose segment(s) you wish to use in your final model\n click 'submit' when you are finished \n\n"
        listbox_dict["fontsize"] = "0.93em"
        listbox_dict["multiple"] = "true"

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Final Seqment Choice "
        my_question["text"] = "<p>Review system segment report below</p><br><hr><br>"
        my_question["buttons"] = ["submit"]
        my_question["fields"] = [ prep_dict, label_dict, listbox_dict ]

        self.log.info(json.dumps(prep_dict))

        answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        
        return answer

    def create_segment_data(self, other_self):

        my_values = []
        my_returns = []

        i = 0

        for row in other_self.resid_descriptions:
            my_values.append('{0:7s} {1:>6} {2:>10s} {3:>6s} {4:>12s}'.format(
                            row[0], row[2], row[3], row[4], row[5]))

  

        return answer 

class SegnameEditor():
    """
    Interface to allow users to split, join and rename segments in
    SasMol objects
    """

    def __init__(self, other_self, segnames, pdbscan_report, log):
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
        self.json_variables = other_self.json_variables
        self.log = log

        # Get initial locations of segment name changes in description list
        other_self.starting_breaks = numpy.where(
            other_self.resid_descriptions[:-1, 0] != other_self.resid_descriptions[1:, 0])[0]

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

    def create_segment_data(self, other_self):

        my_values = []
        my_returns = []

        i = 0

        for row in other_self.resid_descriptions:
            my_values.append('{0:7s} {1:>6} {2:>10s} {3:>6s} {4:>12s}'.format(
                            row[0], row[2], row[3], row[4], row[5]))
            my_returns.append(i)
            i += 1

        return my_values, my_returns

    def display_and_query_segments_loop(self, other_self, mol, pdbscan_report):

        timeout = 3600

        listbox_dict = {}
        listbox_dict["id"] = "segment_list_box"
        listbox_dict["type"] = "listbox"
        listbox_dict["fontfamily"] = "monospace"

        my_values, my_returns = self.create_segment_data(other_self)

        h_list = ['Segname', 'Resid', 'Resname', 'Chain', 'Moltype']
        header = '{0:<7s} {1:>7s} {2:>10s}   {3:<8s} {4:<8s}'.format(h_list[0], h_list[1], h_list[2], h_list[3], h_list[4])

        listbox_dict["values"] = my_values
#        listbox_dict["returns"] = my_returns
        listbox_dict["size"] = 10 
        listbox_dict["help"] = "select a row and choose an option below"
        listbox_dict["header"] = header                          
        listbox_dict["fontsize"] = "0.93em"

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Seqment Editor"
        my_question["text"] = "<p>Edit Segments Definitions Below</p><br><hr><br>"
        my_question["buttons"] = ["split", "join", "rename", "accept"]
        my_question["fields"] = [listbox_dict]

        self.log.info(json.dumps(listbox_dict))

        self.answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        
        return 
   
    def query_new_segname(self):

        timeout = 3600

        segname_dict = {}
        segname_dict["id"] = "new_segname"
        segname_dict["type"] = "text"
        segname_dict["help"] = "enter new segname : 1 to 4 characters only"
        segname_dict["norow"] = "true"
        segname_dict["label"] = "new segname"
        segname_dict["default"] = "A"

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "PDB Rx Seqment Editor"
        my_question["text"] = "<p>Enter New Segname Below</p><br><hr><br>"
        my_question["buttons"] = ["accept"]
        my_question["fields"] = [segname_dict]

        answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        
        return answer
   
    
