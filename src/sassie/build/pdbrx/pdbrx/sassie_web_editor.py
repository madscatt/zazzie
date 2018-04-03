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

class FastaEditor():
    """
    Interface to for input/ouput of questions related to sequence
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
        my_question["text"] = "<p>Review system seguence report below</p><br><hr><br>"
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

	self.answer = self.ask_question_edit_segmentation(pdbscan_report)

        return 
    
    def ask_question_edit_segmentation(self, pdbscan_report):

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

        answer = communication.tcpquestion(self.json_variables, my_question, timeout);
        
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
   
    
