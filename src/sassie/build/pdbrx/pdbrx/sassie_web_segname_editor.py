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

from __future__ import division  # You don't need this in Python3

from math import ceil
import json
import sys
from StringIO import StringIO
import numpy as np

import sassie.util.communication as communication


class SegnameEditor():
    """
    Curses interface to allow users to split, join and rename segments in
    SasMol objects
    """

    def __init__(self, segnames, resid_descriptions, json_variables, pdbscan_report, log):
        """
        Setup the curses environment and display to show residue infromation and
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
        self.starting_breaks = np.where(
            self.resid_descriptions[:-1, 0] != self.resid_descriptions[1:, 0])[0]

	#self.answer = self.ask_question_edit_segmentation(pdbscan_report)
	self.answer = self.ask_question_edit_segmentation_dict(pdbscan_report)

        return 
    
    def ask_question_edit_segmentation_dict(self, pdbscan_report):

        pdbscan_dict = {} 
        pdbscan_dict["id"] =  "text_1"
        pdbscan_dict["type"] = "text"
        pdbscan_dict["label"] = pdbscan_report

        my_question = {}
        my_question["id"] = "q1"
        my_question["title"] = "title"
        my_question["title"] = "title"
        my_question["text"] = "<p>Do you want to edit the system segment definitions?</p><hr>"
        my_question["buttons"] = ["yes", "no"]
        my_question["fields"] = [ pdbscan_dict ]

        #string_my_question = ''.join('{} : {}'.format(key, val) for key, val in my_question.items())
        #my_question = json.loads(string_my_question)

        self.log.info(json.dumps(pdbscan_dict))

        answer = communication.tcpquestion(self.json_variables, my_question);
        
        return answer

    def ask_question_edit_segmentation(self, pdbscan_report):

        my_question = '''
{
    "id" : "q1"
    ,"title" : "Segment Edit Query"
    ,"text" : "<p>Do you want to edit the system segment definitions?</p><hr>"
    ,"buttons"  : [
                    "yes"
                    ,"no"
                   ]
    ,"fields" : [
                   {
                   "type"    : "label"
                   ,"id"      : "segment_edit_list_box"
                   ,"label"    : "select yes or no"
                   ,"fontfamily" :"monospace"
                   } 
    ]
}
'''.strip()

        answer = communication.tcpquestion(self.json_variables, my_question);
        
        return answer

    def valid_segname(self, segname):
        """
        Check that the input segment name is valid to use for a new segment,
        i.e is 4 characters long or less and not an existing segment name.

        @type segname :  str
        @param segname:  Proposed segment name
        @rtype :  bool
        @return:  Is the input segname valid
        """

        valid = False

        if len(segname) <= 4 and segname not in self.segnames:
            valid = True

        return valid

    def create_display_lines(self):
        """
        Format residue information for display

        @return:
        """

        input_data = self.resid_descriptions

        menu_input = []

        for row in input_data:
            menu_input.append('{0:7s} {1:>6} {2:7s} {3:5s} {4:8s}'.format(
                row[0], row[2], row[3], row[4], row[5]))

        return menu_input

    def split_segnames(self, ndx, new_segname):
        """
        Split an existing segment and name the newly created segment.

        @type ndx :  int
        @param ndx:  Index of the residue selected for segment break (in list
                     of residue descritions)
        @type new_segname :  str
        @param new_segname:  Name to be applied to the newly created segment
        @return:
        """

        resid_desc = self.resid_descriptions

        last_ndx = len(resid_desc) - 1

        current_segname = resid_desc[ndx][0]

        if ndx != 0:
            previous_segname = resid_desc[ndx - 1][0]
        else:
            previous_segname = ''

        if previous_segname == current_segname:

            updated_data = []

            for i in range(len(resid_desc)):

                line = self.resid_descriptions[i]

                if i >= ndx and self.resid_descriptions[i][0] == current_segname:
                    line[0] = new_segname

                updated_data.append(line)

            self.resid_descriptions = np.array(updated_data)

            self.segnames.append(new_segname)

        return

    def rename_segment(self, ndx, new_segname):
        """
        Change the name of selected segment (the one including the selected
        residue).

        @type ndx :  int
        @param ndx:  Index of the user selected residue
        @type new_segname :  str
        @param new_segname:  New name for segment
        @return:
        """

        target_segname = self.resid_descriptions[ndx][0]

        updated_data = []

        for line in self.resid_descriptions:
            if line[0] == target_segname:
                line[0] = new_segname
            updated_data.append(line)

        self.resid_descriptions = np.array(updated_data)

        self.segnames = [x if (x != target_segname)
                         else new_segname for x in self.segnames]

        return

    def join_segnames(self, ndx):
        """
        Join segment containing ndx-th residue to the previous segment.

        @type ndx:          integer
        @param ndx:         Index of the residue that starts segment to join
                            previous segment
        @return:            Updated array of lines containing: segnames, indices,
                            resids, resnames, chains, moltypes
        """

        resid_desc = self.resid_descriptions

        last_ndx = len(resid_desc) - 1

        current_segname = resid_desc[ndx][0]
        current_moltype = resid_desc[ndx][-1]

        if ndx != 0:

            previous_segname = resid_desc[ndx - 1][0]
            previous_moltype = resid_desc[ndx - 1][-1]

            moltype_match = (previous_moltype == current_moltype)
            resid_match = (resid_desc[ndx - 1][2] < resid_desc[ndx][2])

        else:
            previous_segname = ''
            # No previous segment, so joining makes no sense
            moltype_match = False
            resid_match = False

        segname_mismatch = (previous_segname != current_segname)

        acceptable_join = moltype_match and resid_match and segname_mismatch

        error = ''

        if acceptable_join:

            updated_data = []

            for i in range(len(resid_desc)):

                line = resid_desc[i]

                if i >= ndx and resid_desc[i][0] == current_segname:
                    line[0] = previous_segname

                updated_data.append(line)

            self.resid_descriptions = np.array(updated_data)

            self.segnames.remove(current_segname)

        else:

            if not segname_mismatch:
                error = 'Segments with the same name cannot be joined'
            elif not resid_match:
                error = 'Joined segment must start with higher resid'
            else:
                error = 'Joined segments must have same moltype'

        return error

    def get_segment_starts(self):
        """
        Get indicies where the resid descriptions change segment name.

        @return:
        """

        resid_desc = self.resid_descriptions

        new_breaks = np.where(resid_desc[:-1, 0] != resid_desc[1:, 0])[0]

        if (new_breaks != self.starting_breaks).any():

            new_breaks += 1
            new_breaks = np.append([0], new_breaks)

            start_segnames = {}

            # Note residue descriptions give last atom in that residue
            # account for that here
            for start_ndx in new_breaks:
                if start_ndx == 0:
                    sasmol_index = 0
                else:
                    sasmol_index = int(resid_desc[start_ndx - 1][1]) + 1

                start_segnames[sasmol_index] = resid_desc[start_ndx][0]

        else:
            start_segnames = {}

        return json.dumps(start_segnames)


def get_input_variables_json(input_json):
    """
    Parse input JSON to produce list segnames, residue descritions and
    expected maximum number of screen rows.

    @type input_json :
    @param input_json:
    @rtype :  list, list, int
    @return:  List of segment names.
              List of tuples describing segname,
              first atomic index, resid, resname, chain
              and moltype for each residue.
              Maximum number of screen lines
    """

    json_stingio = StringIO(input_json)
    json_variables = json.load(json_stingio)

    segnames = json_variables['segnames']

    resid_descriptions = json_variables['resid_descriptions']

    dtypes = np.dtype('a10, int, int, a10, a10, a10')

    for index, info in enumerate(resid_descriptions):
        resid_descriptions[index] = tuple(info)

    resid_descriptions = np.array(resid_descriptions, dtypes)

    if 'max_row' in json_variables:
        max_row = json_variables['max_row']
    else:
        max_row = 10

    return segnames, resid_descriptions, max_row


def main():
    """
    Read JSON definition of segments, residues annd screen size from argv.
    Run segmentation editing and provide updated segmentation information as
    output JSON.

    @return:
    """

    input_json = StringIO(sys.argv[1])
    segnames, resid_descriptions, max_row = get_input_variables_json(
        input_json)

    edited_segments = SegnameEditor(
        segnames, resid_descriptions, max_row).get_segment_starts()

    print json.dumps({'segname_starts': edited_segments})

    return


if __name__ == "__main__":
    # execute only if run as a script
    main()
