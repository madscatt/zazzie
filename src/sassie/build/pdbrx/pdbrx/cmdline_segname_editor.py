#!/usr/bin/env python
from __future__ import division  # You don't need this in Python3

import curses
import curses.wrapper

from math import *
import numpy as np
import json
import sys
from StringIO import StringIO


class SegnameEditor():

    def __init__(self, segnames, resid_descriptions, max_row=10):

        self.segnames = segnames
        self.resid_descriptions = resid_descriptions

        self.starting_breaks = np.where(
            self.resid_descriptions[:-1, 0] != self.resid_descriptions[1:, 0])[0]

        self.display_lines = self.create_display_lines()

        self.screen = None

        self.curses_start()

        if max_row + 8 >= curses.LINES:
            self.max_row = 5
        else:
            self.max_row = max_row

        self.row_num = len(self.display_lines)

        self.pages = int(ceil(self.row_num / self.max_row))
        self.position = 1
        self.page = 1

        self.screen_setup()
        curses.wrapper(self.curses_loop)
        self.curses_stop()

        return

    def curses_start(self):

        self.screen = curses.initscr()
        curses.noecho()
        curses.cbreak()
        curses.start_color()

        self.screen.keypad(1)

        return

    def screen_setup(self):

        curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
        self.highlightText = curses.color_pair(1)
        self.normalText = curses.A_NORMAL

        self.screen.border(0)
        curses.curs_set(0)

        max_row = self.max_row

        column_heads = '{0:7s} {1:>6s} {2:7s} {3:5s} {4:8s}'.format(
            'Segname', 'Resid', 'Resname', 'Chain', 'Moltype')

        self.screen.addstr(1, 3, column_heads)

        if max_row + 12 < curses.LINES - 1:

            self.screen.addstr(max_row + 6, 5, "Scroll using Up/Down arrows")
            self.screen.addstr(max_row + 8, 5, "Available commands:")
            self.screen.addstr(max_row + 9, 5, "s: Split")
            self.screen.addstr(max_row + 10, 5, "j: Join")
            self.screen.addstr(max_row + 11, 5, "r: Rename")
            self.screen.addstr(
                max_row + 12, 5, "a: Accept current segmentation")

        else:
            self.screen.addstr(max_row + 6, 5, "Scroll using Up/Down arrows")
            self.screen.addstr(
                max_row + 7, 5, "Commands: (s)plit, (j)oin, (r)ename")
            self.screen.addstr(
                max_row + 8, 5, "          (a)ccept segmentation")

        self.box = curses.newwin(max_row + 2, 64, 2, 1)
        self.box.box()

    def curses_stop(self):

        self.screen.clear()

        curses.nocbreak()
        self.screen.keypad(0)
        curses.echo()
        curses.endwin()

        return

    def redraw(self):

        max_row = self.max_row
        page = self.page
        row_num = self.row_num
        display_lines = self.display_lines
        position = self.position

        self.box.erase()
        self.screen.border(0)
        self.box.border(0)

        for i in range(1 + (max_row * (page - 1)), max_row + 1 + (max_row * (page - 1))):
            if row_num == 0:
                self.box.addstr(1, 1, "There aren't strings",
                                self.highlightText)
            else:
                if (i + (max_row * (page - 1)) == position + (max_row * (page - 1))):
                    self.box.addstr(i - (max_row * (page - 1)),
                                    2, display_lines[i - 1], self.highlightText)
                else:
                    self.box.addstr(i - (max_row * (page - 1)),
                                    2, display_lines[i - 1], self.normalText)
                if i == row_num:
                    break

        self.last_occupied = i

        self.screen.refresh()
        self.box.refresh()

        return

    def curses_loop(self, stdscr):

        max_row = self.max_row
        pages = self.pages

        while 1:
            self.redraw()

            x = stdscr.getch()

            if x == curses.KEY_DOWN:
                if self.page == 1:
                    if self.position < self.last_occupied:
                        self.position = self.position + 1
                    else:
                        if pages > 1:
                            self.page += + 1
                            self.position = 1 + (max_row * (self.page - 1))
                elif self.page == pages:
                    if self.position < self.row_num:
                        self.position = self.position + 1
                else:
                    if self.position < max_row + (max_row * (self.page - 1)):
                        self.position = self.position + 1
                    else:
                        self.page += 1
                        self.position = 1 + (max_row * (self.page - 1))

            elif x == curses.KEY_UP:
                if self.page == 1:
                    if self.position > 1:
                        self.position = self.position - 1
                else:
                    if self.position > (1 + (max_row * (self.page - 1))):
                        self.position = self.position - 1
                    else:
                        self.page -= 1
                        self.position = max_row + (max_row * (self.page - 1))

            elif x == curses.KEY_LEFT:
                if self.page > 1:
                    self.page -= 1
                    self.position = 1 + (max_row * (self.page - 1))

            elif x == curses.KEY_RIGHT:
                if self.page < pages:
                    self.page += 1
                    self.position = (1 + (max_row * (self.page - 1)))

            elif x in [ord('s'), ord('S')]:

                ndx = self.position - 1

                curses.echo()
                self.screen.addstr(max_row + 4, 3, "Name new segment: ")
                new_segname = self.screen.getstr(max_row + 4, 21)
                curses.noecho()

                self.screen.addstr(
                    max_row + 4, 3, "                                                         ")
                self.screen.refresh()

                if self.valid_segname(new_segname):

                    self.split_segnames(ndx, new_segname)
                    self.display_lines = self.create_display_lines()
                    self.box.refresh()

            elif x in [ord('j'), ord('J')]:

                ndx = self.position - 1

                error = self.join_segnames(ndx)

                if error:
                    self.screen.addstr(max_row + 4, 3, error)

                self.display_lines = self.create_display_lines()
                self.box.refresh()

            elif x in [ord('r'), ord('R')]:

                ndx = self.position - 1
                current_segment = self.resid_descriptions[ndx][0]

                curses.echo()
                self.screen.addstr(
                    max_row + 4, 3, "Rename " + current_segment + " to: ")
                new_segname = self.screen.getstr(max_row + 4, 19)
                curses.noecho()

                self.screen.addstr(
                    max_row + 4, 3, "                                                         ")
                self.screen.refresh()

                if self.valid_segname(new_segname):

                    self.rename_segment(ndx, new_segname)
                    self.display_lines = self.create_display_lines()
                    self.box.refresh()

            elif x in [ord('a'), ord('A')]:
                break

        return

    def valid_segname(self, segname):

        valid = False

        if len(segname) <= 4 and segname not in self.segnames:
            valid = True

        return valid

    def create_display_lines(self):

        input_data = self.resid_descriptions

        menu_input = []

        for row in input_data:

            menu_input.append('{0:7s} {1:>6} {2:7s} {3:5s} {4:8s}'.format(
                row[0], row[2], row[3], row[4], row[5]))

        return menu_input

    def split_segnames(self, ndx, new_segname):

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

        @type masked_data : np.array
        @param masked_data: Array of lines containing: segnames, indices, resids,
                            resnames, chains, moltypes
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

    json_stingio = StringIO(input_json)
    json_variables = json.load(json_stingio)

    segnames = json_variables['segnames']

    resid_descriptions = json_variables['resid_descriptions']

    dt = np.dtype('a10, int,int,a10,a10,a10')

    for ind, l in enumerate(resid_descriptions):
        resid_descriptions[ind] = tuple(l)

    resid_descriptions = np.array(resid_descriptions, dt)

    if 'max_row' in json_variables:
        max_row = json_variables['max_row']
    else:
        max_row = 10

    return segnames, resid_descriptions, max_row


def main():

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
