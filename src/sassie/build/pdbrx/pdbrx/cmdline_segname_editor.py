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

from . import segname_utils as segname_utils

import curses
import curses.wrapper

from math import ceil
import json
import sys
from StringIO import StringIO
import numpy as np


class SegnameEditor():
    """
    Curses interface to allow users to split, join and rename segments in
    SasMol objects
    """

    def __init__(self, other_self, mol, max_row=10):
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

        # Get initial locations of segment name changes in description list
        other_self.starting_breaks = np.where(
            other_self.resid_descriptions[:-1, 0] != other_self.resid_descriptions[1:, 0])[0]

        self.display_lines = self.create_display_lines(other_self)

        self.screen = None

        self.curses_start()

        if max_row + 8 >= curses.LINES:
            self.max_row = 5
        else:
            self.max_row = max_row

        self.row_num = len(self.display_lines)

        # Residue descriptions shown in a scrollable window
        # Pages define list of simultaneously viewable lines
        self.pages = int(ceil(self.row_num / self.max_row))

        # Initial viewing and selection location
        self.position = 1
        self.page = 1

        self.screen_setup()
        curses.wrapper(self.curses_loop, other_self, mol)
        self.curses_stop()

        return

    def curses_start(self):
        """
        Initial setup of curses environment.

        @return:
        """

        self.screen = curses.initscr()
        curses.noecho()
        curses.cbreak()
        curses.start_color()

        self.screen.keypad(1)

        return

    def screen_setup(self):
        """
        Define the viewable areas in curses and the constant content (such as
        user instructions).

        @return:
        """

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
        """
        Stop curses environment and tidy up.

        @return:
        """

        self.screen.clear()

        curses.nocbreak()
        self.screen.keypad(0)
        curses.echo()
        curses.endwin()

        return

    def redraw(self):
        """
        Redraw screen after change in state.

        @return:
        """

        max_row = self.max_row
        page = self.page
        row_num = self.row_num
        display_lines = self.display_lines
        position = self.position

        self.box.erase()
        self.screen.border(0)
        self.box.border(0)

        for i in range(1 + (max_row * (page - 1)),
                       max_row + 1 + (max_row * (page - 1))):

            if row_num == 0:

                self.box.addstr(1, 1, "There aren't strings",
                                self.highlightText)
            else:

                paged_lines = max_row * (page - 1)

                if i + paged_lines == position + paged_lines:

                    self.box.addstr(i - paged_lines, 2,
                                    display_lines[i - 1], self.highlightText)
                else:

                    self.box.addstr(i - paged_lines, 2,
                                    display_lines[i - 1], self.normalText)

                if i == row_num:
                    break

        self.last_occupied = i

        self.screen.refresh()
        self.box.refresh()

        return

    #def curses_loop(stdscr, self, other_self, mol):
    def curses_loop(self, stdscr, other_self, mol):
        """
        Main loop to obtain and interpret user input

        @param stdscr:
        @return:
        """

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

                if segname_utils.valid_segname(new_segname, mol.segnames()):
                    segname_utils.split_segnames(other_self, mol, ndx, new_segname)
                    self.display_lines = self.create_display_lines(other_self)
                    self.box.refresh()

            elif x in [ord('j'), ord('J')]:

                ndx = self.position - 1

                error = segname_utils.join_segnames(other_self, mol, ndx)

                if error:
                    self.screen.addstr(max_row + 4, 3, error)

                self.display_lines = self.create_display_lines(other_self)
                self.box.refresh()

            elif x in [ord('r'), ord('R')]:

                ndx = self.position - 1
                current_segment = other_self.resid_descriptions[ndx][0]

                curses.echo()
                self.screen.addstr(
                    max_row + 4, 3, "Rename " + current_segment + " to: ")
                new_segname = self.screen.getstr(max_row + 4, 19)
                curses.noecho()

                self.screen.addstr(
                    max_row + 4, 3, "                                                         ")
                self.screen.refresh()

                if segname_utils.valid_segname(new_segname, mol.segnames()):
                    segname_utils.rename_segment(other_self, mol, ndx, new_segname)
                    self.display_lines = self.create_display_lines(other_self)
                    self.box.refresh()

            elif x in [ord('a'), ord('A')]:
                break

        return

    def create_display_lines(self, other_self):
        """
        Format residue information for display

        @return:
        """

        input_data = other_self.resid_descriptions

        menu_input = []

        for row in input_data:
            menu_input.append('{0:7s} {1:>6} {2:7s} {3:5s} {4:8s}'.format(
                row[0], row[2], row[3], row[4], row[5]))

        return menu_input


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

    #edited_segments = SegnameEditor(
    #    segnames, resid_descriptions, max_row).get_segment_starts()
    
    edited_segments = SegnameEditor(
        other_self).get_segment_starts()

    print json.dumps({'segname_starts': edited_segments})

    return


if __name__ == "__main__":
    # execute only if run as a script
    main()
