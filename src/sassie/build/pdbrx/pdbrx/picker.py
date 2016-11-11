#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Curses Picker

Based on the code from:
https://github.com/pp19dd/picker

The original licensing text read:
This code is available for use under CC0 (Creative Commons 0 - universal).
You can copy, modify, distribute and perform the work, even for commercial
purposes, all without asking permission. For more information, see LICENSE.md or
https://creativecommons.org/publicdomain/zero/1.0/

usage:
opts = Picker(
   title = 'Delete all files',
   options = ["Yes", "No"]
).get_selected()


returns a simple list
cancel returns False

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

import curses
import curses.wrapper


class Picker:
    """Allows you to select from a list with curses"""
    stdscr = None
    win_pick = None
    title = ""
    arrow = ""
    footer = ""
    more = ""
    c_selected = ""
    c_empty = ""

    cursor = 0
    offset = 0
    selected = 0
    selcount = 0
    aborted = False

    window_height = 5
    window_width = 79
    window_y = 2
    window_x = 4

    all_options = []
    length = 0

    win_info = None
    info_height = 8
    info_width = 79
    info_y = 2
    info_x = 4
    info_cursor = 0
    info_offset = 0

    def curses_start(self):
        """
        Setup curses environment
        @return:
        """

        self.stdscr = curses.initscr()
        curses.noecho()
        curses.cbreak()

        nlines = curses.LINES
        ncols = curses.COLS

        self.width = curses.COLS - 1
        self.info_width = self.width

        n_win_lines = (nlines - 8) / 2

        self.window_height = n_win_lines
        self.info_height = n_win_lines

        self.window_y = 3 + self.info_height

        self.win_pick = curses.newwin(
            5 + self.window_height,
            self.window_width,
            self.window_y,
            self.window_x
        )

        if self.info:
            self.win_info = curses.newwin(
                3 + self.info_height,
                self.info_width,
                self.info_y,
                self.info_x
            )

    def curses_stop(self):
        """
        Clean up curses environemnt after exit
        @return:
        """

        curses.nocbreak()
        self.stdscr.keypad(0)
        self.stdscr.erase()
        curses.echo()
        curses.endwin()

    def get_selected(self):
        """
        Return input after user makes choice

        @return:
        """

        if self.aborted:
            return (False)

        ret_s = filter(lambda x: x["selected"], self.all_options)
        ret = map(lambda x: x["label"], ret_s)
        
        return (ret)

    def redraw_info(self):
        """
        Redraw information section after user has moved viewing window.

        @return:
        """

        info_offset = self.info_offset

        line_no = 0

        visible_lines = self.info[info_offset:info_offset + self.info_height]

        for line in visible_lines:
            self.win_info.addstr(line_no, 0, line)
            line_no += 1

        return

    def redraw(self):
        """
        Redraw screen after update.

        @return:
        """

        self.win_pick.clear()
        self.win_pick.border(
            self.border[0], self.border[1],
            self.border[2], self.border[3],
            self.border[4], self.border[5],
            self.border[6], self.border[7]
        )

        self.win_pick.addstr(
            self.window_height + 4, 5, " " + self.footer + " "
        )

        position = 0

        available = self.all_options[self.offset:self.offset + self.window_height + 1]

        for option in available:

            if option["selected"]:

                line_label = self.c_selected + " "

            else:

                line_label = self.c_empty + " "

            self.win_pick.addstr(position + 2, 5, line_label + option["label"])
            position = position + 1

        # hint for more content above
        if self.offset > 0:
            self.win_pick.addstr(1, 5, self.more)

        # hint for more content below
        if self.offset + self.window_height <= self.length - 2:
            self.win_pick.addstr(self.window_height + 3, 5, self.more)

        self.win_pick.addstr(0, 5, " " + self.title + " ")
        self.win_pick.addstr(
            0, self.window_width - 8,
            " " + str(self.selcount) + "/" + str(self.length) + " "
        )
        self.win_pick.addstr(self.cursor + 2, 1, self.arrow)

        if self.info:
            self.win_info.clear()

            self.win_info.addstr(
                self.info_height + 2, 15, " Scroll information with Page Up/Page Down"
            )

            self.redraw_info()
            self.win_info.refresh()

        self.win_pick.refresh()

    def check_info_cursor(self):
        """
        Ensure cursor is within the list of options.

        @return:
        """

        if self.info_cursor < 0:

            if self.info_offset > 0:
                self.info_offset -= 1

        elif self.info_cursor > 0:
            self.info_offset += 1

            if self.info_offset >= len(self.info):
                self.info_offset -= 1

        self.info_cursor = 0

    def check_cursor_up(self):
        """
        Update cursor position after user scrolls up

        @return:
        """

        if self.cursor < 0:
            self.cursor = 0
            if self.offset > 0:
                self.offset = self.offset - 1

    def check_cursor_down(self):
        """
        Update cursor position after user scrolls down
        @return:
        """

        if self.cursor >= self.length:
            self.cursor = self.cursor - 1

        if self.cursor > self.window_height:
            self.cursor = self.window_height
            self.offset = self.offset + 1

            if self.offset + self.cursor >= self.length:
                self.offset = self.offset - 1

    def select_mutually_exclusive(self):
        """
        Select option if those in list are mutually exclusive - i.e. selection
        of one item inverts the selection of the others.

        @return:
        """

        old_value = self.all_options[self.selected]["selected"]
        new_value = not old_value

        for i in range(len(self.all_options)):
            if i == self.selected:
                self.all_options[i]["selected"] = new_value
            else:
                self.all_options[i]["selected"] = old_value

    def curses_loop(self, stdscr):
        """
        Main loop

        @param stdscr:
        @return:
        """
        while 1:
            self.redraw()
            c = stdscr.getch()

            if c == ord('q') or c == ord('Q'):
                self.aborted = True
                break

            elif c == curses.KEY_UP:
                self.cursor = self.cursor - 1

            elif c == curses.KEY_DOWN:
                self.cursor = self.cursor + 1

            elif c == ord(' '):

                if self.mutually_exclusive:
                    self.select_mutually_exclusive()
                else:
                    self.all_options[self.selected]["selected"] = \
                        not self.all_options[self.selected]["selected"]

            elif c == ord('i') or c == ord('I'):

                for option in self.all_options:
                    option['selected'] = not option['selected']

            elif c == curses.KEY_NPAGE:
                self.info_cursor += 1

            elif c == curses.KEY_PPAGE:
                self.info_cursor -= 1

            elif c == 10:
                break

            # deal with interaction limits
            self.check_cursor_up()
            self.check_cursor_down()

            self.check_info_cursor()

            # compute selected position only after dealing with limits
            self.selected = self.cursor + self.offset

            temp = self.get_selected()
            self.selcount = len(temp)

    def __init__(
            self,
            options,
            title='Select',
            arrow="-->",
            footer="Space = toggle, Enter = accept, i = invert, q = cancel",
            more="...",
            border="||--++++",
            c_selected="[X]",
            c_empty="[ ]",
            info=[],
            mutually_exclusive=False
    ):

        """ Setup screen for display"""

        self.title = title
        self.arrow = arrow
        self.footer = footer
        self.more = more
        self.border = border
        self.c_selected = c_selected
        self.c_empty = c_empty

        self.mutually_exclusive = mutually_exclusive

        self.info = info

        self.all_options = []

        for option in options:
            self.all_options.append({
                "label": option,
                "selected": False
            })
            self.length = len(self.all_options)

        if self.mutually_exclusive:
            self.all_options[0]["selected"] = True

        if self.length < self.window_height - 1:
            self.window_height = self.length + 1

        if info:
            if len(info) < 15:
                self.info_height = len(info) + 2
            else:
                self.info_height = 17

            self.window_y = self.info_y + self.info_height + 5

        self.curses_start()
        curses.wrapper(self.curses_loop)
        self.curses_stop()
