#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Get user to choose segments to be included in model from curses interface

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

import picker

def select_segnames(segname_list, system_description):
    '''
    Give user a choice of segments to include in the final model

    @type segname_list :  list
    @param segname_list:  Available segment names
    @type  system_description:  list
    @param system_description:  Text lines describing the segments contents
    @return:
    '''

    title_txt = 'Segments to include in model'

    chosen_segnames = picker.Picker(
        title = title_txt,
        options = segname_list,
        info = system_description
        ).getSelected()

    return chosen_segnames