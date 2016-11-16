#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Present user with a choice of which BIOMT records to apply to the input
structure.

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

from . import picker


def select_biomt(available_biomt, biomt_description):
    '''
    Provide a ncurses interface to pick which BIOMT transforms to apply to the
    structure.

    @type available_biomt :  list
    @param available_biomt:  Identifiers of the BIOMT records
    @type biomt_description :  list
    @param biomt_description:  Contains lines of description for each BIOMT
    @rtype :  list
    @return:  Identifiers of selected BIOMT records
    '''

    biomt_list = [str(x) for x in available_biomt]

    title_txt = 'Pick BIOMTs to apply to model'

    chosen_biomt = picker.Picker(
        title=title_txt,
        options=biomt_list,
        info=biomt_description
        ).get_selected()

    return chosen_biomt
