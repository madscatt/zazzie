'''
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

'''
import os
import string
import locale
import numpy
import io

import sasmol.system as system


def read_structure_files(other_self):

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    log.debug('in read_structure_files')
    pgui('in read_structure_files')
    avars = other_self.asaxs_variables

    mol = system.Molecule()

    return



