# -*- coding: utf-8 -*-
"""
Preprocessor to finalize system description after PDB Scan, this is the first
step in PDB Rx

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

from __future__ import print_function

from . import segname_utils as segname_utils
from . import terminal_utils as terminal_utils
from . import sassie_web_utils as sassie_web_utils

def user_input(other_self, mol, pdbscan_report):
    """
    Get user input from terminal or other source. ONLY TERMINAL CURRENTLY

    @return:
    """
    
    other_self.resid_descriptions = segname_utils.create_residue_descriptions(mol)
    
    if other_self.mvars.user_interface == 'terminal':
            
        terminal_utils.handle_terminal_user_input(other_self, mol)

    elif other_self.mvars.user_interface == 'sassie_web':

        sassie_web_utils.handle_sassie_web_user_input(other_self, mol, pdbscan_report)

    return
