#!/usr/bin/env python
"""
AltLoc Choice
Provides a commandline UI for selecting an AltLoc to use for a particular
residue.
"""

from . import picker


def select_altloc(segname, resid, altloc_list):
    '''
    Prompt user to select an altloc to use for aresidue for which multiple
    conformations are present in the PDB

    @type segname :  string
    @param segname:  Name of the segment for which the choice is to be applied
    @type resid :    int
    @param resid:    Residue number for which altloc is being selected
    @type altloc_list :  List
    @param altloc_list:  List of altloc identifiers present in the residue
    @rtype :  string
    @return:  Identifier for the selected altloc
    '''

    title_txt = 'AltLocs for seg {0:s} resid {1:d}'.format(segname, resid)

    info = ['Choose one of the altloc labels detected in the structure']

    chosen_altloc = picker.Picker(
        title=title_txt,
        options=altloc_list,
        info=info,
        mutually_exclusive=True
    ).getSelected()

    return chosen_altloc
