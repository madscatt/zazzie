#!/usr/bin/env python

import picker

def select_altloc(segname, resid, altloc_list):

    title_txt = 'AltLocs for seg {0:s} resid {1:d}'.format(segname, resid)

    info = ['Choose one of the altloc labels detected in the structure']

    chosen_altloc = picker.Picker(
        title = title_txt,
        options = altloc_list,
        info = info,
        mutually_exclusive = True
        ).getSelected()

    return chosen_altloc