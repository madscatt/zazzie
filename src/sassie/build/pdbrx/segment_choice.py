#!/usr/bin/env python

import picker

def select_segnames(segname_list, system_description):

    title_txt = 'Segments to include in model'

    chosen_segnames = picker.Picker(
        title = title_txt,
        options = segname_list,
        info = system_description
        ).getSelected()

    return chosen_segnames