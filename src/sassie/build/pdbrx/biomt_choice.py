#!/usr/bin/env python

import picker

def select_biomt(available_biomt, biomt_description):

    biomt_list = [str(x) for x in available_biomt]

    title_txt = 'Pick BIOMTs to apply to model'

    chosen_biomt = picker.Picker(
        title = title_txt,
        options = biomt_list,
        info = biomt_description
        ).getSelected()

    return chosen_biomt