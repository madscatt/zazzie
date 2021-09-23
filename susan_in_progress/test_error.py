# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import locale

error = []
def test_float(variable_value):

    print('variable value before try: ', variable_value)

    try:
        dum = locale.atof(variable_value)
        print('variable value: ', dum)
    except:
        error.append('variable is not a floating number')
    

#            raise
#            error.append('%s is not a floating number', %(variable_string))

variable_value = u'1.7a'

test_float(variable_value)

print('error = ', error)