
'''
GUI_MIMIC for build_utilities
'''

import sys

import tool_methods
import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
pdbfile = '../../../developer_files_for_testing/torsion_angle_monte_carlo/hiv1_gag.pdb'

seed = '0, 123'  # set this to '1,123' if you want to set the seed or '0,123' if not

#### end user input ####
#### end user input ####
#### end user input ####


svariables['runname'] = (runname, 'string')
svariables['pdbfile'] = (pdbfile, 'string')
svariables['seed'] = (seed, 'int_array')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()

#import pprint; pprint.pprint(variables); exit()

txtQueue = multiprocessing.JoinableQueue()
build_utilities = tool_methods.build_utilities()
build_utilities.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)

#print 'in GUI and txtOutput = ', this_text, '\n'





