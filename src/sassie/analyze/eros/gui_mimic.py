
FAST = True
FAST = False

'''
GUI_MIMIC for eros
'''

import sys
import os

DEBUG = True

if not DEBUG:
    import sassie.analyze.eros.eros as eros
else:
    sys.path.append('./')
    import eros as eros

import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####


runname = 'run_0'
iq_data_path = os.path.join('rsv_data', 'post_fusion', 'run_0', 'sascalc', 'xray')
iq_data_path = os.path.join('rsv_data', 'cluster_500_prefusion', 'run_0', 'sascalc', 'xray')
iq_data_path = os.path.join('rsv_data', 'post_fusion_all', 'spectra')
iq_data_path = os.path.join('second_rsv_data', 'prefusion', 'com_dist_62.0_spectra')
iq_data_path = os.path.join('second_rsv_data', 'pre_plus_post_62')
#second_rsv_data]$ ls prefusion/
#com_dist_56.0_spectra/ com_dist_62.0_spectra/ com_dist_70.0_spectra/

goal_iq_data_file = os.path.join('rsv_data', 'rsvf_0d_saxs_int.dat')

io = '1.0'

number_of_files_to_use = '1469'
number_of_monte_carlo_steps = '1000'
weight_step_size_fraction = '0.1'
theta = '0.2'
theta = '0.0'
beta = '10000.'
beta = '100.'

reduced_x2 = '1' # 1 == reduced_x2, #0 == chi-square distribution, #2 Pearson's X2, #3 R-factor

local_bokeh_server = True



seed = '0, 123'  # set this to '1,123' if you want to set the seed or '0,123' if not

#### end user input ####
#### end user input ####
#### end user input ####

if FAST:
    number_of_files_to_use = '100'
    number_of_monte_carlo_steps = '100'



svariables['runname'] = (runname, 'string')
svariables['iq_data_path'] = (iq_data_path, 'string')
svariables['goal_iq_data_file'] = (goal_iq_data_file, 'string')
svariables['io'] = (io, 'float')
svariables['local_bokeh_server'] = (local_bokeh_server, 'boolean')
svariables['seed'] = (seed, 'int_array')

svariables['number_of_files_to_use'] = (number_of_files_to_use, 'int')
svariables['number_of_monte_carlo_steps'] = (number_of_monte_carlo_steps, 'int')
svariables['weight_step_size_fraction'] = (weight_step_size_fraction, 'float')
svariables['theta'] = (theta, 'float')
svariables['beta'] = (beta, 'float')
svariables['reduced_x2'] = (reduced_x2, 'int')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()

#import pprint; pprint.pprint(variables); exit()

txtQueue = multiprocessing.JoinableQueue()
my_eros = eros.eros()
my_eros.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)

print 'in GUI and txtOutput = ', this_text, '\n'
