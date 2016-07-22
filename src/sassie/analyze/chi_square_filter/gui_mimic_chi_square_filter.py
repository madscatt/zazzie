from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

'''
Driver method to run the chi_square filter module
'''

import os
import sys
import string

global_install = True

if global_install:

    import sassie.analyze.chi_square_filter.chi_square_filter as chi_square_filter
    import sassie.interface.input_filter as input_filter
    import sassie.interface.chi_square_filter.chi_square_filter_filter as chi_square_filter_filter

else:

    import chi_square_filter as chi_square_filter

    sys.path.append(os.path.join('..', '..', 'interface'))
    import input_filter as input_filter
    sys.path.append(os.path.join('..', '..', 'interface', 'chi_square_filter'))
    import chi_square_filter_filter as chi_square_filter_filter

import multiprocessing


#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
saspath = os.path.join('..', '..', '..', 'developer_files_for_testing',
                       'chi_square_filter', 'run_1', 'sascalc', 'neutron_D2Op_100')
sasintfile = os.path.join(
    '..', '..', '..', 'developer_files_for_testing', 'chi_square_filter', 'sans_data.dat')
io = '0.19'
number_of_weight_files = '0'
basis_string = ''
weight_file_names = ''
x2highcut = '10.0'
x2highweightfile = 'x2highweights.txt'
x2lowcut = '1.0'
x2lowweightfile = 'x2lowweights.txt'
rghighcut = '60.0'
rghighweightfile = 'rghighweights.txt'
rglowcut = '40.0'
rglowweightfile = 'rglowweights.txt'
sastype = '0'
reduced_x2 = '1'
plotflag = '0'
folder_flag = False

path = '.'
data_path = path

#### end user input ####
#### end user input ####
#### end user input ####


svariables = {}
svariables['runname'] = (str(runname), 'string')
svariables['saspath'] = (str(saspath), 'string')
svariables['sasintfile'] = (str(string.strip(sasintfile, " ")), 'string')
svariables['io'] = (str(io), 'float')
svariables['number_of_weight_files'] = (str(number_of_weight_files), 'int')
svariables['basis_string'] = (str(basis_string), 'string')
svariables['weight_file_names'] = (str(weight_file_names), 'string')
svariables['x2highcut'] = (str(x2highcut), 'float')
svariables['x2highweightfile'] = (str(x2highweightfile), 'string')
svariables['x2lowcut'] = (str(x2lowcut), 'float')
svariables['x2lowweightfile'] = (str(x2lowweightfile), 'string')
svariables['rghighcut'] = (str(rghighcut), 'float')
svariables['rghighweightfile'] = (str(rghighweightfile), 'string')
svariables['rglowcut'] = (str(rglowcut), 'float')
svariables['rglowweightfile'] = (str(rglowweightfile), 'string')
svariables['sastype'] = (str(sastype), 'int')
svariables['reduced_x2'] = (str(reduced_x2), 'int')
svariables['plotflag'] = (str(plotflag), 'int')

svariables['path'] = (path, 'string')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print('error = ' + error)
    sys.exit()

error = chi_square_filter_filter.check_chi_square_filter(
    variables, no_file_check="true")
if len(error) > 0:
    print('error = ' + error)
    sys.exit()

txtQueue = multiprocessing.JoinableQueue()
process = multiprocessing.Process(
    target=chi_square_filter.find_best, args=(variables, txtQueue))
process.start()
