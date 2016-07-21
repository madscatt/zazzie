'''
Driver method to run the Altens module
'''

import sys

sys.path.append('./')
import altens

#import sassie.calculate.sascalc as sascalc
import sassie.interface.input_filter as input_filter
#import sassie.interface.sascalc_filter as sascalc_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
rdc_input_file = "RDC.txt"
nh_vector_coordinate_file = "vNH.txt"
active_residue_list_file = "reslist.txt"
number_of_monte_carlo_steps = "500"


#### end user input ####
#### end user input ####
#### end user input ####

svariables={}

svariables['runname'] = (runname,'string')
svariables['rdc_input_file'] = (rdc_input_file,'string')
svariables['nh_vector_coordinate_file'] = (nh_vector_coordinate_file,'string')
svariables['active_residue_list_file'] = (active_residue_list_file,'string')
svariables['number_of_monte_carlo_steps'] = (number_of_monte_carlo_steps,'int')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()
else:
    pass
#    error=sascalc_filter.check_sascalc(variables)
#    if(len(error) != 0):
#        print 'error = ',error
#        sys.exit()


import time; start = time.time()
txtQueue = multiprocessing.JoinableQueue()
altens = altens.altens()
altens.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)
print "time used: ",time.time()-start

#print 'in GUI and txtOutput = ', this_text, '\n'


