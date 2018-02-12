'''
Driver method to run the Altens module
'''

import sys

sys.path.append('./')
import altens

import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
rdc_input_file = "RDC_D_K63.txt"
pdbfile='1D3Z_mod1_new.pdb'
residue_list_file = "reslist.txt"
mcon=True
number_of_monte_carlo_steps = '500' 
seed = '1,123'

#### end user input ####
#### end user input ####
#### end user input ####

svariables={}

svariables['runname'] = (runname,'string')
svariables['rdc_input_file'] = (rdc_input_file,'string')
svariables['pdbfile'] = (pdbfile,'string')
svariables['residue_list_file'] = (residue_list_file,'string')
svariables['mcon'] =(mcon,'boolean')
svariables['number_of_monte_carlo_steps'] = (number_of_monte_carlo_steps,'int')
svariables['seed'] = (seed, 'int_array')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()
else:
    pass

import time; start = time.time()
txtQueue = multiprocessing.JoinableQueue()
altens = altens.altens()
altens.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)
print ("time used: ",time.time()-start)

#print 'in GUI and txtOutput = ', this_text, '\n'


