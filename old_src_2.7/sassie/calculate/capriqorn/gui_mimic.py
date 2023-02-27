'''
Driver method to run the Capriqorn module
'''

import sys, os

import sassie.calculate.capriqorn.capriqorn as capriqorn
import sassie.interface.input_filter as input_filter
import sassie.calculate.capriqorn.capriqorn_utils as capriqorn_utils
import sassie.interface.capriqorn_filter as capriqorn_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
pdbfile = 'protein.pdb'
#dcdfile = 'protein.crdbox'
dcdfile = 'protein_10frames_from_crdbox.dcd'
xstfile = 'protein_10frames_from_crdbox.xst'
aliasfile= 'alias.dat' # the atomic type dictionary that is needed by capriqorn

number_q_values='21'
q_max='0.2'

create_alias_flag = False

core_path = os.getcwd()

#### end user input ####
#### end user input ####
#### end user input ####


svariables={}

svariables['runname'] = (runname,'string')
svariables['pdbfile'] = (pdbfile,'string')
svariables['dcdfile'] = (dcdfile,'string')
if dcdfile[-3:]=='dcd':
    svariables['xstfile'] = (xstfile,'string')
svariables['create_alias_flag'] = (create_alias_flag,'boolean')
svariables['aliasfile'] = (aliasfile,'string')
svariables['number_q_values'] = (number_q_values,'int')
svariables['q_max'] = (q_max,'float')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()

else:
    error=capriqorn_filter.check_capriqorn(variables)
    if(len(error) != 0):
        print 'error = ',error
        sys.exit()

#import pprint; pprint.pprint(variables); exit()

import time; start = time.time()
txtQueue = multiprocessing.JoinableQueue()
capriqorn = capriqorn.capriqorn()
capriqorn.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)
print "time used: ",time.time()-start

#print 'in GUI and txtOutput = ', this_text, '\n'


