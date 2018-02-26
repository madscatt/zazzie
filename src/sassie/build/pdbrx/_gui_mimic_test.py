import sys
import logging

sys.path.append('./')

import test

import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname='run_0'
#pdbfile='testing/data/3OVO.pdb'
pdbfile='testing/data/header_coord_mismatch.pdb'

#### end user input ####
#### end user input ####
#### end user input ####

logging.basicConfig()

svariables['runname'] = (runname,'string')
svariables['pdbfile'] = (pdbfile,'string')

error,variables = input_filter.type_check_and_convert(svariables)
if(len(error)>0):
    print 'error = ',error
    sys.exit()
                                 
txtQueue=multiprocessing.JoinableQueue()   

scan = test.PDBScan()
scan.main(variables,txtQueue)

this_text = txtQueue.get(True, timeout=0.1)
