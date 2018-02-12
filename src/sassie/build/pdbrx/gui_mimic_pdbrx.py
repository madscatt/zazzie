import sys, os
import logging

sys.path.append('./')

import sassie.build.pdbrx.pdb_rx as pdb_rx
import sassie.util.sasconfig as sasconfig

import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
pdbfile = 'testing/data/5E3L.pdb'
#pdbfile = '3rki.pdb'
#pdbfile = 'testing/data/3MWY.pdb'
topfile = os.path.join(sasconfig.__bin_path__,'toppar','top_all27_prot_na.inp')
use_defaults = False
gui = 'terminal'
#gui = 'sassie-web'

#### end user input ####
#### end user input ####
#### end user input ####

logging.basicConfig()

svariables['runname'] = (runname,'string')
svariables['pdbfile'] = (pdbfile,'string')
svariables['topfile'] = (topfile,'string')
svariables['defaults'] = (use_defaults,'boolean')
svariables['gui'] = (gui,'string')

error,variables = input_filter.type_check_and_convert(svariables)
if(len(error)>0):
    print 'error = ',error
    sys.exit()
                                 
txtQueue = multiprocessing.JoinableQueue()

scan = pdb_rx.PDBRx()
scan.main(variables,txtQueue)

this_text = txtQueue.get(True, timeout=0.1)
