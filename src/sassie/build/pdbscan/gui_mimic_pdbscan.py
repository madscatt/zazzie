import sys
import logging

sys.path.append('./')

import sassie.build.pdbscan.pdb_scan as pdb_scan

import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname='run_0'
pdbfile='testing/data/3MWY.pdb'
pdbfile='testing/data/5712.pdb'
pdbfile='testing/data/286.pdb'
pdbfile='testing/data/609.pdb'
pdbfile='testing/data/324.pdb'
pdbfile='testing/data/0.pdb'
pdbfile='testing/data/1984.pdb'
pdbfile='testing/data/7361.pdb'
pdbfile='testing/data/732.pdb'
pdbfile='testing/data/8068.pdb'
pdbfile='testing/data/hiv1_gag.pdb'
pdbfile='testing/data/small.pdb'
#pdbfile='testing/data/5E3L.pdb'

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

scan = pdb_scan.PDBScan()
scan.main(variables,txtQueue)

this_text = txtQueue.get(True, timeout=0.1)
