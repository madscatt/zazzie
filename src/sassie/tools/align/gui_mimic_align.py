'''
Driver method to run the align module
'''

import sys

#import sassie.tools.align as align
import align as align
import sassie.interface.input_filter as input_filter
import sassie.interface.align_filter as align_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
path = './'
#infile = 'hiv1_gag_200_frames.dcd'
infile = 'hiv1_gag.pdb'
pdbmol1 = 'hiv1_gag.pdb'
pdbmol2 = 'hiv1_gag.pdb'
ofile = 'aligned_hiv1_gag.pdb'
#ofile = 'aligned_hiv1_gag.dcd'
basis1 = 'CA'
basis2 = 'CA'
lowres1 = '145'
lowres2 = '145'
highres1 = '350'
highres2 = '350'
ebasis1 = 'None'
ebasis2 = 'None'


#### end user input ####
#### end user input ####
#### end user input ####


svariables={}

svariables['runname'] = (runname,'string')
svariables['path'] = (path,'string')
svariables['infile'] = (infile,'string')
svariables['pdbmol1'] = (pdbmol1,'string')
svariables['pdbmol2'] = (pdbmol2,'string')
svariables['ofile'] = (ofile,'string')
svariables['basis1'] = (basis1,'string')
svariables['basis2'] = (basis2,'string')
svariables['lowres1'] = (lowres1,'int')
svariables['lowres2'] = (lowres2,'int')
svariables['highres1'] = (highres1,'int')
svariables['highres2'] = (highres2,'int')
svariables['ebasis1'] = (ebasis1,'string')
svariables['ebasis2'] = (ebasis2,'string')



error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()
else:
    error=align_filter.check_align(variables)
    if(len(error) != 0):
        print 'error = ',error
        sys.exit()

import time; start = time.time()
txtQueue = multiprocessing.JoinableQueue()

align = align.align()
align.main(variables, txtQueue)

this_text = txtQueue.get(True, timeout=0.1)
print "time used: ",time.time()-start




