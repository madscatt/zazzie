from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sasmol.sasmol as sasmol
from compiledUtils._cLoops import ffi, lib
from lennard_gofr import *
import matplotlib as mpl
from GV import *
import cPickle as pickle
import datetime
from multiprocessing import Pool


print(datetime.datetime.now())
startTime = datetime.datetime.now()

pdb_NAME = 'Data/small/final-arg-5.pdb'
dcd_NAME = 'Data/small/run_1.dcd'

startFrame = 2
endFrame   =  -1 #-1 = use all
NUM_Q      = 250
START_Q    = -1
END_Q      = 1.6
NUM_GV       = 27


mol = sasmol.SasMol(0)
mol.read_pdb(pdb_NAME)

mol.read_dcd(dcd_NAME)


print('here')
gv = GV(NUM_GV).get_golden_vectors()
Q_list = np.logspace(START_Q,END_Q,NUM_Q)


coor=mol.coor()
num = len(coor[0])

if endFrame==-1:
    endFrame=len(coor)-1
frames = np.arange(startFrame,endFrame,1)
print('here')
def cast_matrix(matrix, ffi):
    ap = ffi.new("double* [%d]" % (matrix.shape[0]))
    ptr = ffi.cast("double *", matrix.ctypes.data)
    for i in range(matrix.shape[0]):
        ap[i] = ptr + i*matrix.shape[1]
    return ap
def process_frame(frame):
    I = np.zeros(len(Q_list))
    for i,Q in enumerate(Q_list):
        I_tmp = 0
        for g in gv:
            q=g*Q
            cast_coor = cast_matrix(coor[frame],ffi)
            cast_q = ffi.cast('double*',q.ctypes.data)
            I_tmp += lib.sQ(cast_coor,cast_q,num,num)
        I[i] = I_tmp/len(gv)
    return I
pool = Pool(processes=75)              # process per core
I_mp=pool.map(process_frame, frames)

def getOutputDir():
    '''
    Returns a path to a time stamped directory under the Outputs folder as a string
    '''
    import os
    from time import strftime
    directory = "Outputs/"+strftime("%Y-%m-%d_%H-%M")
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory+"/"

endTime = datetime.datetime.now()
outdir = getOutputDir()
f = open(outdir+'EXP_info.txt','w')
f.write('Exp Method')
f.write('\n\nSTART_Q,\t'+str(START_Q))
f.write('\nEND_Q,\t'+str(END_Q))
f.write('\nNUM_Q,\t'+str(NUM_Q))
f.write('\nNUM_GV,\t'+str(NUM_GV))

f.write('\n\nstartFrame,\t'+str(startFrame))
f.write('\nendFrame,\t'+str(endFrame))
f.write('\n\npdb_NAME,\t'+str(pdb_NAME))
f.write('\ndcd_NAME,\t'+str(dcd_NAME))

f.write('\n\nstartTime,\t'+str(startTime))
f.write('\nendTime,\t'+str(endTime))
f.close()
pickle.dump( I_mp, open( outdir+"multiFrame-"+str(frames[0])+'-'+str(frames[-1])+'_'+str(START_Q)+
		'-'+str(END_Q)+'_'+str(NUM_Q), "wb" ) )
np.save(outdir + 'outPutI-Q'+str(NUM_Q),I_mp)
print(datetime.datetime.now())
