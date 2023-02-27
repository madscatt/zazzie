'''
Full exponential calculator.
This one takes only one pdb/dcd pair. Has the same inputs as debye.py
output is I[frame][Q] as a .npy file, which can be loaded with np.load
Also outputs Q_list.npy
O(num GV * N^2)
'''
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sasmol.sasmol as sasmol
import numexpr as ne
from compiledUtils._cLoops import ffi, lib
from lennard_gofr import *
import matplotlib as mpl
from GV import *
import cPickle as pickle
import datetime
from multiprocessing import Pool, cpu_count


print(datetime.datetime.now())
startTime = datetime.datetime.now()

pdb_NAME = '/home/data/sascalc_pbc/pdbs/run_0.pdb'
dcd_NAME = '/home/data/sascalc_pbc/pdbs/run_1.dcd'

startFrame =995 
endFrame   =  -1 #-1 = use all
NUM_Q      = 100
START_Q    = -1 #units in logscale. startq = 10^(START_Q)
END_Q      = 1.6
NUM_GV       = 27
sigma=3.4

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
    '''
    process_frame by using the C code from cLoops.py
    Doesn't break memory. Could be much more optimized by passing more thing to C code that remain the same over Q's
    '''
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
def process_frame_new(frame):
    """
    modified to be purely numpy/numexpr
    Can overallocate memory but slightly nicer than passing things to the C code.
    """
    I = np.zeros(len(Q_list))
    disp = np.zeros((len(coor[frame]),len(coor[frame]),3))
    for i in xrange(len(coor[frame])):
        for j in xrange(len(coor[frame])):
            disp[i][j] = coor[frame][i]-coor[frame][j]
    reshaped = disp.reshape((mol.natoms()**2,3))*sigma # sigma from reduced units
    dotted = np.inner(gv,reshaped)
    for i,Q in enumerate(Q_list):
        inner = dotted*Q
        I[i] = ne.evaluate("sum(cos(inner))")
    I /= len(gv)
    return I
pool = Pool(processes=cpu_count())              # process per core
#I_mp=pool.map(process_frame, frames)
I_mp=pool.map(process_frame_new, frames)

def getOutputDir(base=''):
    '''
    Returns a path to a time stamped directory under the Outputs folder as a string
    '''
    import os
    from time import strftime
    directory = base+"Outputs/"+strftime("%Y-%m-%d_%H-%M")
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory+"/"

endTime = datetime.datetime.now()
outdir = getOutputDir('/home/data/')
print('output: '+outdir)
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
np.save(outdir + 'outPutI-Q' + str(NUM_Q), I_mp)
np.save(outdir + 'Q_list', Q_list)

np.save('outPutI-Q' + str(NUM_Q), I_mp)
np.save('Q_list', Q_list)
print(datetime.datetime.now())
