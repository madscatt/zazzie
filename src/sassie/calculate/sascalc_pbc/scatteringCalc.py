from __future__ import division
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import sasmol.sasmol as sasmol
from lib._cLoops import ffi, lib
from lennard_gofr import *
import matplotlib as mpl
from GV import *
import cPickle as pickle

# Style plots
plt.style.use('ggplot')
mpl.rcParams['figure.figsize'] = (16, 9)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 25
mpl.rcParams['axes.titlesize'] = 25
mpl.rcParams['figure.titlesize'] = 25
mol = sasmol.SasMol(0)
mol.read_pdb('Data/run_0.pdb')
mol.read_dcd('Data/run_1.dcd')

gc = gofr_calc(mol)


def cast_matrix(matrix, ffi):
    ap = ffi.new("double* [%d]" % (matrix.shape[0]))
    ptr = ffi.cast("double *", matrix.ctypes.data)
    for i in range(matrix.shape[0]):
        ap[i] = ptr + i * matrix.shape[1]
    return ap

startFrame = 1
endFrame = -1  # -1 = use all
NUM_Q = 100
START_Q = -1
END_Q = 1
N_GV = 35
gv = GV(N_GV).get_golden_vectors()
Q_list = np.logspace(START_Q, END_Q, NUM_Q)

import datetime
print(datetime.datetime.now())

from multiprocessing import Pool
coor = mol.coor()[startFrame:endFrame]
num = len(coor[0])
N_Q = 100


def process_frame(frame):
    I = np.zeros(len(Q_list))
    for i, Q in enumerate(Q_list):
        I_tmp = 0
        for g in gv:
            q = g * Q
            cast_coor = cast_matrix(coor[frame], ffi)
            cast_q = ffi.cast('double*', q.ctypes.data)
            I_tmp += lib.sQ(cast_coor, cast_q, num, num)
        I[i] = I_tmp / len(gv)
    return I
pool = Pool(processes=70)              # process per core
# frames = np.arange(600,1000,1)
frames = np.arange(600, 680, 1)
I_mp = pool.map(process_frame, frames)
pickle.dump(I_mp, open("multiFrame-" + str(NUM_Q) + "-allFrames", "wb"))
print(datetime.datetime.now())
