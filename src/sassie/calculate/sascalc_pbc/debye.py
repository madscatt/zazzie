from __future__ import division
import numpy as np
import sasmol.sasmol as sasmol
import cPickle as pickle
import numexpr as ne
from multiprocessing import Pool, cpu_count
import datetime


def pairwise_numpy(X):
    return np.sqrt(((X[:, None, :] - X) ** 2).sum(-1))


def neSinc(x):
    a = ne.evaluate("sin(x)/x")
    a[np.isnan(a)] = 1
    return a


def process_frame(frame):
    I = np.zeros(len(Q_list))
    pw = pairwise_numpy(coor[frame])
    for i, Q in enumerate(Q_list):
        I[i] = np.sum(neSinc(Q * pw))
    return I
# Define parameters
folderName = '/home/data/sascalc_pbc/ellipsoids_simulation/simulations/'
pdb_NAME = 'LJ_sphere_monomer_2095bar/final.pdb'
dcd_NAME = 'LJ_sphere_monomer_2095bar/run_1.dcd'

startFrame = 1
endFrame = -1  # -1 = use all
NUM_Q = 250
START_Q = -1
END_Q = 1.6

# Load pdb + dcd
mol = sasmol.SasMol(0)
mol.read_pdb(pdb_NAME)
mol.read_dcd(dcd_NAME)
coor = mol.coor()
Q_list = np.logspace(START_Q, END_Q, NUM_Q)
if endFrame == -1:
    endFrame = len(coor) - 1
frames = np.arange(startFrame, endFrame, 1)
# frames = np.arange(600,604,1) # for debugging


# Do calculation
startTime = datetime.datetime.now()
print(startTime)
pool = Pool(processes=cpu_count())              # process per core
I_mp = pool.map(process_frame, frames)
endTime = datetime.datetime.now()
print(endTime)
minutes = (endTime - startTime).seconds / 60
print(minutes)

# Ouput Results


def getOutputDir():
    '''
    Returns a path to a time stamped directory under the Outputs folder as a string
    '''
    import os
    from time import strftime
    directory = "Outputs/" + strftime("%Y-%m-%d_%H-%M")
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory + "/"

outdir = getOutputDir()
f = open(outdir + 'DEBYE_info.txt', 'w')
f.write('DEBYE METHOD')
f.write('\n\nSTART_Q,\t' + str(START_Q))
f.write('\nEND_Q,\t' + str(END_Q))
f.write('\nNUM_Q,\t' + str(NUM_Q))

f.write('\n\nstartFrame,\t' + str(startFrame))
f.write('\nendFrame,\t' + str(endFrame))
f.write('\n\npdb_NAME,\t' + str(pdb_NAME))
f.write('\ndcd_NAME,\t' + str(dcd_NAME))

f.write('\n\nstartTime,\t' + str(startTime))
f.write('\nendTime,\t' + str(endTime))
f.write('\nminutes,\t' + str(minutes))
f.close()
pickle.dump(I_mp, open(outdir + "multiFrame-" + str(frames[0]) + '-' + str(frames[-1]) + '_' + str(START_Q) +
                       '-' + str(END_Q) + '_' + str(NUM_Q), "wb"))
np.save(outdir + 'outPutI-Q' + str(NUM_Q), I_mp)
np.save(outdir + 'Q_list', Q_list)
