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

def process_Q(Q):
    return np.sum(neSinc(Q*pw))
# Define parameters
overAllStart = datetime.datetime.now()

folder = '/home/data/sascalc_pbc/'
pdb_NAME = folder + 'pdbs/run_0.pdb' 
dcd_NAME = folder + 'pdbs/run_1.dcd'

NUM_Q = 100
START_Q = -1
END_Q = 1.6

mol = sasmol.SasMol(0)
mol.read_pdb(pdb_NAME)
mol.read_dcd(dcd_NAME)
coor = mol.coor()
pw = pairwise_numpy(coor[0])
Q_list = np.logspace(START_Q, END_Q, NUM_Q)

# Do calculation
startTime = datetime.datetime.now()
print(startTime)
pool = Pool(processes=cpu_count())              # process per core
I_mp = pool.map(process_Q, Q_list)
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

f.write('\n\npdb_NAME,\t' + str(pdb_NAME))
f.write('\ndcd_NAME,\t' + str(dcd_NAME))

f.write('\n\nstartTime,\t' + str(startTime))
f.write('\nendTime,\t' + str(endTime))
f.write('\nminutes,\t' + str(minutes))
f.close()
np.save(outdir + 'outPutI-Q' + str(NUM_Q), I_mp)
np.save(outdir + 'Q_list', Q_list)
import time
f = open("Outputs/RunEnding-" + time.strftime("%Y-%m-%d_%H-%M"), 'w')
f.write("Outdirs:\n")
end = datetime.datetime.now()
f.write('\nStart:\t{}'.format(overAllStart))
f.write('\nEnd:\t{}'.format(end))
minutes = (end - overAllStart).seconds / 60
f.write("\nminutes:\t{}".format(minutes))
f.close()
