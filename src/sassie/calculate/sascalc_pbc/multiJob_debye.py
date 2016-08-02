'''
same as debye.py but modified to loop over multiple input pdbs/dcds.
'''
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
folder = '/home/data/sascalc_pbc/ellipsoids_simulation/simulations/'

"""
#First Big run. Did July 22.
conds = [
        ('LJ_sphere_monomer_high_pressure/','final.pdb','run_1.dcd'),
        ('LJ_sphere_monomer_low_pressure/','final-arg-low_p_002.pdb','run_low_p_002.dcd'),
        ('LJ_sphere_monomer_p0.9/','final.pdb','run_1.dcd'),
        ('LJ_sphere_monomer_p1/','final.pdb','run_1.dcd'),
        ('LJ_sphere_monomer_2095bar/','final.pdb','run_1.dcd'),
        ('LJ_sphere_monomer_2A/','final-arg-2.pdb','run_2.dcd'),
        ('LJ_sphere_monomer_4A/','final-arg-4.pdb','run_4.dcd'),
        ('LJ_sphere_monomer_5A/','final-arg-5.pdb','run_1.dcd'),
        ('LJ_sphere_monomer_8A/','final-arg-8.pdb','run_8.dcd'),
        ('LJ_sphere_monomer_12A/','final-arg-12.pdb','run_12.dcd'),
        ('LJ_sphere_monomer_13.4A/','final.pdb','run_13.dcd') # this may be the wrong pdb
        ]
"""
conds = [
    ('LJ_sphere_monomer_p1_p0.0213_d0.8_n5000/',
     'n5000_fcc0.80001.pdb', 'run_1.dcd'),
    ('LJ_sphere_monomer_p1_p0.0213_d3.82_n5000/', 'n5000_fcc3.82.pdb', 'run_1.dcd'),
    ('LJ_sphere_monomer_p1_p0.0213_d0.8/', 'n2457_fcc0.80001.pdb', 'run_1.dcd'),
    ('LJ_sphere_monomer_p1_p0.0213_d3.82/', 'n2457_fcc3.82.pdb', 'run_1.dcd'),
]
outDirs = []
overAllStart = datetime.datetime.now()
for cond in conds:
    print(cond[0])
    pdb_NAME = folder + cond[0] + cond[1]
    dcd_NAME = folder + cond[0] + cond[2]

    startFrame = 500
    endFrame = -1  # -1 = use all
    NUM_Q = 200
    START_Q = -1
    END_Q = 1.6

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_NAME)
    mol.read_dcd(dcd_NAME)
    # Load pdb + dcd mol = sasmol.SasMol(0) mol.read_pdb(pdb_NAME)
    # mol.read_dcd(dcd_NAME)
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
    outDirs.append(outdir)
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
    np.save(outdir + 'outPutI-Q' + str(NUM_Q), I_mp)
    np.save(outdir + 'Q_list', Q_list)
import time
f = open("Outputs/RunEnding-" + time.strftime("%Y-%m-%d_%H-%M"), 'w')
f.write("Outdirs:\n")
for i in outDirs:
    f.write(i + '\n')
end = datetime.datetime.now()
f.write('\nStart:\t{}'.format(overAllStart))
f.write('\nEnd:\t{}'.format(end))
minutes = (end - overAllStart).seconds / 60
f.write("\nminutes:\t{}".format(minutes))
f.close()
