'''
Calculates the scattering profile using the Debye formula.  Requires a pdb
and dcd respectively in the pdb_fname and dcd_fname variables. Performs the
Debye calculation in parallal over frames via the debye.process_frame
function, which acutally does the work (parallelized over the frames of
a given dcd). Only takes one dcd/pdb input pair.
'''
import logging
import multiprocessing
import os
import time
import numpy as np
from functools import partial

import debye
import sasmol.sasmol as sasmol

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

# Ouput Results
def getOutputDir():
    '''
    Returns a path to a time stamped directory
    under the `outputs/` folder as a string
    '''
    from time import strftime
    directory = "outputs/{}".format(strftime("%Y-%m-%d_%H-%M"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory + "/"


# Define parameters
path = '/home/data/sascalc_pbc/ellipsoids_simulation/simulations/'
pdb_fname = 'LJ_sphere_monomer_13.4A/final.pdb'
dcd_fname = 'LJ_sphere_monomer_13.4A/run_13.dcd'

start_frame = 0 # minimum of ...?
end_frame = -1  # -1 = use all
n_q = 100
q_min = 1 # 1/A
q_max = 3 # 1/A

sigma = 3.405

# Load pdb + dcd
mol = sasmol.SasMol(0)
mol.read_pdb(path + pdb_fname)
mol.read_dcd(path + dcd_fname)
coor = mol.coor() * sigma
q = np.logspace(np.log10(q_min), np.log10(q_max), n_q)
if end_frame == -1:
    end_frame = len(coor) - 1
input_coors = coor[start_frame:end_frame]


# Do calculation
tic = time.time()
# I_mp = pool.map(simple_debye, frames)
# I_mp = pool.map(glatter_kratky, (q, coor[frames])
debug = True

if debug:
    # I = debye.glatter_kratky(q, coor[690]) # single frame debug
    # np.savetxt('iq_690.dat', I)
    I = debye.simple_debye(q, coor[690]) # single frame debug
    np.savetxt('iq_ian_690.dat', I)


else:
    # setup a worker pool with 1 process/core
    func = partial(debye.glatter_kratky, q)
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    I_all = pool.map(func, input_coors)
    I_all = np.array(I_all)
    logging.debug('seconds: {} \nminutes:{}'.format(toc, toc/60))

    # print(I_all)
    I_mean = I_all.mean(axis=0)

    I_fname = 'iq_{}to{}.dat'.format(start_frame, end_frame)
    np.savetxt(I_fname, I_mean)

toc = time.time() - tic

do_plot = False
do_save = False

if do_plot:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(I_mean[:, 0], I_mean[:, 1], label='mean I')
    plt.plot(I_all[0, :, 0], I_all[0, :, 1], label='all I')
    plt.xlabel(r'Q ($\AA$)')
    plt.ylabel(r'I(Q)')
    plt.legend()
    plt.savefig('I_comparison.png', dpi=400, bbox_inches='tight')
    plt.close('all')

if do_save:
    outdir = getOutputDir()
    with open(outdir + 'DEBYE_info.txt', 'w') as f:
        f.write('DEBYE METHOD')
        f.write('\n\nSTART_Q,\t' + str(q_min))
        f.write('\nEND_Q,\t' + str(q_max))
        f.write('\nNUM_Q,\t' + str(n_q))

        f.write('\n\nstartFrame,\t' + str(start_frame))
        f.write('\nendFrame,\t' + str(end_frame))
        f.write('\n\npdb_NAME,\t' + str(pdb_fname))
        f.write('\ndcd_NAME,\t' + str(dcd_fname))

        f.write('\n\nstartTime,\t' + str(tic))
        f.write('\nendTime,\t' + str(tic + toc))
        f.write('\nminutes,\t' + str(minutes))

    # pickle.dump(I_mp, open(outdir + "multiFrame-" + str(frames[0]) + '-' + str(frames[-1]) + '_' + str(START_Q) +
                           # '-' + str(END_Q) + '_' + str(NUM_Q), "wb"))
    np.save(os.path.join(outdir, 'outPutI-Q{}'.format(n_q)), I_all)
    np.save(os.path.join(outdir, 'Q_list'), q)

logging.debug('\m/ >.< \m/')