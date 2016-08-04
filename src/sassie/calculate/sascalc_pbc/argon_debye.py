import datetime
from multiprocessing import Pool, cpu_count

# Define parameters
path = '/home/data/sascalc_pbc/ellipsoids_simulation/simulations/'
pdb_fname = 'LJ_sphere_monomer_13.4A/final.pdb'
dcd_fname = 'LJ_sphere_monomer_13.4A/run_13.dcd'

start_frame = 690
end_frame = 691  # -1 = use all
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
frames = np.arange(start_frame, end_frame, 1)
# frames = np.arange(600,604,1) # for debugging


# Do calculation
tic = datetime.datetime.now()
print(tic)

pool = Pool(processes=cpu_count())              # process per core

# I_mp = pool.map(process_frame, frames)
# I_mp = pool.map(process_frame, frames)

I = glatter_kratky(q, coor, 690)

print(I)



toc = datetime.datetime.now()
print(toc)
minutes = (toc - tic).seconds / 60
print(minutes)


do_save = False

if do_save:

    outdir = getOutputDir()
    f = open(outdir + 'DEBYE_info.txt', 'w')
    f.write('DEBYE METHOD')
    f.write('\n\nSTART_Q,\t' + str(q_min))
    f.write('\nEND_Q,\t' + str(q_max))
    f.write('\nNUM_Q,\t' + str(n_q))

    f.write('\n\nstartFrame,\t' + str(start_frame))
    f.write('\nendFrame,\t' + str(end_frame))
    f.write('\n\npdb_NAME,\t' + str(pdb_fname))
    f.write('\ndcd_NAME,\t' + str(dcd_fname))

    f.write('\n\nstartTime,\t' + str(tic))
    f.write('\nendTime,\t' + str(toc))
    f.write('\nminutes,\t' + str(minutes))
    f.close()
    # pickle.dump(I_mp, open(outdir + "multiFrame-" + str(frames[0]) + '-' + str(frames[-1]) + '_' + str(START_Q) +
                           # '-' + str(END_Q) + '_' + str(NUM_Q), "wb"))
    np.save(outdir + 'outPutI-Q' + str(n_q), I_mp)
    np.save(outdir + 'Q_list', q)