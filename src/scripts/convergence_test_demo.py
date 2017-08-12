pdb_fname = '/home/schowell/data/scratch/convergence_debug/ref.pdb'
dcd_fname = '/home/schowell/data/scratch/convergence_debug/run_0.dcd'

from convergence_test import calc_spatial_convergence_all

calc_spatial_convergence_all(pdb_fname, [dcd_fname])
