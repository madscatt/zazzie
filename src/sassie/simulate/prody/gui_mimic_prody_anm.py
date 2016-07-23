import sys

import sassie.simulate.prody.prody_anm as prody_anm
import sassie.interface.input_filter as input_filter
import sassie.interface.prody.prody_filter as prody_filter
import multiprocessing

svariables = {}

# BEGIN USER EDIT
# BEGIN USER EDIT
# BEGIN USER EDIT

runname = 'run_0'
pdbfile = '../../../developer_files_for_testing/prody/hivr.pdb'
number_modes = '5'  # number of normal modes to compute
number_conformations_samp = '50'  # number of conformations to generate by random sampling of modes
number_steps_traverse = '10'  # number of steps to tranverse each mode in both diretcions
rmsd_conformations_samp = '1.0'  # average RMSD of randomly sampled conformations with respect to initial conformation
rmsd_traverse = '1.5'  # maximum RMSD of conformations for trajectory from traversed mode with respect to initial conformation
advanced_usage = '0'  # advanced usage option: 0=no; 1=yes
advanced_usage_cmd = ' '  # user supplied ProDy command if advanced_usage=1

# END USER EDIT
# END USER EDIT
# END USER EDIT

svariables['runname'] = (runname, 'string')
svariables['pdbfile'] = (pdbfile, 'string')
svariables['number_modes'] = (number_modes, 'int')
svariables['number_conformations_samp'] = (number_conformations_samp, 'int')
svariables['number_steps_traverse'] = (number_steps_traverse, 'int')
svariables['rmsd_conformations_samp'] = (rmsd_conformations_samp, 'float')
svariables['rmsd_traverse'] = (rmsd_traverse, 'float')
svariables['advanced_usage'] = (advanced_usage, 'int')
svariables['advanced_usage_cmd'] = (advanced_usage_cmd, 'string')

error, variables = input_filter.type_check_and_convert(svariables)
if(len(error) > 0):
    print 'error = ', error
    sys.exit()

error = prody_filter.check_prody(variables)
if(len(error) > 0):
    print 'error = ', error
    sys.exit()

txtQueue = multiprocessing.JoinableQueue()

simulation = prody_anm.simulation()
simulation.main(variables, txtQueue)

# the following line is for the 1.0 execution only

this_text = txtQueue.get(True, timeout=0.1)
