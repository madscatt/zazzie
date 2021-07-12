'''
Driver method to run the monte_carlo module
'''

import sys

# import sassie.simulate.monte_carlo.monte_carlo as monte_carlo
import monte_carlo
import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'
dcdfile = 'run_0.dcd'
pdbfile = 'systems_for_docs/gag_matrix/output_building/hiv1_gag_ma.pdb'
psffile = 'systems_for_docs/gag_matrix/output_building/hiv1_gag_ma.psf'
psf_flag = True
# psf_flag = False # sh: making this false to prevent crash

max_steps = '5000'
energy_convergence = '1.0'
step_size = '0.002'

number_of_flexible_regions = '1'

basis_string_array = []
basis_string_array.append('(resid>116 and resid < 139)')
delta_theta_array = '30.0'
rotation_type_array = ['protein_backbone_torsion']
rotation_direction_array = ['forward']
overlap_basis = 'heavy'

post_basis_string_array = []
#post_basis_string_array.append('None')
post_basis_string_array.append('resid 139')

temperature = '300.0'
trial_steps = '300'
goback = '10'

low_rg_cutoff = '0'
high_rg_cutoff = '400.0'

z_flag = False
z_cutoff = '0.0'

constraint_flag = False
constraint_file = 'constraints.txt'

directed_mc = '0'

nonbondflag = '0'
seed = '0, 123'  # set this to '1,123' if you want to set the seed or '0,123' if not

#### end user input ####
#### end user input ####
#### end user input ####


svariables['runname'] = (runname, 'string')
svariables['dcdfile'] = (dcdfile, 'string')
svariables['pdbfile'] = (pdbfile, 'string')
svariables['psffile'] = (psffile, 'string')
svariables['psf_flag'] = (psf_flag, 'boolean')

svariables['max_steps'] = (max_steps, 'int')
svariables['energy_convergence'] = (energy_convergence, 'float')
svariables['step_size'] = (step_size, 'float')

svariables['number_of_flexible_regions'] = (number_of_flexible_regions, 'int')
svariables['basis_string_array'] = (basis_string_array, 'string')
svariables['delta_theta_array'] = (delta_theta_array, 'float_array')
svariables['rotation_type_array'] = (rotation_type_array, 'string')
svariables['rotation_direction_array'] = (rotation_direction_array, 'string')
svariables['overlap_basis'] = (overlap_basis, 'string')
svariables['post_basis_string_array'] = (post_basis_string_array, 'string')
svariables['temperature'] = (temperature, 'float')
svariables['trial_steps'] = (trial_steps, 'int')
svariables['goback'] = (goback, 'int')
svariables['directed_mc'] = (directed_mc, 'float')


svariables['low_rg_cutoff'] = (low_rg_cutoff, 'float')
svariables['high_rg_cutoff'] = (high_rg_cutoff, 'float')

svariables['z_flag'] = (z_flag, 'boolean')
svariables['z_cutoff'] = (z_cutoff, 'float')

svariables['constraint_flag'] = (constraint_flag, 'boolean')
svariables['constraint_file'] = (constraint_file, 'string')

svariables['nonbondflag'] = (nonbondflag, 'int')
svariables['seed'] = (seed, 'int_array')



error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()

#import pprint; pprint.pprint(variables); exit()

txtQueue = multiprocessing.JoinableQueue()
simulation = monte_carlo.simulation()
simulation.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)


#print 'in GUI and txtOutput = ', this_text, '\n'





