'''
Driver method to run the monte_carlo module using mixed DNA and protein moves
with simple options for a nucleosomecore particle (NCP)
$Id: gui_mimic_simplest_ncp.py 3103 2016-04-26 20:15:12Z schowell $
'''
import sys
import multiprocessing
import sassie.interface.input_filter as input_filter
import sassie.simulate.torsion_angle_monte_carlo.monte_carlo as monte_carlo

svariables = {}

################################# user input ##################################
################################# user input ##################################
################################# user input ##################################

# input files
pdbfile='../../../developer_files_for_testing/torsion_angle_monte_carlo/c36_w601_ncp_min.pdb'
psffile='../../../developer_files_for_testing/torsion_angle_monte_carlo/c36_w601_ncp.psf'

# output file
dcdfile='ncp_test.dcd'

# run parameters
runname = 'run_ncp_test'
trial_steps = '10'
goback = '1'
temperature = '300.0'
n_flex_regions = 2
number_of_flexible_regions = str(n_flex_regions)
rotation_direction_array = ['reverse'] * n_flex_regions # irrelevant for DNA
delta_theta_array = '10.0, 30.0'
overlap_basis = 'heavy'

# setup flexible regions
basis_string_array = []
post_basis_string_array = []
basis_string_array.append(
    '(segname DNA1 and resid > 150 and resid < 161) or (segname DNA2 and resid > 193 and resid < 204)'
)
basis_string_array.append(
    'segname 1H3 and resid < 40 and resid > 1'
)
post_basis_string_array.append(
    '(segname DNA1 and resid 161) or (segname DNA2 and resid 193)'
)
post_basis_string_array.append(
    'segname 1H3 and resid <= 1'
)
rotation_type_array = ['double_stranded_nucleic_torsion',
                       'protein_backbone_torsion']

# hard coded parameters
psf_flag = True

# openmm parameters no longer used but still required (irrelevant)
max_steps = '5000'
energy_convergence = '1.0'
step_size = '0.002'

# advanced input (may not be fully implemented)
low_rg_cutoff = '0'
high_rg_cutoff = '400.0'

z_flag = False
z_cutoff = '0.0'

constraint_flag = False
constraint_file = 'constraints.txt'

directed_mc = '0'
nonbondflag = '0' # not sure what this is for
seed = '0,123'   # set this to '1,123' ('0,123') to (not) set the seed

############################### end user input ################################
############################### end user input ################################
############################### end user input ################################

svariables['runname'] = (runname, 'string')
svariables['dcdfile'] = (dcdfile, 'string')
svariables['pdbfile'] = (pdbfile, 'string')
svariables['psffile'] = (psffile, 'string')
svariables['psf_flag'] = (psf_flag, 'string')

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
svariables['seed'] = (seed,  'int_array')

error, variables = input_filter.type_check_and_convert(svariables)

assert not error, 'ERROR: %s' % error

txtQueue=multiprocessing.JoinableQueue()
simulation = monte_carlo.simulation()
simulation.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)

# perform alignment
try:
    import os
    import subprocess
    import sassie.util.file_utils as file_utils
    import sassie.tools.align_driver as align_driver

    dcd = os.path.join(runname, 'monte_carlo', dcdfile)
    assert os.path.exists(dcd), 'no such file: %s' % dcd
    align_basis = (
        '((name[i] == "CA") and (segname[i] == "1H2A") and (resid[i] > 105) and (resid[i] < 115))'
    )
    inputs = align_driver.inputs()
    inputs.path = ''
    inputs.goal_filter = align_basis
    inputs.move_filter = align_basis
    inputs.goal = pdbfile
    inputs.ref = pdbfile
    inputs.move = dcd
    inputs.out = dcd.replace('.dcd', '_al.dcd')
    file_utils.mkdir_p(os.path.split(inputs.out)[0])
    align_driver.align(inputs)
    cmd = 'mv %s %s' % (inputs.out, inputs.move)
    return_code = subprocess.call(cmd, shell=True)
    if return_code:
        print('Failed to move output: %s' % cmd)
except:
    print('Aligment of NCP failed')
