'''
Driver method to run the monte_carlo module
using DNA moves with 60 bp of dsDNA
$Id: gui_mimic_BDNA.py 3282 2016-07-15 05:08:59Z schowell $
'''

import sys
import sassie.simulate.torsion_angle_monte_carlo.monte_carlo as monte_carlo
import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

# input files
pdbfile = '../../../developer_files_for_testing/torsion_angle_monte_carlo/c36_dsDNA60_min.pdb'
psffile = '../../../developer_files_for_testing/torsion_angle_monte_carlo/c36_dsDNA60.psf'

# output files
dcdfile = 'dsDNA60.dcd'

# minimization parameters
psf_flag = False # openmm minimization does not work with dsDNA
max_steps = '5000'
energy_convergence = '1.0'
step_size = '0.002'

# setup flexible regions
basis_string_array = []
post_basis_string_array = []
case = 3
overlap_basis = 'heavy'
if 1 == case:
    runname = 'run_dsDNA_one_flexible_regions'
    number_of_flexible_regions = '1'

    basis_string_array.append(
        '(segname DNA1 and resid >= 11 and resid <= 20) or '
        '(segname DNA2 and resid >= 101 and resid <= 110)'
    )
    delta_theta_array = '10.0'
    rotation_type_array = ['double_stranded_nucleic_torsion']
    rotation_direction_array = ['forward'] # won't use

    post_basis_string_array.append(
        '(segname DNA1 and resid > 20) or '
        '(segname DNA2 and resid < 101)'
    )
elif -1 == case:
    runname = 'run_dsDNA_one_reverse_flexible_regions'
    number_of_flexible_regions = '1'

    basis_string_array.append(
        '(segname DNA2 and resid >= 101 and resid <= 110) or '
        '(segname DNA1 and resid >= 11 and resid <= 20)'
    )
    delta_theta_array = '10.0'
    rotation_type_array = ['double_stranded_nucleic_torsion']
    rotation_direction_array = ['forward'] # won't use

    post_basis_string_array.append(
        '(segname DNA2 and resid > 110) or '
        '(segname DNA1 and resid < 11)'
    )
elif 2 == case:
    runname = 'run_dsDNA_two_flexible_regions'
    number_of_flexible_regions = '2'

    basis_string_array.append('(segname DNA1 and resid > 10 and resid < 21) or '
                              '(segname DNA2 and resid > 100 and resid < 111)')
    basis_string_array.append('(segname DNA1 and resid > 30 and resid < 46) or '
                              '(segname DNA2 and resid > 75 and resid < 91)')
    delta_theta_array = '10.0, 10.0'
    rotation_type_array = ['double_stranded_nucleic_torsion',
                           'double_stranded_nucleic_torsion']
    rotation_direction_array = ['forward', 'forward'] # won't use

    post_basis_string_array.append('(segname DNA1 and resid > 20) or '
                                   '(segname DNA2 and resid < 101)')
    post_basis_string_array.append('(segname DNA1 and resid > 45) or '
                                   '(segname DNA2 and resid < 76)')

elif 3 == case:
    runname = 'run_dsDNA_all_flexible'
    number_of_flexible_regions = '1'

    basis_string_array.append(
        '(segname DNA1 and resid > 1 and resid < 60) or (segname DNA2 and resid > 61 and resid < 120)'
    )
    delta_theta_array = '10.0'
    rotation_type_array = ['double_stranded_nucleic_torsion']
    rotation_direction_array = ['forward']

    post_basis_string_array.append(
        '(segname DNA1 and resid > 59) or (segname DNA2 and resid < 62)'
    )

elif 4 == case:
    runname = 'run_dsDNA_test_one_flexible'
    number_of_flexible_regions = '1'

    basis_string_array.append('(segname DNA1 and resid > 10 and resid < 21) or '
                              '(segname DNA2 and resid > 100 and resid < 111)')
    delta_theta_array = '10.0'
    rotation_type_array = ['double_stranded_nucleic_torsion']
    rotation_direction_array = ['forward']

    post_basis_string_array.append('(segname DNA1 and resid > 20) or '
                                   '(segname DNA2 and resid < 101)')

temperature = '300.0'
trial_steps = '999'
goback = '1'

low_rg_cutoff = '0'
high_rg_cutoff = '400.0'

z_flag = False
z_cutoff = '0.0'

constraint_flag = False
constraint_file = 'constraints.txt'

directed_mc = '0'

nonbondflag = '0' # not sure what this is for
seed = '0, 123'  # set this to '1,123' if you want to set the seed or '0,123' if not

#### end user input ####
#### end user input ####
#### end user input ####


svariables['runname']                    = (runname, 'string')
svariables['dcdfile']                    = (dcdfile, 'string')
svariables['pdbfile']                    = (pdbfile, 'string')
svariables['psffile']                    = (psffile, 'string')
svariables['psf_flag']                   = (psf_flag, 'string')

svariables['max_steps']                  = (max_steps, 'int')
svariables['energy_convergence']         = (energy_convergence, 'float')
svariables['step_size']                  = (step_size, 'float')

svariables['number_of_flexible_regions'] = (number_of_flexible_regions, 'int')
svariables['basis_string_array']         = (basis_string_array, 'string')
svariables['delta_theta_array']          = (delta_theta_array, 'float_array')
svariables['rotation_type_array']        = (rotation_type_array, 'string')
svariables['rotation_direction_array']   = (rotation_direction_array, 'string')
svariables['overlap_basis']              = (overlap_basis, 'string')
svariables['post_basis_string_array']    = (post_basis_string_array, 'string')
svariables['temperature']                = (temperature, 'float')
svariables['trial_steps']                = (trial_steps, 'int')
svariables['goback']                     = (goback, 'int')
svariables['directed_mc']                = (directed_mc, 'float')

svariables['low_rg_cutoff']              = (low_rg_cutoff, 'float')
svariables['high_rg_cutoff']             = (high_rg_cutoff, 'float')

svariables['z_flag']                     = (z_flag, 'boolean')
svariables['z_cutoff']                   = (z_cutoff, 'float')

svariables['constraint_flag']            = (constraint_flag, 'boolean')
svariables['constraint_file']            = (constraint_file, 'string')

svariables['nonbondflag']                = (nonbondflag, 'int')
svariables['seed']                       = (seed,  'int_array')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print('error = ', error)
    sys.exit()

txtQueue=multiprocessing.JoinableQueue()
simulation = monte_carlo.simulation()
simulation.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)

