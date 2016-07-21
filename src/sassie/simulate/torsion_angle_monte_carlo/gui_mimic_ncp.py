'''
Driver method to run the monte_carlo module using mixed DNA and protein moves
with several options for a nucleosomecore particle (NCP)
$Id: gui_mimic_ncp.py 3249 2016-06-27 03:44:03Z schowell $
'''
import sys
import multiprocessing
import sassie.interface.input_filter as input_filter
# import sassie.simulate.monte_carlo.monte_carlo as monte_carlo
import monte_carlo_ue as monte_carlo

svariables = {}

################################# user input ##################################
################################# user input ##################################
################################# user input ##################################

# input files
pdbfile='c36_w601_ncp_min.pdb'
psffile='c36_w601_ncp.psf'

# output files
dcdfile='w601_ncp.dcd'

# minimization parameters
psf_flag = False # openmm minimization does not work with dsDNA
max_steps = '5000'
energy_convergence = '1.0'
step_size = '0.002'

# setup flexible regions
basis_string_array = []
post_basis_string_array = []
# flex_option = 'mixed'
# flex_option = 'tails'
# flex_option = '40bp'
flex_option = '40arm1'
if flex_option == '10bp':
    runname = 'run_' + pdbfile[:-4] + '_10bp'
    basis_string_array.append(
        '(segname DNA1 and resid > 15 and resid < 26) or (segname DNA2 and resid > 328 and resid < 339)'
    )
    basis_string_array.append(
        '(segname DNA1 and resid > 150 and resid < 161) or (segname DNA2 and resid > 193 and resid < 204)'
    )

    post_basis_string_array.append(
        '(segname DNA1 and resid > 25) or (segname DNA2 and resid < 329) or segname 1H2A or segname 2H2A or segname 1H2B or segname 2H2B or segname 1H3 or segname 2H3 or segname 1H4 or segname 2H4'
    )
    post_basis_string_array.append(
        '(segname DNA1 and resid 161) or(segname DNA2 and resid 193)'
    )

    rotation_type_array = ['double_stranded_nucleic_torsion',
                           'double_stranded_nucleic_torsion']

elif flex_option == '10arm1':
    runname = 'run_' + pdbfile[:-4] + '_arm1_10bp'
    basis_string_array.append(
        '(segname DNA1 and resid > 15 and resid < 26) or (segname DNA2 and resid > 328 and resid < 339)'
    )

    post_basis_string_array.append(
        '(segname DNA1 and resid > 25) or (segname DNA2 and resid < 329) or segname 1H2A or segname 2H2A or segname 1H2B or segname 2H2B or segname 1H3 or segname 2H3 or segname 1H4 or segname 2H4'
    )

    rotation_type_array = ['double_stranded_nucleic_torsion',
                           'double_stranded_nucleic_torsion']

elif flex_option == '40arm1':
    runname = 'run_' + pdbfile[:-4] + '_arm1_40bp'
    basis_string_array.append(
        '(segname DNA1 and resid >= 15 and resid <= 54) or (segname DNA2 and resid >= 300 and resid <= 339)'
    )

    post_basis_string_array.append(
        '(segname DNA1 and resid > 54) or (segname DNA2 and resid < 300) or segname 1H2A or segname 2H2A or segname 1H2B or segname 2H2B or segname 1H3 or segname 2H3 or segname 1H4 or segname 2H4'
    )

    rotation_type_array = ['double_stranded_nucleic_torsion',
                           'double_stranded_nucleic_torsion']

elif flex_option == '40bp':
    runname = 'run_' + pdbfile[:-4] + '_40bp'
    basis_string_array.append(
        '(segname DNA1 and resid >= 15 and resid <= 54) or (segname DNA2 and resid >= 300 and resid <= 339)'
    )
    basis_string_array.append(
        '(segname DNA1 and resid >= 121 and resid <= 160) or (segname DNA2 and resid >= 194 and resid <= 233)'
    )

    post_basis_string_array.append(
        '(segname DNA1 and resid > 54) or (segname DNA2 and resid < 300) or segname 1H2A or segname 2H2A or segname 1H2B or segname 2H2B or segname 1H3 or segname 2H3 or segname 1H4 or segname 2H4'
    )
    post_basis_string_array.append(
        '(segname DNA1 and resid 161) or (segname DNA2 and resid 193)'
    )

    rotation_type_array = ['double_stranded_nucleic_torsion',
                           'double_stranded_nucleic_torsion']

elif flex_option == 'H2A':
    runname = 'run_' + pdbfile[:-4] + '_tails'
    basis_string_array.append(
        'segname 1H2A and resid < 20 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H2A and resid < 20 and resid > 1'
    )
    post_basis_string_array.append(
        'segname 1H2A and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H2A and resid <= 1'
    )
    rotation_type_array = ['protein_backbone_torsion',
                           'protein_backbone_torsion']


elif flex_option == 'H2B':
    # pivot is found to be inside a loop
    runname = 'run_' + pdbfile[:-4] + '_tails'
    basis_string_array.append(
        'segname 1H2B and resid < 16 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H2B and resid < 16 and resid > 1'
    )
    post_basis_string_array.append(
        'segname 1H2B and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H2B and resid <= 1'
    )
    rotation_type_array = ['protein_backbone_torsion',
                           'protein_backbone_torsion']

elif flex_option == 'H3':
    # pivot is found to be inside a loop
    runname = 'run_' + pdbfile[:-4] + '_tails'
    # basis_string_array.append('segname 1H3 and resid < 10 and resid > 1')
    basis_string_array.append(
        'segname 1H3 and resid < 16 and resid > 1'
    )
    # basis_string_array.append('segname 1H3 and resid < 19 and resid > 15')
    basis_string_array.append(
        'segname 2H3 and resid < 16 and resid > 1'
    )
    post_basis_string_array.append(
        'segname 1H3 and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H3 and resid <= 1'
    )
    rotation_type_array = ['protein_backbone_torsion',
                           'protein_backbone_torsion']
    # rotation_type_array = ['protein_backbone_torsion']

elif flex_option == 'H4':
    runname = 'run_' + pdbfile[:-4] + '_tails'
    basis_string_array.append(
        'segname 1H4 and resid < 31 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H4 and resid < 31 and resid > 1'
    )
    post_basis_string_array.append(
        'segname 1H4 and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H4 and resid <= 1'
    )
    rotation_type_array = ['protein_backbone_torsion',
                           'protein_backbone_torsion']

elif flex_option == 'tails':
    runname = 'run_' + pdbfile[:-4] + '_tails'
    basis_string_array.append(
        'segname 1H2A and resid < 20 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H2A and resid < 20 and resid > 1'
    )
    basis_string_array.append(
        'segname 1H2B and resid < 27 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H2B and resid < 27 and resid > 1'
    )
    basis_string_array.append(
        'segname 1H3 and resid < 40 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H3 and resid < 40 and resid > 1'
    )
    basis_string_array.append(
        'segname 1H4 and resid < 31 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H4 and resid < 31 and resid > 1'
    )

    post_basis_string_array.append(
        'segname 1H2A and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H2A and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 1H2B and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H2B and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 1H3 and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H3 and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 1H4 and resid <= 1'
    )
    post_basis_string_array.append(
        'segname 2H4 and resid <= 1'
    )

    rotation_type_array = ['protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion']

elif flex_option == 'mixed':
    runname = 'run_' + pdbfile[:-4] + '_mixed_MC'
    basis_string_array.append(
        '(segname DNA1 and resid > 14 and resid < 55) or (segname DNA2 and resid > 299 and resid < 340)'
    )
    basis_string_array.append(
        '(segname DNA1 and resid > 120 and resid < 161) or (segname DNA2 and resid > 193 and resid < 234)'
    )
    basis_string_array.append(
        'segname 1H2A and resid < 20 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H2A and resid < 20 and resid > 1'
    )
    basis_string_array.append(
        'segname 1H2B and resid < 27 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H2B and resid < 27 and resid > 1'
    )
    basis_string_array.append(
        'segname 1H3 and resid < 40 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H3 and resid < 40 and resid > 1'
    )
    basis_string_array.append(
        'segname 1H4 and resid < 31 and resid > 1'
    )
    basis_string_array.append(
        'segname 2H4 and resid < 31 and resid > 1'
    )

    post_basis_string_array.append(
        '(segname DNA1 and resid > 54) or (segname DNA2 and resid < 300) or segname 1H2A or segname 2H2A or segname 1H2B or segname 2H2B or segname 1H3 or segname 2H3 or segname 1H4 or segname 2H4'
    )
    post_basis_string_array.append(
        '(segname DNA1 and resid 161) or (segname DNA2 and resid 193)'
    )
    post_basis_string_array.append(
        'segname 1H2A and resid < 2'
    )
    post_basis_string_array.append(
        'segname 2H2A and resid < 2'
    )
    post_basis_string_array.append(
        'segname 1H2B and resid < 2'
    )
    post_basis_string_array.append(
        'segname 2H2B and resid < 2'
    )
    post_basis_string_array.append(
        'segname 1H3 and resid < 2'
    )
    post_basis_string_array.append(
        'segname 2H3 and resid < 2'
    )
    post_basis_string_array.append(
        'segname 1H4 and resid < 2'
    )
    post_basis_string_array.append(
        'segname 2H4 and resid < 2'
    )

    rotation_type_array = ['double_stranded_nucleic_torsion',
                           'double_stranded_nucleic_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion',
                           'protein_backbone_torsion']
    delta_theta_array = '10.0, 10.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0'
else:
    assert False, 'ERROR: No such flex_option. Review script!'

n_flex_regions = len(basis_string_array)
number_of_flexible_regions = str(n_flex_regions)
rotation_direction_array = ['reverse'] * n_flex_regions # irrelevant

dta = delta_theta_array = '10.0'
for i in xrange(n_flex_regions):
    dta += ', ' + delta_theta_array
delta_theta_array = dta

if not delta_theta_array:
    delta_theta_array = delta_theta = '10.0'
    for i in xrange(n_flex_regions - 1):
        delta_theta_array += ', ' + delta_theta

overlap_basis = 'heavy'
temperature = '300.0'

trial_steps = '99'
goback = '50'

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
        print 'Failed to move output: %s' % cmd
except:
    print 'Aligment of NCP failed'