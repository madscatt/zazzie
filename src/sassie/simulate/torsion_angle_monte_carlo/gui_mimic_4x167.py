#!/usr/bin/env python
#!/share/apps/bin/python
#
# Author:  Steven C. Howell
# Purpose: MC generation of 4x167_k010 structures
# Created: 24 June 2015
#
# $Id: dna_mc_driver.py 121 2015-02-05 15:15:38Z xraylab $
#
#000000001111111111222222222233333333334444444444555555555566666666667777777777
#234567890123456789012345678901234567890123456789012345678901234567890123456789

import sassie.interface.input_filter as input_filter
import multiprocessing
import numpy as np
# import sassie.simulate.monte_carlo.monte_carlo as monte_carlo
import monte_carlo_ue as monte_carlo

# pdbfile = '4x167_min.pdb'
pdbfile = '4x167_sasmol.pdb'
psffile = '4x167.psf'
dcdfile = '4x167.dcd'
trial_steps = '99'
goback = '1'
case = 'all'
# case = 3

# minimization parameters
psf_flag = False # openmm minimization does not work with dsDNA
max_steps = '5000'
energy_convergence = '1.0'
step_size = '0.002'

# setup flexible regions
bps = np.array([np.linspace(0, 693, 694), np.linspace(694, 1, 694)]).T
flex = 20
if flex == 20:
    w601 = 23 + 167 * np.arange(4)
    link = w601 + 146
elif flex == 40:
    w601 = 33 + 167 * np.arange(4)
    link = w601 + 126
elif flex == 80:
    w601 = 53 + 167 * np.arange(4)
    link = w601 + 86
else:
    assert False, 'ERROR: flex selection not yet implemented'

basis_string_array = []
post_basis_string_array = []
flex_dna = ('(segname DNA1 and resid > %d and resid < %d) or '
            '(segname DNA2 and resid < %d and resid > %d)')
post_dna = '(segname DNA1 and resid >= %d) or (segname DNA2 and resid <= %d)'
octamer1 = (' or segname 1H2A or segname 1H3'
            ' or segname 1H2B or segname 1H4'
            ' or segname 2H2A or segname 2H3'
            ' or segname 2H2B or segname 2H4')
octamer2 = (' or segname 3H2A or segname 3H3'
            ' or segname 3H2B or segname 3H4'
            ' or segname 4H2A or segname 4H3'
            ' or segname 4H2B or segname 4H4')
octamer3 = (' or segname 5H2A or segname 5H3'
            ' or segname 5H2B or segname 5H4'
            ' or segname 6H2A or segname 6H3'
            ' or segname 6H2B or segname 6H4')
octamer4 = (' or segname 7H2A or segname 7H3'
            ' or segname 7H2B or segname 7H4'
            ' or segname 8H2A or segname 8H3'
            ' or segname 8H2B or segname 8H4')
flex1 = (flex_dna % (bps[1,0], bps[w601[0],0],
                     bps[1,1], bps[w601[0],1]))
flex2 = (flex_dna % (bps[link[0],0], bps[w601[1],0],
                     bps[link[0],1], bps[w601[1],1]))
flex3 = (flex_dna % (bps[link[1],0], bps[w601[2],0],
                     bps[link[1],1], bps[w601[2],1]))
flex4 = (flex_dna % (bps[link[2],0], bps[w601[3],0],
                     bps[link[2],1], bps[w601[3],1]))
flex5 = (flex_dna % (bps[link[3],0], bps[-1,0],
                     bps[link[3],1], bps[-1,1]))

post1 = (post_dna % (bps[w601[0],0], bps[w601[0],1]))
post2 = (post_dna % (bps[w601[1],0], bps[w601[1],1]))
post3 = (post_dna % (bps[w601[2],0], bps[w601[2],1]))
post4 = (post_dna % (bps[w601[3],0], bps[w601[3],1]))
post5 = (post_dna % (bps[-1,0], bps[-1,1]))
post1 += octamer1 + octamer2 + octamer3 + octamer4
post2 += octamer2 + octamer3 + octamer4
post3 += octamer3 + octamer4
post4 += octamer4

if 'all' == case:
    runname = '4x167_all'
    basis_string_array.append(flex1)
    basis_string_array.append(flex2)
    basis_string_array.append(flex3)
    basis_string_array.append(flex4)
    basis_string_array.append(flex5)

    post_basis_string_array.append(post1)
    post_basis_string_array.append(post2)
    post_basis_string_array.append(post3)
    post_basis_string_array.append(post4)
    post_basis_string_array.append(post5)

elif 3 == case:
    runname = '4x167_l3'
    basis_string_array.append(flex3)

    post_basis_string_array.append(post3)
else:
    assert False, 'ERROR: No such flex_option. Review script!'

n_flex_regions = len(basis_string_array)
number_of_flexible_regions = str(n_flex_regions)
rotation_type_array = ['double_stranded_nucleic_torsion'] * n_flex_regions
rotation_direction_array = ['forward'] * n_flex_regions # irrelevant

dta = delta_theta_array = '10.0'
for i in xrange(n_flex_regions):
    dta += ', ' + delta_theta_array
delta_theta_array = dta

overlap_basis = 'heavy'
temperature = '300.0'

low_rg_cutoff = '0'
high_rg_cutoff = '1000.0'

z_flag = False
z_cutoff = '0.0'

constraint_flag = False
constraint_file = 'constraints.txt'

directed_mc = '0'

nonbondflag = '0' # not sure what this is for
seed = '0, 123'  # set this to '1,123' ('0,123') to (not) set the seed

svariables = {}

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
assert len(error) == 0, 'ERROR: %s' % error

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
        '((name[i] == "CA") and (segname[i] == "3H2A") and '
        '(resid[i] > 5) and (resid[i] < 115))'
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