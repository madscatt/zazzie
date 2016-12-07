
'''
GUI_MIMIC for build_utilities
'''

import sys

import tool_methods
import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'

pdb_utilities_flag = True
pdbfile = '../../../developer_files_for_testing/torsion_angle_monte_carlo/hiv1_gag.pdb'

renumber_flag = True
renumber_flag = False
renumber_output_filename = 'renumbered.pdb'
renumber_indices_flag = True
renumber_resids_flag = True

first_index = '7'
first_resid = '9'

pdb_constraints_flag = False
number_of_constraint_files = '1'
constraint_options = 'protein'
constraint_options = 'nucleic'
constraint_options = 'solute'
constraint_options = 'backbone'
constraint_options = 'heavy'
constraint_fields = 'occupancy'
constraint_fields = 'beta'
number_of_constraint_files = '2'


constraint_options = 'heavy, backbone'
constraint_fields = 'beta, beta'
constraint_resets = 'True, True'
constraint_filenames = 'constrain_heavy.pdb, constrain_backbone.pdb'

translation_rotation_flag = False
translation_rotation_output_filename = 'trans_rot.pdb'
pre_center_flag = False
translation_array = '0.0, 0.0, 0.0'
rotation_axes = 'pmi'
rotation_axes = 'cardinal'
rotation_order = 'xyz'
#rotation_order = 'xzy'
#rotation_order = 'yxz'
#rotation_order = 'yzx'
#rotation_order = 'zxy'
#rotation_order = 'zyx'
rotation_array = '0.0, 0.0, 0.0'


### i am working here

align_pmi_on_axis_flag = True
pmi_eigenvector = '3'
''' values: x, y, z '''
alignment_vector_axis = 'z'
settle_on_plane = True
''' values: x, y, z '''
plane = 'z' 

### i stop working here

fasta_utilities_flag = False
fasta_input_option = 'file'
fasta_input_option = 'sequence'
fasta_input_sequence = 'RRAGMPSCYLK'
fasta_input_file = 'test_fasta.txt'

seed = '0, 123'  # set this to '1,123' if you want to set the seed or '0,123' if not

#### end user input ####
#### end user input ####
#### end user input ####


svariables['runname'] = (runname, 'string')
svariables['pdbfile'] = (pdbfile, 'string')

svariables['pdb_utilities_flag'] = (pdb_utilities_flag, 'boolean')

svariables['renumber_flag'] = (renumber_flag, 'boolean')
svariables['renumber_output_filename'] = (renumber_output_filename, 'string')

svariables['renumber_indices_flag'] = (renumber_indices_flag, 'boolean')
svariables['first_index'] = (first_index, 'int')

svariables['renumber_resids_flag'] = (renumber_resids_flag, 'boolean')
svariables['first_resid'] = (first_resid, 'int')


svariables['pdb_constraints_flag'] = (pdb_constraints_flag, 'boolean')
svariables['number_of_constraint_files'] = (number_of_constraint_files, 'int')
svariables['constraint_options'] = (constraint_options, 'string')
svariables['constraint_filenames'] = (constraint_filenames, 'string')
svariables['constraint_fields'] = (constraint_fields, 'string')
svariables['constraint_resets'] = (constraint_resets, 'string')


svariables['translation_rotation_flag'] = (translation_rotation_flag, 'boolean')
svariables['translation_rotation_output_filename'] = (translation_rotation_output_filename, 'string')
svariables['pre_center_flag'] = (pre_center_flag, 'boolean')
svariables['translation_array'] = (translation_array, 'float_array')
svariables['rotation_axes'] = (rotation_axes, 'string')
svariables['rotation_order'] = (rotation_order, 'string')
svariables['rotation_array'] = (rotation_array, 'float_array')

svariables['align_pmi_on_axis_flag'] = (align_pmi_on_axis_flag, 'boolean')
svariables['pmi_eigenvector'] = (pmi_eigenvector, 'int')
svariables['alignment_vector_axis'] = (alignment_vector_axis, 'string')
svariables['settle_on_plane'] = (settle_on_plane, 'string')
svariables['plane'] = (plane, 'string')


svariables['fasta_utilities_flag'] = (fasta_utilities_flag, 'boolean')

svariables['fasta_input_option'] = (fasta_input_option, 'string')
svariables['fasta_input_sequence'] = (fasta_input_sequence, 'string')
svariables['fasta_input_file'] = (fasta_input_file, 'string')


svariables['seed'] = (seed, 'int_array')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()

#import pprint; pprint.pprint(variables); exit()

txtQueue = multiprocessing.JoinableQueue()
build_utilities = tool_methods.build_utilities()
build_utilities.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)

#print 'in GUI and txtOutput = ', this_text, '\n'





