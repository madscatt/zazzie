
'''
GUI_MIMIC for build_utilities
'''

import sys

import sassie.build.build_utilities.build_methods as build_methods
import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'

pdb_utilities_flag = False
pdb_utilities_flag = True
input_pdbfile = '../../../developer_files_for_testing/torsion_angle_monte_carlo/hiv1_gag.pdb'
input_pdbfile = 'hiv1_gag.pdb'

renumber_flag = True
renumber_flag = False
renumber_output_filename = 'renumbered.pdb'
renumber_indices_flag = False
renumber_resids_flag = False

first_index = '17'
first_resid = '19'

pdb_constraints_flag = False
constraint_options = 'protein'
constraint_options = 'nucleic'
constraint_options = 'solute'
constraint_options = 'backbone'
constraint_options = 'heavy'
constraint_fields = 'occupancy'
constraint_fields = 'beta'

number_of_constraint_files = '1'
number_of_constraint_files = '2'
constraint_options = 'heavy, backbone'
constraint_fields = 'beta, beta'
constraint_resets = 'True, True'
constraint_filenames = 'constrain_heavy.pdb, constrain_backbone.pdb'

modify_fields_flag = True
modify_fields_output_filename = 'modified_fields.pdb'
number_of_fields_to_modify = '2'
field_selections = 'index < 10, (index > 10 and index < 20)'
field_options = 'beta,segname'
field_values = '20.0,PLS'

translation_rotation_flag = False
translation_rotation_output_filename = 'trans_rot.pdb'
pre_center_flag = True
translation_array = '0.0, 0.0, 0.0'
rotation_flag = True
rotation_type = 'cardinal'
rotation_axes_order = 'xyz'
user_vector_1 = 'None'
#rotation_axes_order = 'xzy'
#rotation_axes_order = 'yxz'
#rotation_axes_order = 'yzx'
#rotation_axes_order = 'zxy'
#rotation_axes_order = 'zyx'
rotation_theta = '90.0, 0.0, 0.0'

align_pmi_on_axis_flag = True
align_pmi_on_axis_flag = False
align_pmi_on_cardinal_axes_flag = False
align_pmi_on_cardinal_axes_flag = True
align_pmi_output_filename = 'aligned_on_pmi.pdb'
''' values: 0, 1, or 2 for I1, I2, or I3 '''
pmi_eigenvector = '2'
''' values: x, y, z '''
alignment_vector_axis = 'z'
user_vector_2 = 'None'
settle_on_surface_flag = True
settle_on_surface_flag = False
''' values: x, y, z '''
surface_plane = 'z' 
invert_along_axis_flag = True
invert_along_axis_flag = False

fasta_utilities_flag = True
fasta_utilities_flag = False
fasta_output_filename = 'sequence_from_fasta.pdb'
fasta_input_option = 'file'
fasta_input_file = 'test_fasta.txt'
#fasta_input_option = 'sequence'
fasta_input_sequence = 'RRAGMPSCYLK'
fasta_moltype = 'protein'

seed = '0, 123'  # set this to '1,123' if you want to set the seed or '0,123' if not

#### end user input ####
#### end user input ####
#### end user input ####


svariables['runname'] = (runname, 'string')
svariables['input_pdbfile'] = (input_pdbfile, 'string')

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

svariables['modify_fields_flag'] = (modify_fields_flag, 'boolean')
svariables['modify_fields_output_filename'] = (modify_fields_output_filename,
                                               'string')
svariables['number_of_fields_to_modify'] = (number_of_fields_to_modify, 'int')
svariables['field_selections'] = (field_selections, 'string')
svariables['field_options'] = (field_options, 'string')
svariables['field_values'] = (field_values, 'string')

svariables['translation_rotation_flag'] = (translation_rotation_flag, 'boolean')
svariables['translation_rotation_output_filename'] = (translation_rotation_output_filename, 'string')
svariables['pre_center_flag'] = (pre_center_flag, 'boolean')
svariables['translation_array'] = (translation_array, 'float_array')
svariables['rotation_flag'] = (rotation_flag, 'boolean')
svariables['rotation_type'] = (rotation_type, 'string')
svariables['user_vector_1'] = (user_vector_1, 'string')
svariables['rotation_axes_order'] = (rotation_axes_order, 'string')
svariables['rotation_theta'] = (rotation_theta, 'float_array')

svariables['align_pmi_output_filename'] = (align_pmi_output_filename, 'string')
svariables['align_pmi_on_cardinal_axes_flag'] = (align_pmi_on_cardinal_axes_flag, 'boolean')
svariables['align_pmi_on_axis_flag'] = (align_pmi_on_axis_flag, 'boolean')
svariables['pmi_eigenvector'] = (pmi_eigenvector, 'int')
svariables['alignment_vector_axis'] = (alignment_vector_axis, 'string')
svariables['user_vector_2'] = (user_vector_2, 'string')
svariables['settle_on_surface_flag'] = (settle_on_surface_flag, 'boolean')
svariables['surface_plane'] = (surface_plane, 'string')
svariables['invert_along_axis_flag'] = (invert_along_axis_flag, 'boolean')

svariables['fasta_utilities_flag'] = (fasta_utilities_flag, 'boolean')
svariables['fasta_output_filename'] = (fasta_output_filename, 'string')
svariables['fasta_input_option'] = (fasta_input_option, 'string')
svariables['fasta_input_sequence'] = (fasta_input_sequence, 'string')
svariables['fasta_input_file'] = (fasta_input_file, 'string')
svariables['fasta_input_option'] = (fasta_input_option, 'string')
svariables['fasta_moltype'] = (fasta_moltype, 'string')

svariables['seed'] = (seed, 'int_array')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()

#import pprint; pprint.pprint(variables); exit()

txtQueue = multiprocessing.JoinableQueue()
build_utilities = build_methods.build_utilities()
build_utilities.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)

#print 'in GUI and txtOutput = ', this_text, '\n'





