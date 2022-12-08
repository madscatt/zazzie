'''
Driver method to run the torsion_angle_monte_carlo module
'''

import sys

import sassie.simulate.torsion_angle_monte_carlo.monte_carlo as torsion_angle_monte_carlo
import sassie.interface.input_filter as input_filter
import multiprocessing
import time

def user_variables(self, **kwargs):

    #### user input ####
    #### user input ####
    #### user input ####

    self.run_name = 'run_0'
    self.dcdfile = 'run_0.dcd'
    self.pdbfile = '../../../../src_2.7/developer_files_for_testing/torsion_angle_monte_carlo/ten_mer.pdb'
    self.psffile = '../../../../src_2.7/developer_files_for_testing/torsion_angle_monte_carlo/ten_mer.psf'
    self.psf_flag = True # sh: making this false to prevent crash
    #self.psf_flag = False # sh: making this false to prevent crash

    self.max_steps = '5000'
    self.energy_convergence = '1.0'
    self.step_size = '0.002'

    self.number_of_flexible_regions = '1'

    self.basis_string_array = []
    self.basis_string_array.append('resid >= 2 and resid < 5')
    self.delta_theta_array = '30.0'
    self.rotation_type_array = ['protein_backbone_torsion']
    self.rotation_direction_array = ['forward']
    self.post_basis_string_array = ['resid>=5 and resid<=10']
    #self.rotation_direction_array = ['backward']
    #self.post_basis_string_array = ['resid<3']
    self.overlap_basis = 'heavy'
    
    self.temperature = '300.0'
    self.trial_steps = '100'
    self.goback = '1'
    
    self.low_rg_cutoff = '0'
    self.high_rg_cutoff = '400.0'

    self.z_flag = False
    self.z_cutoff = '0.0'
    
    self.constraint_flag = False
    self.constraint_file = 'constraints.txt'

    self.directed_mc = '0'

    self.nonbondflag = '0'
    self.seed = '0, 123'  # set this to '1,123' if you want to set the seed or '0,123' if not


    #### end user input ####
    #### end user input ####
    #### end user input ####

def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['run_name'] = (self.run_name, 'string')
    svariables['dcdfile'] = (self.dcdfile, 'string')
    svariables['pdbfile'] = (self.pdbfile, 'string')
    svariables['psffile'] = (self.psffile, 'string')
    svariables['psf_flag'] = (self.psf_flag, 'boolean')

    svariables['max_steps'] = (self.max_steps, 'int')
    svariables['energy_convergence'] = (self.energy_convergence, 'float')
    svariables['step_size'] = (self.step_size, 'float')

    svariables['number_of_flexible_regions'] = (self.number_of_flexible_regions, 'int')
    svariables['basis_string_array'] = (self.basis_string_array, 'string')
    svariables['delta_theta_array'] = (self.delta_theta_array, 'float_array')
    svariables['rotation_type_array'] = (self.rotation_type_array, 'string')
    svariables['rotation_direction_array'] = (self.rotation_direction_array, 'string')
    svariables['overlap_basis'] = (self.overlap_basis, 'string')
    svariables['post_basis_string_array'] = (self.post_basis_string_array, 'string')
    svariables['temperature'] = (self.temperature, 'float')
    svariables['trial_steps'] = (self.trial_steps, 'int')
    svariables['goback'] = (self.goback, 'int')
    svariables['directed_mc'] = (self.directed_mc, 'float')


    svariables['low_rg_cutoff'] = (self.low_rg_cutoff, 'float')
    svariables['high_rg_cutoff'] = (self.high_rg_cutoff, 'float')

    svariables['z_flag'] = (self.z_flag, 'boolean')
    svariables['z_cutoff'] = (self.z_cutoff, 'float')

    svariables['constraint_flag'] = (self.constraint_flag, 'boolean')
    svariables['constraint_file'] = (self.constraint_file, 'string')

    svariables['nonbondflag'] = (self.nonbondflag, 'int')
    svariables['seed'] = (self.seed, 'int_array')

    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print('error = ', error)
        sys.exit()

    import pprint; pprint.pprint(self.variables); #exit()

    txtQueue = multiprocessing.JoinableQueue()
    simulation = torsion_angle_monte_carlo.simulation()
    simulation.main(self.variables, txtQueue)
    this_text = txtQueue.get(True, timeout=0.1)

    #print 'in GUI and txtOutput = ', this_text, '\n'

class gui_mimic_torsion_angle_monte_carlo():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'torsion_angle_monte_carlo'

    def __init__(self, test, paths):

        if not test:
            user_variables(self)
        else:
            test_variables(self, paths)

        run_module(self)

if __name__ == '__main__':

    test = False  
    paths = None

# We are thinking of defining the install path so the gui mimic can be run from anywhere as long as it is called from that particular python
# That way, the test files will always be available to the user.
    if test:
        pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
        dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
        module_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'interface', 'align') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_torsion_angle_monte_carlo(test, paths)
    print("time used: ", time.time() - start)


