import sys,string,numpy
import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
from sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.tamc_utilities.sample_torsion import sample_torsion 
from sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.tamc_utilities.setup_torsion_parameters import setup_torsion_parameters 

if sasconfig.__level__ == "DEBUG": DEBUG = True

class torsion_variables():
    def __init__(self,parent = None):
        pass
    def get_module(self, this_group_rotation):
        if this_group_rotation == 'protein_backbone_torsion':
            from sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.protein_backbone_torsion import protein_backbone_torsion_parameters as torsion_parameter_module
        elif this_group_rotation == 'single_stranded_nucleic_backbone_torsion':
            from sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.single_stranded_nucleic_backbone_torsion import single_stranded_nucleic_backbone_torsion_parameters as torsion_parameter_module
        elif this_group_rotation == 'isopeptide_bond_torsion':
            from sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.isopeptide_bond_torsion import isopeptide_bond_torsion_parameters as torsion_parameter_module
        self.torsion_parameter_module = torsion_parameter_module


def sample(other_self,group_number):
    '''
    method to choose flexible region and then sample the move
    '''
   
    log = other_self.log
    
    mcvars = other_self.mcvars
    
    mcvars.trial_accepted = False
            
    group_molecule = other_self.group_molecules[group_number]
    group_mask = other_self.group_masks[group_number]
    group_variables = other_self.group_variables[group_number]

    sample_torsion(other_self,group_number)

    return 


def setup(other_self,group_number):

    ''' this method prepares a group for a flexible isopeptide bond region'''

    mvars = other_self.mvars
    log = other_self.log
    pgui = other_self.run_utils.print_gui
    mol = other_self.group_molecules[group_number]
    
    this_group_rotation = mvars.rotation_type_array[group_number]

    if 'pvars' not in vars(mvars):
        pvars = torsion_variables()
    else:
        pvars = mvars.pvars

    log.debug('in setup')

    ''' create an instance of isopeptide bond variables '''

    frame = 0

    flexible_basis = mvars.basis_string_array[group_number]

    log.info('initializing group = '+str(group_number)+' : '+flexible_basis)

    pvars.get_module(this_group_rotation)
    pvars.flexible_basis = flexible_basis
    setup_torsion_parameters(other_self, mol, group_number, pvars)

    mvars.pvars = pvars

    return

