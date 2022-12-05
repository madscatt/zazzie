import sys,os
import sassie.util.sasconfig as sasconfig
import sasmol.system as system
import sasmol.linear_algebra as linear_algebra
import sassie.simulate.energy.readpsf as readpsf
import sassie.simulate.energy.readparam as readparam
import sassie.simulate.torsion_angle_monte_carlo.group_psf as group_psf
import numpy
import pprint

import sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.tamc_utilities.graph_traversal as traversal

def assign_atom_types(mol, atoms):

    names = mol.name()
    natoms = mol.natoms()

    i = 0
    charmm_type_list = []
    for this_name in names:
        if this_name == atoms[i][0]:
            charmm_type_list.append(atoms[i][1])
        else:
            message = 'atom name / type misalignment: mol name = '+this_name+', atom name = '+atoms[i][0]
            log.error(message)

        i += 1

    mol.setCharmm_type(charmm_type_list)

    return

def clean_up_dihedral_parameters(pdihedrals):

    new_dihedrals = []
    for parm in pdihedrals:
        atom_list = []
        parameter_list = []
        i = 0
        for item in parm:
            if i < 4:
                atom_list.append(item)
            else:
                parameter_list.append(item)
            i += 1

        new_dihedrals.append([atom_list,parameter_list])

    return new_dihedrals

def identify_dihedral_match(this_dihedral, new_dihedrals):

    found = False

#    print 'this_dihedral = ', this_dihedral
    td0 = this_dihedral[0] ; td1 = this_dihedral[1]
    td2 = this_dihedral[2] ; td3 = this_dihedral[3]
    this_dihedral_backwards = [td3, td2, td1, td0]
    x_dihedral = ['X', td1, td2, 'X']
    x_dihedral_backwards = ['X', td2, td1, 'X'] ;

    for library_dihedral in new_dihedrals:
        if this_dihedral == library_dihedral[0] or this_dihedral_backwards == library_dihedral[0] \
            or x_dihedral == library_dihedral[0] or x_dihedral_backwards == library_dihedral[0]:
            found = True
            return found, library_dihedral[1]

    return found, None

def compare_two_lists(a, b):
    '''
    determine if each element of lists a & b are identical (order and identity)
    '''

    if len(a) != len(b):
        return False
    else:
        i = 0
        for this_a in a:
            if this_a != b[i]:
                return False
            i += 1

    return True

def grab_single_basis_elements(basis_string):

#(resid[i] == 9 and name[i] == "C") or (resid[i] == 10 and name[i] == "N") or (resid[i] == 10 and name[i] == "CA") or (resid[i] == 10 and name[i] == "C")

    new_basis_strings = []
    for char in basis_string:

        if char == '(':
            this_basis = ''
        elif char == ')':
            new_basis_strings.append(this_basis)
        else:
            this_basis += char

    return new_basis_strings


def get_subset_mask_and_indices_string_order(self,basis_filter):
    '''
    method to instead of sasmol get_subset_mask due to order of atoms being
    incorrect in some PDB files.  This method uses the order of the basis_string
    elements instead of the atom ordering
    '''

    index 	    = self.index()
    name 		= self.name()
    loc		    = self.loc()
    resname	    = self.resname()
    chain		= self.chain()
    resid		= self.resid()
    rescode	    = self.rescode()
    occupancy   = self.occupancy()
    beta		= self.beta()
    segname	    = self.segname()
    element	    = self.element()
    charge	    = self.charge()
    moltype	    = self.moltype()

    error = []

    indices = []

    new_basis_filter = grab_single_basis_elements(basis_filter)

    mask_array = numpy.zeros(self.natoms(),numpy.int32)

    for this_basis_filter in new_basis_filter:
        found = False
        for i in range(self.natoms()):
            try:
                if(eval(this_basis_filter)):
                    mask_array[i] = 1
                    #indices.append(index[i])
                    indices.append(i+1)
                    found = True
            except:
                message = 'failed to evaluate basis filter = '+this_basis_filter+' for atom '+str(i)
                return error, mask_array

        if not found:
            message = 'found no atoms using filter selection '+str(new_basis_filter)
            return error, mask_array

    return error, mask_array, indices

def collect_dihedrals_from_psf_inputs(mol, pivot_indices, atoms, dihedrals, new_dihedrals):
    '''
    method to assign dihedral parameters from PSF and parameter files
    '''

    ''' first determine if all pivots exist in PSF file '''

    count = 0

    main_pivots_parameters = {}

    for this_pivot in pivot_indices:

        found = False
        for this_dihedral in dihedrals:
            this_int_dihedral = [int(i) for i in this_dihedral]
            #### OPEN:  left over data type from legacy method ... should change
            if compare_two_lists(this_int_dihedral, this_pivot) or compare_two_lists(this_int_dihedral, this_pivot[::-1]):
               # print 'found one!'
               # print 'this_pivot = ', this_pivot
               # print 'this_dihedral = ', this_int_dihedral
               # print
                found = True

        if not found:
            p = this_pivot
            atom_st = 'atoms = '+str(atoms[p[0]-1])+' '+str(atoms[p[1]-1])+' '+str(atoms[p[2]-1])+' '+str(atoms[p[3]-1])
            message = 'no valid pivot found for pivot definition = '+str(this_pivot)
            message += '\n'+atom_st
            log.error(message)

        this_pivot_type = []
        for atom in this_pivot:
            this_pivot_type.append(mol.charmm_type()[atom - 1])

        found, correct_dihedral = identify_dihedral_match(this_pivot_type, new_dihedrals)

        if found:
            #print 'correct_dihedral['+str(count)+'] = ', correct_dihedral
            main_pivots_parameters[count] = correct_dihedral
            count += 1

        elif not found:

            print('count = ', count)
            print('this_pivot_type = ', this_pivot_type)

            message = 'no valid pivot found for pivot type definition = '+str(this_pivot_type)
            log.error(message)

    return main_pivots_parameters

def setup_main_pivot_force_field_parameters(mol, main_pivots_indices, main_pivots_masks, residue_main_pivots, input_psf_file):
    '''
    method to assign torsion parameters to each pivot
    '''
    #print '>>> assigning force-field parameters to main_pivots'

    bin_path = sasconfig.__bin_path__

    parameter_file = os.path.join(os.path.sep,bin_path,"toppar","par_all27_prot_na.inp")

    pbonds,pangles,pdihedrals,pimpropers,pnonbond = readparam.getparms(parameter_file)

    new_dihedrals = clean_up_dihedral_parameters(pdihedrals)

    segments, atoms, charge, mass, bonds, angles, dihedrals, impropers = readpsf.getpsf(input_psf_file)

    assign_atom_types(mol,atoms)

    main_pivots_parameters = collect_dihedrals_from_psf_inputs(mol, main_pivots_indices, atoms, dihedrals, new_dihedrals)


    return main_pivots_parameters



def calculate_initial_torsion_variables(other_self, main_pivots_masks):

    mol = other_self.mol

    frame = 0

    torsion_angles = []

    for mask in main_pivots_masks:
        error, coor = mol.get_coor_using_mask(frame,mask)
        if len(error) > 0:
            log.error('ERROR: '+str(error))
        else:
            c = coor[0]
            torsion_angles.append(linear_algebra.dihedral_angle(c[0],c[1],c[2],c[3]))
            #print torsion_angles[-1]

    return

def get_post(group_molecule, group_flexible_mask, direction, main_pivots_indices, main_pivots, input_psf_file=None):
    '''
    get the flexible post
    '''

    group_natoms = len(group_flexible_mask)
    group_flexible_natoms = len(numpy.nonzero(group_flexible_mask)[0])
    main_pivots_post_masks=[]

    DEBUG = False
    if input_psf_file==None:
        if (DEBUG): print("flexible post mask generation not implemented for brute force yet!") ## @NOTE to ZHL: not implemented
        exit(0)
    else:
        '''get the bond list from the external psf file'''
        if (DEBUG): print('\n'+'='*100+'\nGetting the bond list from the external psf file...')
        psf_data = group_psf.psf_data()
        group_psf.parse_psf_file(open(input_psf_file).readlines(), psf_data)
        bonds = traversal.get_bond_list(psf_data)
        if (DEBUG): print('bond list:'); pprint.pprint(bonds)

        '''convert atomic indices to addresses''' ## @NOTE to ZHL: hardwired for the pdb file
        for bond in bonds:
            bond[0] -= 1
            bond[1] -= 1
        if (DEBUG): print('reorganized bond list:'); pprint.pprint(bonds)

        '''Build the graph'''
        if (DEBUG): print('\n'+'='*100+'\nBuilding the graph...')
        graph = traversal.setup_graph(bonds)
        if (DEBUG): print('graph:'); pprint.pprint(graph)

        for indices in main_pivots_indices:
            if (DEBUG): print('selected mask',mask)
            flexible_post_mask = numpy.zeros(group_flexible_natoms)
            if direction=="forward":
                pivot_axis = [indices[1]-1, indices[2]-1] ## @NOTE to ZHL: is this always right?
            elif direction=="backward" or direction=="reverse":
                pivot_axis = [indices[2]-1, indices[1]-1] ## @NOTE to ZHL: is this always right?
            else:
                print("ERROR: expected forward/reverse direction, received", direction)
                exit(1)
            if (DEBUG): print('pivot axis: '),pivot_axis
            traversal.traverse_graph_for_pivot(graph, group_molecule,  main_pivots['graph_traversal_exception'], flexible_post_mask, pivot_axis)
            #traversal.traverse_graph_for_pivot(graph, group_molecule,  flexible_post_mask, pivot_axis)
            if (DEBUG): print(('flexible_post_mask after traversal: '),flexible_post_mask)

            post_mask = numpy.ones(group_natoms)
            count=0
            for i in range(group_natoms):
                if group_flexible_mask[i]:
                    post_mask[i] = flexible_post_mask[count]
                    count += 1
            main_pivots_post_masks.append(post_mask)
    if (DEBUG): pprint.pprint([item.nonzero() for item in main_pivots_post_masks])

    return main_pivots_post_masks

def setup_torsion_parameters(other_self, group_molecule, group_number, pvars):
    '''
    method to initialize dihedral energy parameters for the group
    '''

    log = other_self.log
    mvars = other_self.module_variables
    direction = mvars.rotation_direction_array[group_number]
    torsion_parameter_module = pvars.torsion_parameter_module
    
    group_flexible_mask = other_self.group_flexible_masks[group_number]
    group_post_mask = other_self.group_post_masks[group_number]

    group_flexible_molecule = system.Molecule(0)
    error = group_molecule.copy_molecule_using_mask(group_flexible_molecule,group_flexible_mask,0)

    pvars.main_pivots = torsion_parameter_module.define_main_pivots(direction)

    ## @NOTE to ZHL: need to make sure this is the right host for the list
    number_of_groups = mvars.number_of_flexible_regions
    if 'main_pivots_masks' not in vars(pvars):
        pvars.main_pivots_masks = [[]]*number_of_groups
    if 'residue_main_pivots' not in vars(pvars):
        pvars.residue_main_pivots = [[]]*number_of_groups
    if 'main_pivots_post_masks' not in vars(pvars):
        pvars.main_pivots_post_masks = [[]]*number_of_groups
    if 'main_pivots_parameters' not in vars(pvars):
        pvars.main_pivots_parameters = [[]]*number_of_groups

    main_pivots_indices, pvars.main_pivots_masks[group_number], pvars.residue_main_pivots[group_number] = torsion_parameter_module.assign_main_pivots(group_flexible_molecule, pvars)

    # residue_sidechain_indices, pvars.residue_sidechain_masks, residue_other_indices, pvars.residue_other_masks = assign_sidechain_and_other_atoms(group_molecule, pvars.main_pivots) ## @NOTE to ZHL: sidechain pivot not implemented yet

    group_flexible_psf_name = os.path.join(other_self.run_path, "group_flexible_"+str(group_number)+".psf")
    #pvars.main_pivots_post_masks[group_number] = get_post(group_molecule, group_flexible_mask, direction, pvars.main_pivots_masks[group_number], group_flexible_psf_name)
    #pvars.main_pivots_post_masks[group_number] = get_post(group_molecule, group_flexible_mask, direction, main_pivots_indices, pvars.main_pivots, group_flexible_psf_name)
    pvars.main_pivots_post_masks[group_number] = get_post(group_flexible_molecule, group_flexible_mask, direction, main_pivots_indices, pvars.main_pivots, group_flexible_psf_name)

    group_psf_name = os.path.join(other_self.run_path, "group_flexible_"+str(group_number)+".psf")
    pvars.main_pivots_parameters[group_number] = setup_main_pivot_force_field_parameters(group_flexible_molecule, main_pivots_indices, pvars.main_pivots_masks[group_number], pvars.residue_main_pivots[group_number], group_psf_name)

    #calculate_initial_torsion_variables(other_self, pvars.main_pivots_masks)


    #print 'pvars.main_pivots_masks: '
    #print str(pvars.main_pivots_masks)
    #print

    return
