import sys,os,copy
import sassie.util.sasconfig as sasconfig
import sasmol.system as system
import sassie.simulate.energy.readpsf as readpsf
import sassie.simulate.energy.readparam as readparam
import sassie.simulate.torsion_angle_monte_carlo.group_psf as group_psf
import sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.tamc_utilities.setup_torsion_parameters as setup_torsion_parameters
import numpy
import pprint

def log_error(message):

    try:
        raise ValueError(message)
    except(Exception, err):
        sys.stderr.write('\n\n')
        sys.stderr.write('ERROR: '+str(err))
        sys.stderr.write('\n\nQUITTING NOW\n\n')
        sys.exit()

def get_residue_sequence(mol):

    '''
    NOTE: numpy SORTS the array when unique is called : so out of order
    sequences are going to be wrong
    '''
    #### OPEN:  need to deal with this either by contract or fix ordering somehow

    resid = mol.resid()
    resname = mol.resname()
    resids, unique_indices = numpy.unique(resid, return_index=True)
    resnames = numpy.take(resname,unique_indices)

    return resids,resnames

def get_unique_atoms(list):
    '''
    utility method to get unique elements from a nested list
    '''
    unique_list = []
    for i in range(len(list)):
        this_i = list[i]
        for j in range(len(list[i])):
            this_j = this_i[j]
            if this_j not in unique_list:
                unique_list.append(this_j)

    return unique_list

def check_if_all_atoms_captured(mol,unique_basis_atoms,residue_sidechain_atoms):

    frame = 0
    resids,resnames = get_residue_sequence(mol)
    other_atoms_forward, other_atoms_backward = get_other_atoms(mol)

    i = 0

    missing_dict = {}

    for this_resid in resids:
        error,mask = mol.get_subset_mask('resid[i] == '+str(this_resid))
        temp_mol = system.Molecule(0)
        mol.copy_molecule_using_mask(temp_mol,mask,frame)
        this_resname = temp_mol.resname()[0]
        this_names = temp_mol.name()

        test_names = []
        for name in unique_basis_atoms:
            test_names.append(name)

        for name in residue_sidechain_atoms[i]:
            test_names.append(name)

        for name in other_atoms_forward[this_resname]:
            test_names.append(name)

        for name in other_atoms_backward[this_resname]:
            test_names.append(name)

        missing_names = []
        for name in this_names:
            if name not in test_names:
                missing_names.append(name)

        if len(missing_names) > 0:
            missing_dict[this_resid] = missing_names

        i += 1

    return missing_dict

def build_basis_strings(resids,residue_main_pivots,residue_main_pivots_outside,residue_main_pivots_list):

    '''
    method to create basis strings for each main pivot in each residue
    '''

    residue_number = 0
    residue_bases = []

    for resid in resids:

        this_residue_pivots = residue_main_pivots[residue_number]

        this_pivot_counter = 0
        for this_pivot in this_residue_pivots:
            this_pivot_number = residue_main_pivots_list[residue_number][this_pivot_counter]

            this_residue_main_pivots_outside = []

            #for k, v in d.items() ##python 3
            for key, value in residue_main_pivots_outside[residue_number].items():   # iterate over dictionary
                if key == this_pivot_number:
                    this_residue_main_pivots_outside = residue_main_pivots_outside[residue_number][key]
                    break

            i = 0
            residue_basis = ''

            for this_name in this_pivot:
                keyword = False
                keyword_array = False
                if len(this_residue_main_pivots_outside) > 0:
                    if len(this_residue_main_pivots_outside) == 2:  #['previous', [0]] or ['previous', [0, 1]]
                        if len(this_residue_main_pivots_outside[1]) == 1:
                            if i in this_residue_main_pivots_outside[1]:
                                keyword = this_residue_main_pivots_outside[0]
                        elif len(this_residue_main_pivots_outside[1]) == 2:
                            if i in this_residue_main_pivots_outside[1]:
                                keyword = this_residue_main_pivots_outside[0]
                        elif len(this_residue_main_pivots_outside[1]) > 2:

                            message = 'program written to handle two atoms outside residue, you defined \
                                        : '+str(len(this_residue_main_pivots_outside))+ ' atoms not in resid: '\
                                        + str(resid) + '\n'

                            log_error(message)

                if i != len(this_pivot):
                    if not keyword:
                        if i != len(this_pivot) - 1:
                            residue_basis += '(resid[i] == '+str(resid)+' and name[i] == "'+this_name+'") or '
                        else:
                            residue_basis += '(resid[i] == '+str(resid)+' and name[i] == "'+this_name+'")'
                    elif keyword == "previous":
                        residue_basis += '(resid[i] == '+str(resid-1)+' and name[i] == "'+this_name+'") or '
                    elif keyword == "next":
                        residue_basis += '(resid[i] == '+str(resid+1)+' and name[i] == "'+this_name+'") or '
                i += 1

            ''' clean up last "or" statement at end of basis string '''

            if residue_basis[-3:-1] == 'or':
               residue_basis = residue_basis[:-3]

            residue_bases.append(residue_basis)

            this_pivot_counter += 1

        residue_number += 1

    return residue_bases

def assign_main_pivots(group_flexible_molecule, pvars):
    ''' method cycles through structure and determines index number for each residue '''


    main_pivots = pvars.main_pivots
    resids,resnames = get_residue_sequence(group_flexible_molecule)
    first_resid = resids[0] ; last_resid = resids[-1]
    residue_main_pivots = []
    residue_main_pivots_outside = []

    residue_main_pivots_list = []

    for resid in resids:

### OPEN : RE-DO THIS TO USE TERMINAL INFORMATION EXPLICITLY NOT BY CALLING OTHER DICTIONARY TO GET ELEMENTS
### OPEN : RE-DO THIS TO COMBINE INTO A GENERAL METHOD WITHOUT MOLTYPE DESIGNATION

        #main_pivots["outside"] = { 0: ["previous", [0]], 1: ["next", [3]] }

        #main_pivots["terminals"] = { "terminal_first": [1], "terminal_last": [0] }


        if resid == first_resid:
            residue_main_pivots.append([main_pivots["basis"][i] for i in main_pivots["terminals"]["terminal_first"]])
            residue_main_pivots_list.append(main_pivots["terminals"]["terminal_first"][:])
            #residue_main_pivots_outside.append({1: main_pivots["outside"][main_pivots["terminals"]["terminal_first"][0]]})
            tmp_dict = copy.copy(main_pivots["outside"])
            for idx in range(len(main_pivots["basis"])):
                if idx not in main_pivots["terminals"]["terminal_first"]:
                    del tmp_dict[idx]
            residue_main_pivots_outside.append(tmp_dict)
        elif resid == last_resid:
            residue_main_pivots.append([main_pivots["basis"][i] for i in main_pivots["terminals"]["terminal_last"]])
            residue_main_pivots_list.append(main_pivots["terminals"]["terminal_last"][:])
            #residue_main_pivots_outside.append({0: main_pivots["outside"][main_pivots["terminals"]["terminal_last"][0]]})
            tmp_dict = copy.copy(main_pivots["outside"])
            for idx in range(len(main_pivots["basis"])):
                if idx not in main_pivots["terminals"]["terminal_last"]:
                    del tmp_dict[idx]
            residue_main_pivots_outside.append(tmp_dict)
        else:
            residue_main_pivots.append(main_pivots["basis"])
            residue_main_pivots_outside.append(main_pivots["outside"])
            residue_main_pivots_list.append(range(len(main_pivots["basis"])))


        #main_pivots["outside"] = { 0: ["previous", [0,1]], 1: ["previous", [0]], \
        #                           4: ["next", [3]], 5: ["next", [2,3]] }

        #main_pivots["terminals"] = { "terminal_first": [2, 3, 4, 5], "terminal_last": [0, 1, 2, 3] }


        residue_main_pivots_post = main_pivots["post"]

    residue_basis_strings = build_basis_strings(resids,residue_main_pivots,residue_main_pivots_outside,residue_main_pivots_list)

    main_pivots_indices = []
    main_pivots_masks = []

    #print 'len(resids) = ',len(resids)
    #print 'len(residue_basis_strings) = ', len(residue_basis_strings)

    #print 'residue_basis_strings[-1] = ',residue_basis_strings[-1]

    i = 1
    for basis in residue_basis_strings:
        error, this_mask, this_indices = setup_torsion_parameters.get_subset_mask_and_indices_string_order(group_flexible_molecule,basis)
        main_pivots_indices.append(this_indices)
        #main_pivots_indices.append(numpy.nonzero(this_mask)[0]+1) ## @NOTE to ZHL: temporary fix for a bug
        main_pivots_masks.append(this_mask)
        #print 'this_indices ['+str(i)+'] = ', this_indices
        i += 1

    #print 'main_pivots_indices[-1] = ', main_pivots_indices[-1]
    #print 'main_pivots_masks[-1] = ', main_pivots_masks[-1]

    return main_pivots_indices, main_pivots_masks, residue_main_pivots
