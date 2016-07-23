import sys,os
import sassie.util.sasconfig as sasconfig
import sasmol.sasmol as sasmol
import sasmol.sasmath as sasmath
import sassie.simulate.energy.readpsf as readpsf
import sassie.simulate.energy.readparam as readparam
import sassie.simulate.torsion_angle_monte_carlo.group_psf as group_psf
import sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.tamc_utilities.setup_torsion_parameters as setup_torsion_parameters
import numpy


def define_main_pivots(direction):
    '''
    pivots are defined per the input moltype

    two types of pivots are main_pivots and secondary_pivots

    "basis" sets of four atoms definining the dihedral torsions in each pivot

    "outside" defines atoms not in the current residue that are needed for each torsion
    where needed

    "terminals" define which basis torsions are sampled at the ends of the chain

    "post" define the atoms in the residue that are affected by rotation, this
    can be appended with definitions for the group outside of this
    definition
    '''

    main_pivots = {}
    main_pivots["basis"] = [
                            ["N", "CA", "C", "NZ"], 
                            #["N", "CA", "C", "O"], 
                            #["CA", "C", "NZ", "CE"], 
                            ["C", "NZ", "CE", "CD"], 
                            ["NZ", "CE", "CD", "CG"], 
                            ["CE", "CD", "CG", "CB"], 
                            ["CD", "CG", "CB", "CA"], 
                            #["CG", "CB", "CA", "C"] 
                            ["CG", "CB", "CA", "N"] 
                           ]

    ''' "outside" defines atoms outside residue for each basis type '''
    '''
    note: for example the second element is 1 --> second basis [n,ca,c,n]
          and it will require an atom from the next residue which is the 3
          element ... "n"
    '''    
    
    main_pivots["outside"] = {
                               0: ["next", [0,1,2]],
                               1: ["next", [0]],
                               2: [],
                               3: [],
                               4: [],
                               5: []
                             }

    ''' "terminals" define allowable basis move for each terminal case : both angles are available at termini '''
   
    ## @NOTE to ZHL: should not be needed for isopeptide
    main_pivots["terminals"] = {}

    '''
    post defines the atoms in the residue that will be moved for each basis type
    '''

    ## @NOTE to ZHL: only graph traversal is currenlty implemented, and post is not needed for graph traversal
    if direction == 'forward':
        main_pivots["post"] = [] # [0, "sidechain", "forward", "other_forward"]
    else:
        main_pivots["post"] = [] # [1, "sidechain", "backward", "other_backward"]
    
    '''
    main_pivots["graph_traversal_exception"] is a dictionary that defines the cases (residue types) where a pivot is potentially located inside a loop, and the actions that need be taken when a pivot is found inside a loop.
    '''

    if direction == 'forward':
        main_pivots["graph_traversal_exception"] = {}
    else:
        main_pivots["graph_traversal_exception"] = {}

    return main_pivots

def build_basis_strings(residue_main_pivots,residue_main_pivots_outside,residue_main_pivots_list, isopeptide_basis_string):

    '''
    method to create basis strings for each main pivot in each residue
    '''
    cter = isopeptide_basis_string.split('or')[0] ## @NOTE to ZHL: assume contract
    lys = isopeptide_basis_string.split('or')[1] ## @NOTE to ZHL: assume contract

    cter_list = cter.split()
    cter_resid = int(cter_list[cter_list.index('resid')+1])
    cter_segname = cter_list[cter_list.index('segname')+1]
    lys_list = lys.split()
    lys_resid = int(lys_list[lys_list.index('resid')+1])
    lys_segname = lys_list[lys_list.index('segname')+1]
    #print cter_resid,cter_segname
    #print lys_resid,lys_segname

    residue_main_pivot = residue_main_pivots[0]
    residue_main_pivot_outside = residue_main_pivots_outside[0]
    residue_main_pivot_list = residue_main_pivots_list[0]
    #for pivot,pivot_list in zip(residue_main_pivot,residue_main_pivot_list):
    #    print pivot,pivot_list

    residue_bases = []

    #for resid in resids:
    if True: ## @NOTE only  1 primary residue in isopeptide, which is LYS
        resid = lys_resid

        this_residue_pivots = residue_main_pivots[0]
        residue_number = 0

        this_pivot_counter = 0
        for this_pivot in this_residue_pivots:
            this_pivot_number = residue_main_pivots_list[residue_number][this_pivot_counter]

            this_residue_main_pivots_outside = []
            for key, value in residue_main_pivots_outside[residue_number].iteritems():   # iterate over dictionary
                if key == this_pivot_number:
                    this_residue_main_pivots_outside = residue_main_pivots_outside[residue_number][key]
                    break

            i = 0
            residue_basis = ''

            for this_name in this_pivot:
                keyword = False
                if len(this_residue_main_pivots_outside) > 0:
                    if len(this_residue_main_pivots_outside) == 2:  #['previous', [0]] or ['previous', [0, 1]]
                        if len(this_residue_main_pivots_outside[1]) == 1:
                            if i in this_residue_main_pivots_outside[1]:
                                keyword = this_residue_main_pivots_outside[0]
                        elif len(this_residue_main_pivots_outside[1]) >= 2:
                            if i in this_residue_main_pivots_outside[1]:
                                keyword = this_residue_main_pivots_outside[0]

                if i != len(this_pivot):
                    if not keyword:
                        if i != len(this_pivot) - 1:
                            residue_basis += '(segname[i] == "'+lys_segname+'" and resid[i] == '+str(lys_resid)+' and name[i] == "'+this_name+'") or '
                        else:
                            residue_basis += '(segname[i] == "'+lys_segname+'" and resid[i] == '+str(lys_resid)+' and name[i] == "'+this_name+'")'
                    elif keyword == "previous":
                        print 'isopeptide TAMC along previous direction not implemented yet'
                        print 'Quit now'
                        exit(0)
                    elif keyword == "next":
                        residue_basis += '(segname[i] == "'+cter_segname+'" and resid[i] == '+str(cter_resid)+' and name[i] == "'+this_name+'") or '
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
    residue_main_pivots = []
    residue_main_pivots_outside = []
    residue_main_pivots_list = []

    residue_main_pivots.append(main_pivots["basis"])
    residue_main_pivots_outside.append(main_pivots["outside"])
    residue_main_pivots_list.append(range(len(main_pivots['basis'])))
    residue_main_pivots_post = main_pivots["post"]

    #print residue_main_pivots
    #print residue_main_pivots_outside
    #print residue_main_pivots_list

    residue_basis_strings = build_basis_strings(residue_main_pivots,residue_main_pivots_outside,residue_main_pivots_list, pvars.flexible_basis)

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
