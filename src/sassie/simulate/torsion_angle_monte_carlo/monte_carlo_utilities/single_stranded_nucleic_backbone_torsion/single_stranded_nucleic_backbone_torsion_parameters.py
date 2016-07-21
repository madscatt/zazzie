import sys,os,copy
import sassie.util.sasconfig as sasconfig
import sasmol.sasmol as sasmol
import sasmol.sasmath as sasmath
import sassie.simulate.energy.readpsf as readpsf
import sassie.simulate.energy.readparam as readparam
import sassie.simulate.monte_carlo.group_psf as group_psf
from sassie.simulate.monte_carlo.monte_carlo_utilities.tamc_utilities.general_utilities import *


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
    main_pivots["basis"] = [["O3'", "P", "O5'", "C5'"], ["P", "O5'", "C5'", "C4'"], ["O5'", "C5'", "C4'", "C3'"], ["C5'", "C4'", "C3'", "O3'"], ["C4'", "C3'", "O3'", "P"], ["C3'", "O3'", "P", "O5'"]]  

    '''
     "outside" defines atoms outside residue for each basis type 
    note: for example the second element is 1 --> second basis [n,ca,c,n]
          and it will require an atom from the next residue which is the 3
          element ... "n"
    '''    
    main_pivots["outside"] = { 0: ["previous", [0]], 1: [], 2:[], 3:[], 4: ["next", [3]], 5: ["next", [2,3]] }

    ''' "terminals" define allowable basis move for each terminal case : both angles are available at termini '''
    main_pivots["terminals"] = { "terminal_first": [1, 2, 3, 4, 5], "terminal_last": [0, 1, 2, 3] }

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
