# Hailiang Zhang
# March 2015

import os,sys,locale
import pprint
import copy

import sasmol.sasmol as sasmol

import sassie.simulate.torsion_angle_monte_carlo.group_psf as group_psf

def get_bond_list(psf_data):
    """
    Functionality: get bond list from psf data
    Input: psf_data -- the psf data obtained from group_psf.parse_psf_file
    Return: bonds -- a list of bond pairs
    """
    bonds = []
    for line in psf_data.nbond:
        words = line.split()
        number_of_bonds_this_line = len(words)/2
        if not number_of_bonds_this_line:
            continue
        if words[1]=='!NBOND:':
            number_of_bonds = locale.atoi(words[0])
            continue
        for i in xrange(number_of_bonds_this_line):
            bonds.append([locale.atoi(words[2*i]), locale.atoi(words[2*i+1])])
    return bonds

def setup_graph(bonds):
    """
    Functionality: set up the graph from a bond list
    Input: bonds -- a list of bond pairs (get from get_bond_list)
    Return: graph -- the graph data structure as 'list of successors'
    """
    graph = {}
    for bond in bonds:
        a = bond[0]
        b = bond[1]
        if a==b:
            print('self-looping is not allowed in my graph')
            print('quitting program')
            exit(0)
        if a not in graph:
            graph[a] = []
        graph[a].append(b)
        if b not in graph:
            graph[b] = []
        graph[b].append(a)
    return graph

def traverse_self_util(graph, mol, resid, mask, current, back_watch):
    """
    Functionality: the recursive function for graph traversal
    Input: graph -- the graph data structure as 'list of successors' (get from setup_graph)
           mask -- a mask for the to-be-rotated atoms
           current -- current
           back_watch -- the node that should be watched for situation of pivot-inside-a-loop
    """
    if mol.resid()[current] != resid:
        return
    if mask[current]:
        return
    mask[current] = 1
    for child in graph[current]:
        if child == back_watch:
            return 
        traverse_self_util(graph, mol, resid, mask, child, back_watch)

def traverse_try_util(graph, mol, exceptions, mask, current, back_watch):
    """
    Functionality: the recursive function for graph traversal
    Input: graph -- the graph data structure as 'list of successors' (get from setup_graph)
           mask -- a mask for the to-be-rotated atoms
           current -- current
           back_watch -- the node that should be watched for situation of pivot-inside-a-loop
    """
    if mask[current]:
        return True
    mask[current] = 1
    for child in graph[current]:
        if child == back_watch:
            return False
        return traverse_recursive_watch(graph, mol, exceptions, mask, child, back_watch)

def traverse_recursive_restore(graph, mol, mask, current, back_watch):
    """
    Functionality: the recursive function for graph traversal
    Input: graph -- the graph data structure as 'list of successors' (get from setup_graph)
           mask -- a mask for the to-be-rotated atoms
           current -- current
    """
    print '^^^\ncurrent in traverse_recursive_restore: (',mol.resid()[current],',',mol.resname()[current],',',mol.name()[current],')'
    if mask[current] in [0,1] or current==back_watch:
        print 'return\nvvv'
        return
    elif mask[current]>1:
        print 'reset current'
        mask[current] = 0;

    for child in graph[current]:
        print 'child in traverse_recursive_restore: (',mol.resid()[child],',',mol.resname()[child],',',mol.name()[child],')'
        traverse_recursive_restore(graph, mol, mask, child, back_watch)

    print 'return\nvvv'

    return

def traverse_recursive_try(graph, mol, mask, current, back_watch):
    """
    Functionality: the recursive function for graph traversal
    Input: graph -- the graph data structure as 'list of successors' (get from setup_graph)
           mask -- a mask for the to-be-rotated atoms
           current -- current
    """
    if mask[current]:
        return 
    mask[current] = 2;

    for child in graph[current]:
        if child == back_watch:
            raise Exception
        traverse_recursive_try(graph, mol, mask, child, back_watch)

def traverse_recursive_watch(graph, mol, exceptions, mask, current, back_watch):
    """
    Functionality: the recursive function for graph traversal
    Input: graph -- the graph data structure as 'list of successors' (get from setup_graph)
           mask -- a mask for the to-be-rotated atoms
           current -- current
           back_watch -- the node that should be watched for situation of pivot-inside-a-loop
    """
    if mask[current]:
        return
    mask[current] = 1 

    resname = mol.resname()[current]
    resid = mol.resid()[current]
    name = mol.name()[current]

    flag_exception = False
    if resname in exceptions and name==exceptions[resname][0]:
        flag_exception = True

    for child in graph[current]:

        if flag_exception:
            if not (mol.resid()[child]==resid and mol.name()[child] in exceptions[resname][1]):
                mask_tmp = copy.copy(mask)
                try:
                    traverse_recursive_try(graph, mol, mask_tmp, child, back_watch)
                except:
                    continue
        if child == back_watch:
            #print "privot found inside a loop without exceptions"
            #print vars()
            continue
        
        traverse_recursive_watch(graph, mol, exceptions, mask, child, back_watch)

    return

def traverse_graph_for_pivot(graph, mol, exceptions, mask, pivot, direction='forward'):
    """
    Functionality: graph traversal function customized for TAMC
    Input: graph -- the graph data structure as 'list of successors' (get from setup_graph)
           mask -- a mask for the to-be-rotated atoms
           pivot -- the pivot (a bond pair) that is to be rotated
           direction -- rotating direction ('forward' or 'backword')
    Return:
    main_pivots["graph_traversal_exception"] = {"PRO": ("CA",["N","C"]), "CYS": (["SG",["CB"])}
    """
    """
    '''reset the mask'''
    '''
    for i in xrange(len(mask)):
        mask[i] = 0
    '''
    '''set the pivot atoms'''
    if direction == 'forward':
        back_watch = pivot[0]
        root = pivot[1]
    elif direction == 'backward':
        back_watch = pivot[1]
        root = pivot[0]
    else:
        print('rotating direction needs either to be "forward" or "backward"')
        print('quitting program')
        exit(0)
    '''traverse the graph'''
    mask[root] = 1
    for child in graph[root]:
        if child == back_watch:
            continue
        traverse(graph, mask, child, back_watch)
    """
    '''traverse the graph'''
    mask[pivot[1]] = 1
    for child in graph[pivot[1]]:
        if child == pivot[0]:
            continue
        traverse_recursive_watch(graph, mol, exceptions, mask, child, pivot[0])


'''my test'''
if __name__=='__main__':

    '''set up the sasmol object'''
    print('\n'+'='*100+'\nSetting up the sasmol object...')
    pdb_file = 'ten_mer.pdb'
    m = sasmol.SasMol(0)
    m.read_pdb(pdb_file)
    natoms = m.natoms()
    mask = [0]*natoms

    '''get the bond list from the external psf file'''
    print('\n'+'='*100+'\nGetting the bond list from the external psf file...')
    psf_file = 'ten_mer.psf'
    psf_data = group_psf.psf_data()
    group_psf.parse_psf_file(open(psf_file).readlines(), psf_data)
    bonds = get_bond_list(psf_data)
    #print('bond list:'); pprint.pprint(bonds)

    '''convert atomic indices to addresses''' ## @NOTE to ZHL: hardwired for the pdb file
    for bond in bonds:
        bond[0] -= 1
        bond[1] -= 1
    #print('reorganized bond list:'); pprint.pprint(bonds)

    '''Build the graph'''
    print('\n'+'='*100+'\nBuilding the graph...')
    graph = setup_graph(bonds)
    #print('graph:'); pprint.pprint(graph)

    '''pick up a pivot and rotating direction and traverse the graph'''
    print('\n'+'='*100+'\nTesting...')
    #rotate N-CA bond of residue 1 in forward direction (first residue forward)
    pivot = bonds[0]
    direction = 'forward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of N-CA bond of residue 1 in forward direction (first residue forward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)

    #rotate N-CA bond of residue 1 in forward direction (first residue backward)
    pivot = bonds[0]
    direction = 'backward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of N-CA bond of residue 1 in forward direction (first residue backward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)

    #rotate N-CA bond of residue 3 in backward direction (middle residue forward)
    pivot = [19,21]
    direction = 'forward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of N-CA bond of residue 3 in backward direction (middle residue forward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)

    #rotate N-CA bond of residue 3 in backward direction (middle residue backward)
    pivot = [19,21]
    direction = 'backward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of N-CA bond of residue 3 in backward direction (middle residue backward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)

    #rotate CB-CG bond of residue 5 (Arg) in forward direction (sidechain rotation forward)
    pivot = [48,51]
    direction = 'forward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of CB-CG bond of residue 5 (Arg) in forward direction (sidechain rotation forward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)

    #rotate CB-CG bond of residue 5 (Arg) in backward direction (sidechain rotation backward)
    pivot = [48,51]
    direction = 'backward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of CB-CG bond of residue 5 (Arg) in backward direction (sidechain rotation backward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)

    #rotate CB-CG bond of last residue in backward direction (last residue forward)
    pivot = [135,130]
    direction = 'forward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of CB-CG bond of last residue in backward direction (last residue forward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)

    #rotate CB-CG bond of last residue in backward direction (last residue backward)
    pivot = [135,130]
    direction = 'backward'
    traverse_graph_for_pivot(graph, mask, pivot, direction)
    print('\n'+'.'*100+'\nTesting for rotation of CB-CG bond of last residue in backward direction (last residue backward)...')
    print('mask of the to-be-rotated atoms for pivot (%i, %i) in %s direction'%(pivot[0],pivot[1],direction)); print(mask)
