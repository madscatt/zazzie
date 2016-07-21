#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Align several structures using an align basis
# Created: 6 February 2015
#
# $Id: align_driver.py 3038 2016-03-01 19:55:50Z schowell $
#
'''
This driver serves as a wrapper tool for the align method. This is useful
for aligning a pdb/dcd to the coordinates of a goal pdb structure.
'''

import logging
import sys
import os
import sasmol


class inputs():
    '''
    create an empty inputs holder
    '''
    def __init__(self, parent=None):
        pass


def parse():
    '''
    Returns arguments after parsing command line input
    '''

    parser = argparse.ArgumentParser(
        description='for aligning a pdb/dcd to a goal pdb structure'
    )

    parser.add_argument("-g", "--goal", help="goal pdb")
    parser.add_argument("-r", "--ref",
                        help="pdb with atom info for pdb/dcd that is moving")
    parser.add_argument("-m", "--move", help="pdb/dcd to be moved/aligned")
    parser.add_argument("-o", "--out", help="output dcd file name")
    parser.add_argument("-ms", "--move_seg_chain",
                        help="segname or chain to match")
    parser.add_argument("-gs", "--goal_seg_chain",
                        help="segname or chain to match")
    parser.add_argument("-msc", "--move_seg_or_chain",
                        help="matching a segname or chain")
    parser.add_argument("-gsc", "--goal_seg_or_chain",
                        help="matching a segname or chain")
    parser.add_argument("-p", "--path", help="output path")
    parser.add_argument("-gmn", "--goal_min",
                        help="minimun residue to match on goal molecule")
    parser.add_argument("-gmx", "--goal_max",
                        help="maximum residue to match on goal molecule")
    parser.add_argument("-mmn", "--move_min",
                        help="minimun residue to match on move molecule")
    parser.add_argument("-mmx", "--move_max",
                        help="miximum residue to match on move molecule")
    parser.add_argument("-ba", "--basis_atoms", help="basis_atoms to match")

    return parser.parse_args()


def align(inputs):
    '''
    purpose:
    --------
        method for aligning a pdb/dcd file to another pdb file

    input:
    ------
        inputs: object should contain the following attributes
            goal:    goal pdb file name
            ref:     pdb file with molecule info for the pdb/dcd to be moved
            move:    pdb/dcd file to align
            out:     output dcd file
            path:    output path
            goal_filter:     goal basis filter
            move_filter:     move basis filter

    alternative input:
    ------------------
        inputs:      instead of `inputs.goal_filter` and `inputs.move filter`
            basis_atoms:    atom names to use for alignment, e.g., 'CA'
            goal_seg_or_ch: must be set to either 'segname' or 'chain'
            goal_segname:   segname or chain to match on goal molecule
            goal_res_max:   maximum residue to match on goal molecule
            goal_res_min:   minimum residue to match on goal molecule
            move_seg_or_ch: must be set to either 'segname' or 'chain'
            move_segname:   segname or chain to match on move molecule
            move_res_max:   maximum residue to match on move molecule
            move_res_min:   minimum residue to match on move molecule

    returns:
    --------
        aligned sasmol object

    note:
    -----
        inputs.ref and inputs.move are ofter the same pdb
    '''
    aa_goal_pdb = inputs.goal
    aa_move_pdb = inputs.ref
    aa_move_file = inputs.move
    save_file = inputs.out
    path = inputs.path

    # get the goal_filter
    try:
        goal_filter = inputs.goal_filter
    except:
        basis_atoms = inputs.basis_atoms
        goal_seg_or_ch = inputs.goal_seg_or_chain
        goal_segname = inputs.goal_seg_chain
        goal_res_max = inputs.goal_max
        goal_res_min = inputs.goal_min
        goal_filter = ('((%s[i] == "%s") and (name[i] == "%s") and '
                       '(resid[i] >= %s) and (resid[i] <= %s))' % (
                           goal_seg_or_ch, goal_segname, basis_atoms,
                           goal_res_min, goal_res_max))

    # get the move_filter
    try:
        move_filter = inputs.move_filter
    except:
        basis_atoms = inputs.basis_atoms
        move_seg_or_ch = inputs.move_seg_or_chain
        move_segname = inputs.move_seg_chain
        move_res_max = inputs.move_max
        move_res_min = inputs.move_min
        move_filter = ('((%s[i] == "%s") and (name[i] == "%s") and '
                       '(resid[i] >= %s) and (resid[i] <= %s))' % (
                           move_seg_or_ch, move_segname, basis_atoms,
                           move_res_min, move_res_max))

    # check input files exist
    assert os.path.exists(aa_move_file), ('ERROR: no such file - %s' %
                                          aa_move_file)
    assert os.path.exists(aa_move_pdb), ('ERROR: no such file - %s' %
                                         aa_move_pdb)
    assert os.path.exists(aa_goal_pdb), ('ERROR: no such file - %s' %
                                         aa_goal_pdb)

    # create all the SasMol objects
    sub_goal = sasmol.sasmol.SasMol(0)
    sub_move = sasmol.sasmol.SasMol(0)
    aa_goal = sasmol.sasmol.SasMol(0)
    aa_move = sasmol.sasmol.SasMol(0)

    # populate the all-atom SasMol objects from the pdb files
    aa_goal.read_pdb(aa_goal_pdb)
    aa_move.read_pdb(aa_move_pdb)

    # identify file type to align and the number of frames
    if aa_move_file[-3:] == 'pdb':
        aa_move.read_pdb(aa_move_file)
        n_frames = aa_move.number_of_frames()
        in_type = 'pdb'
    elif aa_move_file[-3:] == 'dcd':
        dcd_file = aa_move.open_dcd_read(aa_move_file)
        n_frames = dcd_file[2]
        in_type = 'dcd'
    else:
        message = "\n~~~ ERROR, unknown input type ~~~\n"
        print_failure(message, txtOutput)
        return

    # setup the output
    out_type = save_file[-3:].lower()
    if 'dcd' == out_type:
        dcd_out_file = aa_move.open_dcd_write(path + save_file)
    elif 'pdb' == out_type:
        dcd_out_file = None
    else:
        raise NotImplementedError

    # get the subset masks
    error, goal_seg_mask = aa_goal.get_subset_mask(goal_filter)
    assert not error, error
    error, move_seg_mask = aa_move.get_subset_mask(move_filter)
    assert not error, error

    # populate the subset SasMol objects
    error = aa_goal.copy_molecule_using_mask(sub_goal, goal_seg_mask, 0)
    assert not error, error
    error = aa_move.copy_molecule_using_mask(sub_move, move_seg_mask, 0)
    assert not error, error

    # calculate the center of mass of the subset of m1
    com_sub_goal = sub_goal.calccom(0)

    # center the m1 coordinates
    sub_goal.center(0)

    # store the m1 centered coordinates
    coor_sub_goal = sub_goal.coor()[0]

    for i in xrange(n_frames):
        if in_type == 'dcd':
            # read in the coordinates for the frame to be moved
            aa_move.read_dcd_step(dcd_file, i)

            # move m2 to be centered at the origin
            aa_move.center(0)
            error, sub_move.coor = aa_move.get_coor_using_mask(0, move_seg_mask)
            assert not error, error
            sub_move.setCoor(sub_move.coor)

            # calculate the center of mass of the subset of m2
            com_sub_move = sub_move.calccom(0)

            # move the subset of m2 to be centered at the origin
            sub_move.center(0)

            # store the new coordinates of the subset of m2
            coor_sub_move = sub_move.coor[0]

            # align m2 using the transformation from sub_m2 to sub_m1
            aa_move.align(0, coor_sub_move, com_sub_move, coor_sub_goal,
                          com_sub_goal)

        elif in_type == 'pdb':
            # move m2 to be centered at the origin
            aa_move.center(i)
            error, sub_move.coor = aa_move.get_coor_using_mask(i, move_seg_mask)
            sub_move.setCoor(sub_move.coor)

            # calculate the center of mass of the subset of m2
            com_sub_move = sub_move.calccom(0)

            # move the subset of m2 to be centered at the origin
            sub_move.center(0)

            # store the new coordinates of the subset of m2
            coor_sub_move = sub_move.coor[0]

            # align m2 using the transformation from sub_m2 to sub_m1
            aa_move.align(i, coor_sub_move, com_sub_move, coor_sub_goal,
                          com_sub_goal)

        # store the result for this frame
        if out_type == 'dcd':
            aa_move.write_dcd_step(dcd_out_file, 0, i + 1)
        elif out_type == 'pdb':
            raise NotImplementedError
        else:
            raise NotImplementedError

    # close file pointers
    if in_type == 'dcd':
        aa_move.close_dcd_read(dcd_file[0])

    if out_type == 'dcd':
        aa_move.close_dcd_write(dcd_out_file)

    print 'alignment complete\n\m/ >.< \m/'


def align_mol(inputs):
    '''
    purpose:
    --------
        method for aligning two sasmol objects

    input:
    ------
        intputs: object that should contain the following attributes
            aa_goal:    goal sasmol object
            aa_move:    sasmol object to align
            goal_basis: goal basis for alignment
            move_basis: move basis for alignment

    returns:
    --------
        aligned sasmol object

    note:
    -----
        inputs.ref and inputs.move are ofter the same pdb
    '''
    aa_goal = inputs.aa_goal
    aa_move = inputs.aa_move
    goal_basis = inputs.goal_basis
    move_basis = inputs.move_basis

    # create the subset SasMol objects for the aligning
    sub_goal = sasmol.sasmol.SasMol(0)
    sub_move = sasmol.sasmol.SasMol(0)

    # get the subset masks
    error, goal_seg_mask = aa_goal.get_subset_mask(goal_basis)
    assert not error, error
    error, move_seg_mask = aa_move.get_subset_mask(move_basis)
    assert not error, error

    # populate the subset SasMol objects
    error = aa_goal.copy_molecule_using_mask(sub_goal, goal_seg_mask, 0)
    assert not error, error
    error = aa_move.copy_molecule_using_mask(sub_move, move_seg_mask, 0)
    assert not error, error

    # calculate the center of mass of the subset of m1
    com_sub_goal = sub_goal.calccom(0)

    # center the m1 coordinates
    sub_goal.center(0)

    # store the m1 centered coordinates
    coor_sub_goal = sub_goal.coor()[0]

    # move m2 to be centered at the origin
    aa_move.center(0)
    error, sub_move.coor = aa_move.get_coor_using_mask(0, move_seg_mask)
    assert not error, error
    sub_move.setCoor(sub_move.coor)

    # calculate the center of mass of the subset of m2
    com_sub_move = sub_move.calccom(0)

    # move the subset of m2 to be centered at the origin
    sub_move.center(0)

    # store the new coordinates of the subset of m2
    coor_sub_move = sub_move.coor[0]

    # align m2 using the transformation from sub_m2 to sub_m1
    aa_move.align(0, coor_sub_move, com_sub_move, coor_sub_goal, com_sub_goal)

    print 'alignment complete/n\m/ >.< \m/'


if __name__ == "__main__":

    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' % __name__, level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    # make ARGS global
    ARGS = parse()

    align(ARGS)
