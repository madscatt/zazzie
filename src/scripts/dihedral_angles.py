# $Id: dihedral_angles.py 2907 2015-11-06 19:20:58Z curtisj $
import numpy
import sys
import sasmol.sasmath as sasmath

def measure(coor, indices, an, this_mask, q0, first_last_resid, molecule_type,
            drude=False):

    ind = numpy.nonzero(this_mask * numpy.arange(1, len(this_mask) + 1))[0]
    this_frame_coor = coor[0, :, :]

    lcoor = numpy.take(this_frame_coor[:, :], ind, 0)

    error = []

    if(molecule_type == 'protein'):

        if(an == 'phi'):
            angle = sasmath.dihedral_angle(lcoor[0, :], lcoor[1, :],
                                           lcoor[2, :], lcoor[3, :])
        elif(an == 'psi'):
            if(q0 == first_last_resid[0]):
                angle = sasmath.dihedral_angle(lcoor[0, :], lcoor[1, :],
                                               lcoor[2, :], lcoor[3, :])
            else:
                angle = sasmath.dihedral_angle(lcoor[1, :], lcoor[2, :],
                                               lcoor[3, :], lcoor[4, :])
        else:
            angle = 0.0
            message = 'error (in rotate.measure): %s angle not in phi/psi' % an
            error.append(message)
            assert not error,  error

    elif(molecule_type == 'rna'):

        if(an == 'alpha'):
            angle = sasmath.dihedral_angle(lcoor[0, :], lcoor[1, :],
                                           lcoor[2, :], lcoor[3, :])
        elif(an == 'beta'):
            if(q0 == first_last_resid[0]):
                angle = sasmath.dihedral_angle(lcoor[0, :], lcoor[1, :],
                                               lcoor[2, :], lcoor[3, :])
            else:
                angle = sasmath.dihedral_angle(lcoor[1, :], lcoor[2, :],
                                               lcoor[3, :], lcoor[4, :])
        elif(an == 'gamma'):
            if(q0 == first_last_resid[0]):
                angle = sasmath.dihedral_angle(lcoor[1, :], lcoor[2, :],
                                               lcoor[3, :], lcoor[4, :])
            else:
                angle = sasmath.dihedral_angle(lcoor[2, :], lcoor[3, :],
                                               lcoor[4, :], lcoor[5, :])
        elif(an == 'delta'):
            if(q0 == first_last_resid[0]):
                angle = sasmath.dihedral_angle(lcoor[2, :], lcoor[3, :],
                                               lcoor[4, :], lcoor[5, :])
            else:
                angle = sasmath.dihedral_angle(lcoor[3, :], lcoor[4, :],
                                               lcoor[5, :], lcoor[6, :])
        elif(an == 'epsilon'):
            if(q0 == first_last_resid[0]):
                angle = sasmath.dihedral_angle(lcoor[3, :], lcoor[4, :],
                                               lcoor[5, :], lcoor[6, :])
            else:
                angle = sasmath.dihedral_angle(lcoor[4, :], lcoor[5, :],
                                               lcoor[6, :], lcoor[7, :])
        elif(an == 'eta'):
            if(q0 == first_last_resid[0]):
                angle = sasmath.dihedral_angle(lcoor[4, :], lcoor[5, :],
                                               lcoor[6, :], lcoor[7, :])
            else:
                angle = sasmath.dihedral_angle(lcoor[5, :], lcoor[6, :],
                                               lcoor[7, :], lcoor[8, :])
        else:
            angle = 0.0
            message = ('error (in rotate.measure): %s angle not in '
                       'alpha/beta/gamma/delta/epsilon/eta' % an)
            error.append(message)
            assert not error,  error

    elif(molecule_type == 'dna'):
        if(first_last_resid[0] == first_last_resid[1]):
            print 'WARNING: only given one DNA base for calculating angles'
            print 'WARNING: change input or revise code'
            print 'WARNING: cowardly refusing to calculate angle, returni'
            return 0
        # order of drude atoms differ from charmm27 (not sure about charmm36)
        if drude:
            i = {"O3'-1": 0,  "P": 1,  "O5'": 2,  "C5'": 3,  "chi3": 4,
                 "chi4": 5,"C4'": 6,  "O4'": 7,  "C1'": 8,  "C3'": 9,
                 "O3'": 10,  "P+1": 11,  "O5'+1": 12}
        else:
            i = {"O3'-1": 0,  "P": 1,  "O5'": 2,  "C5'": 3,  "C4'": 4,
                 "O4'": 5,"C1'": 6,  "chi3": 7,  "chi4": 8,  "C3'": 9,
                 "O3'": 10,  "P+1": 11,  "O5'+1": 12}
        if q0 == first_last_resid[0]:
            if len(lcoor) == 12:
                # missing the O3'-1 atom
                first_last_coor = numpy.zeros((13, 3))
                first_last_coor[1:] = lcoor
            elif len(lcoor) == 11:
                # missing the O3'-1 and P atoms
                first_last_coor = numpy.zeros((13, 3))
                first_last_coor[2:] = lcoor
            else:
                # not sure
                assert False, ('ERROR: cannot calculate dihedral angles from '
                               '%d atoms, not clear which atoms are which) \n'
                               'Please review atoms in PDB and input resids')
            if(an == 'alpha'):
                angle = numpy.nan
                return angle
            elif(an == 'beta' and len(lcoor) == 11):
                angle = numpy.nan
                return angle
            else:
                lcoor = numpy.copy(first_last_coor)

        if(an == 'alpha'):
            angle = sasmath.dihedral_angle(
                lcoor[i["O3'-1"], :], lcoor[i["P"], :],
                lcoor[i["O5'"], :], lcoor[i["C5'"], :])
        elif(an == 'beta'):
            angle = sasmath.dihedral_angle(
                lcoor[i["P"], :], lcoor[i["O5'"], :],
                lcoor[i["C5'"], :], lcoor[i["C4'"], :])
        elif(an == 'gamma'):
            angle = sasmath.dihedral_angle(
                lcoor[i["O5'"], :], lcoor[i["C5'"], :],
                lcoor[i["C4'"], :], lcoor[i["C3'"], :])
        elif(an == 'delta'):
            angle = sasmath.dihedral_angle(
                lcoor[i["C5'"], :], lcoor[i["C4'"], :],
                lcoor[i["C3'"], :], lcoor[i["O3'"], :])
        elif(an == 'epsilon'):
            if(q0 == first_last_resid[1]):
                assert len(lcoor) == 11, ('ERROR: expected 11 atoms for last '
                                          'resid but received %d' % len(lcoor))
                angle = numpy.nan
            else:
                angle = sasmath.dihedral_angle(
                    lcoor[i["C4'"], :], lcoor[i["C3'"], :],
                    lcoor[i["O3'"], :], lcoor[i["P+1"], :])
        elif(an == 'zeta'):
            if(q0 == first_last_resid[1]):
                assert len(lcoor) == 11, ('ERROR: expected 11 atoms for last '
                                          'resid but received %d' % len(lcoor))
                angle = numpy.nan
            else:
                angle = sasmath.dihedral_angle(
                    lcoor[i["C3'"], :], lcoor[i["O3'"], :],
                    lcoor[i["P+1"], :], lcoor[i["O5'+1"], :])
        elif(an == 'chi'):
            angle = sasmath.dihedral_angle(
                lcoor[i["O4'"], :], lcoor[i["C1'"], :],
                lcoor[i["chi3"], :], lcoor[i["chi4"], :])
        else:
            angle = numpy.nan
            message = ('error (in rotate.measure): %s angle not in '
                       'alpha/beta/gamma/delta/epsilon/zeta/chi' % an)
            error.append(message)
            assert not error,  error

    return angle


def get_rotation_indices(m1, molecule_type, flexible_residues, txtOutput):

    print 'getting rotation indices for molecule'

    residue_rotation_indices = {}
    residue_rotation_mask = {}

    if(molecule_type == 'protein'):

        mtype = 0
        mask = m1.get_dihedral_subset_mask(flexible_residues, mtype)

        for i in xrange(len(mask)):
            this_mask = mask[i][:]
            q0 = flexible_residues[i]
            residue_rotation_mask[q0] = this_mask.tolist()
            indices = m1.get_indices_from_mask(this_mask)
            residue_rotation_indices[q0] = indices.tolist()

        # print 'residue_rotation_indices = ', residue_rotation_indices
        # print 'residue_rotation_mask = ', residue_rotation_mask

    elif(molecule_type == 'rna'):

        mtype = 1
        mask = m1.get_dihedral_subset_mask(flexible_residues, mtype)

        for i in xrange(len(mask)):
            this_mask = mask[i][:]
            q0 = flexible_residues[i]
            residue_rotation_mask[q0] = this_mask.tolist()
            indices = m1.get_indices_from_mask(this_mask)
            residue_rotation_indices[q0] = indices.tolist()

            # print 'residue_rotation_indices = ', residue_rotation_indices
            # print 'residue_rotation_mask = ', residue_rotation_mask

    elif(molecule_type == 'dna'):

        mtype = 2
        mask = m1.get_dihedral_subset_mask(flexible_residues, mtype)

        for i in xrange(len(mask)):
            this_mask = mask[i][:]
            q0 = flexible_residues[i]
            residue_rotation_mask[q0] = this_mask.tolist()
            indices = m1.get_indices_from_mask(this_mask)
            residue_rotation_indices[q0] = indices.tolist()

            # print 'residue_rotation_indices = ', residue_rotation_indices
            # print 'residue_rotation_mask = ', residue_rotation_mask

    else:
        message = ('rotation basis set not defined for molecule type  = %s' %
                   molecule_type)
        print_failure(message, txtOutput)

    print 'done getting rotation indices'
    sys.stdout.flush()

    return residue_rotation_indices, residue_rotation_mask


def get_flexible_residues(numranges, reslow, numcont):
    '''
    Method to determine residues that are going to be rotated.
    '''
    flexible_residues = []
    templist = []
    for i in xrange(numranges):
        thisreslow = reslow[i]
        thisnumcont = numcont[i]
        thisrange = numpy.arange(thisnumcont) + thisreslow
        templist.append(thisrange.tolist())

    flexible_residues = [item for sublist in templist for item in sublist]

    print 'flexible_resiudes = ', flexible_residues

    return flexible_residues
