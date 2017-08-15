#       check overlap between segments
#
#       06/26/09        --      initial coding			:       sr/jc
#       11/18/11        --      adapted to sasmol		:	jc
#
import os
import sys
import string
import numpy
#import sassie.simulate.monte_carlo.complex.interres as interres


def moloverlap(coords1, coords2, cutoff, interres=[]):
    atomlist = []
    for i in range(len(coords1[0])):
        x1 = coords1[0, i, 0]
        y1 = coords1[0, i, 1]
        z1 = coords1[0, i, 2]
        for j in range(len(coords2[0])):
            x2 = coords2[0, j, 0]
            y2 = coords2[0, j, 1]
            z2 = coords2[0, j, 2]
            qex = 0
            for k in range(len(interres)):
                if i == interres[k][0] and j == interres[k][1]:
                    qex = 1
            if qex:
                continue
            distance = numpy.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            if distance <= cutoff:
                atomlist.append([i, j])
    '''
	print '\n\n\nZHL ',interres
	print 'cutoff = ', cutoff
	print 'number of overlap = ',len(atomlist),atomlist
	'''
    return atomlist


def nmer_overlap_check(m1, path, pdbfile, cutoff, basis, keyword_basis):

    segname = m1.segname()
    natoms = m1.natoms()

    segs = []
    for tseg in segname:
        if tseg not in segs:
            segs.append(tseg)

    basisall = []
    basisindex = []
    n = len(segs)
    pair = []
    for i in range(n):
        for j in range(i, n):
            seg1 = segs[i]
            seg2 = segs[j]
            if(seg1 != seg2):
                pair.append([segs[i], segs[j]])
                basisindex.append([i, j])
    npairs = len(pair)

    if not keyword_basis:

        for i in range(npairs):
            x = basis[basisindex[i][0]].strip()
            y = basis[basisindex[i][1]].strip()
            basisall.append([x, y])

    inter = []
    total = []

    print '>>> CALCULATING INITIAL OVERLAP ARRAYS'

    for i in range(npairs):
        if not keyword_basis:
            a = str(basisall[i][0])
            b = str(basisall[i][1])
            mol1_basis = 'name[i] == "' + a + \
                '" and segname[i] == "' + pair[i][0] + '"'
            mol2_basis = 'name[i] == "' + b + \
                '" and segname[i] == "' + pair[i][1] + '"'
        else:
            # TODO --> all not attempted yet
            if(basis[0].lower() == 'all'):
                print 'setting up all atom overlap arrays'
            #	m1.set_average_vdw()
            #	npairs = m1.natoms()*(m1.natoms() - 1)/2
            #	cutoff_array = numpy.zeros(npairs,numpy.float)
            #	pairs.pairs(m1.atom_vdw(),cutoff_array)
            #	basis_filter = 'not name[i] == None'
            elif(basis[0].lower() == 'backbone' or basis[0].lower() == 'heavy'):
                mol1_basis = '(not name[i][0] == "H") and segname[i] == "' + \
                    pair[i][0] + '"'
                mol2_basis = '(not name[i][0] == "H") and segname[i] == "' + \
                    pair[i][1] + '"'
            else:
                print 'ERROR IN ASSIGNING BASIS IN NMER_OVERLAP_CHECK'
                sys.exit()

        error, mol1_mask = m1.get_subset_mask(mol1_basis)
        error, coords1 = m1.get_coor_using_mask(0, mol1_mask)
        error, mol2_mask = m1.get_subset_mask(mol2_basis)
        error, coords2 = m1.get_coor_using_mask(0, mol2_mask)

        # TODO --> the following routine needs to be in fortran
        # HERE here
        na1 = len(coords1[0])
        na2 = len(coords2[0])
        np = (na1 + na2) * ((na1 + na2) - 1) / 2
        np = na1 * na2
        atom_list = numpy.zeros(np, numpy.int)
        atomlist = moloverlap(coords1, coords2, cutoff)
        interres.mover(coords1[0], coords2[0], cutoff, atom_list)
        if(len(atomlist) > 0):
            # inter.append([i,atomlist])
            inter.append([basisindex[i], atomlist])
        total.append(len(atomlist))

    print 'inter = ', inter
    print 'total = ', total
    print 'total overlap for a cutoff of ', cutoff, ' angstroms is ', sum(total)
    print 'FORTRAN: total overlap for a cutoff of ', cutoff, ' angstroms is ', numpy.sum(atom_list)

    return inter, sum(total)
