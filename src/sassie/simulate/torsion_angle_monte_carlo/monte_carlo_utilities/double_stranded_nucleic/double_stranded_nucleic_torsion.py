# $Id: double_stranded_nucleic_torsion.py 3046 2016-03-04 16:35:58Z schowell $
import numpy
import time
import inspect
import warnings
import os
import sys

import sasmol.system as system
    
import sassie.util.sasconfig as sasconfig
import sassie.util.basis_to_python as basis_to_python
import sassie.simulate.torsion_angle_monte_carlo.dna_overlap as dna_overlap

# ??? THIS NEEDS TO BE DEBUGGED AND INTEGRATED INTO THE MAIN SASSIE BRANCH ???
# try:
    # # cpp_sasmol specific imports:
    # # place_of_sasmol_cpp = ('/home/schowell/myPrograms/sassie/sassie_2.0/'
        # # 'trunk/sassie/simulate/monte_carlo/monte_carlo_utilities/'
        # # 'double_stranded_nucleic/')
    # # sys.path.append(place_of_sasmol_cpp)
    # from sasmol_as_cpp import create_sasmol_as_cpp
    # from get_bead_masks import get_bead_masks_in_cpp
    # use_sasmol_cpp = True
    # print 'succeeded in loading sasmol_cpp'
# except:
    # use_sasmol_cpp = False
    # print 'failed to load sasmol_cpp'
use_sasmol_cpp = False

if sasconfig.__level__ == "DEBUG": DEBUG = True

class double_stranded_nucleic_variables():
    def __init__(self, parent = None):
        pass


def get_chain_parameters(coordinates):
    """
    method to generate the u-vectors pointing from bead to bead
    returns the u-vectors and their magnitudes
    """

    bead_vectors = coordinates[1:, :]-coordinates[:-1, :] # unit vec btwn beads
    bead_distances = numpy.sqrt(bead_vectors[:, 0]**2 +
                                bead_vectors[:, 1]**2 +
                                bead_vectors[:, 2]**2) # magnitudes of u-vecs
    n_vectors = bead_distances.size
    bead_unit_vectors = bead_vectors/bead_distances.reshape(n_vectors,1)

    return bead_unit_vectors, bead_distances


def get_twist_angles(xyz_vectors):
    """
    1. orient z2 along z1
    2. calculate the angle between x2 and x1
    (should be identical to the angle btwn y2 and y1)

    Parameters
    ----------
    xyz_vectors : Nx3 np.array
        The bead axes

    Returns
    -------
    twist_angles : N-1 x 3 np.array
        The twist angles between dsDNA bps
    """

    xyz_1 = xyz_vectors[:,:-1]
    xyz_2 = xyz_vectors[:,1:]
    xyz_2_aligned = numpy.zeros(xyz_2.shape)

    # get the rotation matrices to rotate Z2 onto Z1
    R = rotate_v2_to_v1(xyz_1[2], xyz_2[2])

    twist_angles = numpy.zeros(len(R))
    for i in xrange(len(R)):
        # align the z-axes
        for j in xrange(3):
            xyz_2_aligned[j,i] = numpy.dot(R[i], xyz_2[j,i])

        # get angle between x-axes, always be between 0 and 180
        angle_i_rad = numpy.arccos(numpy.dot(xyz_2_aligned[0, i], xyz_1[0, i]))
        angle_i = angle_i_rad * 180.0/numpy.pi

        # check if the angle is positive/negative
        cross_term = numpy.cross(xyz_1[0, i], xyz_2_aligned[0, i])
        pos_neg = numpy.dot(xyz_1[2, i], cross_term)
        if pos_neg < 0:
            angle_i *= -1.0

        twist_angles[i] = angle_i

    return twist_angles


def get_twist_energy(xyz_vectors):
    """
    this function determines the twist energy between a chain of dna basepairs
    (Uses Wilma Olson's 3DNA energy for twisting DNA)

    Parameters
    ----------
    xyz_vectors : Nx3 np.array
        The bead axes

    Returns
    -------
    twist_energy : float
        The twist energy for ds_dna
    """

    twist_angles = get_twist_angles(xyz_vectors)

    force_constant = 0.062
    # average_angle = 34.175
    average_angle = 36  # pg 1217 of doi:10.1038/nprot.2008.104
    twist_energy = numpy.sum(0.5 * force_constant *
                             (twist_angles-average_angle)**2)

    return twist_energy


def rotate_v2_to_v1(v1, v2):
    """
    get the rotation matrix that rotate from v2 to v1 along the vector
    orthongal to both

    Parameters
    ----------
    v1 : Nx3 np.array
        The vector/s to rotate to
    v2 : Nx3 np.array
        The vector/s to rotate from


    Returns
    -------
    R : Nx3x3 np.array
        The rotation matrices that rotate from v2 to v1
    """

    if len(v2.shape) != len(v1.shape):
        if len(v2.shape) > len(v1.shape):
            new_v1 = np.zeros(v2.shape)
            new_v1[:] = v1
            v1 = new_v1
        else:
            pass # THIS WILL BREAK

    if len(v2.shape) == 1:
        # make sure v1 and v2 are unit vectors:
        v1 = v1/numpy.sqrt(v1.dot(v1))
        v2 = v2/numpy.sqrt(v2.dot(v2))

        v = numpy.cross(v2, v1) # hinge axis
        c = numpy.dot(v2, v1)
        R = numpy.eye(3)
        s = numpy.sqrt(numpy.dot(v, v))
        if s != 0:
            V = numpy.array([[    0, -v[2],  v[1]],
                             [ v[2],     0, -v[0]],
                             [-v[1],  v[0],    0]])
            R += V + numpy.dot(V, V) * (1 - c) / s**2
    else:
        # make sure v1 and v2 are unit vectors:
        v1 = numpy.einsum('ij,i->ij', v1, 1/numpy.sqrt(
            numpy.einsum('ij,ij->i', v1, v1)))
        v2 = numpy.einsum('ij,i->ij', v2, 1/numpy.sqrt(
            numpy.einsum('ij,ij->i', v2, v2)))

        v = numpy.cross(v2, v1)
        c = numpy.einsum('ij,ij->i', v2, v1)
        s = numpy.sqrt(numpy.einsum('ij,ij->i', v, v))

        n_v = len(v)

        # V = numpy.array([[    0, -v[2],  v[1]],
                         # [ v[2],     0, -v[0]],
                         # [-v[1],  v[0],    0]])
        V = numpy.zeros((n_v, 3, 3))
        V[:,0,1:]  = v[:,2:0:-1]
        V[:,1,::2] = v[:,::-2]
        V[:,2,:2]  = v[:,1::-1]
        V[:,0,1] = -V[:,0,1]
        V[:,1,2] = -V[:,1,2]
        V[:,2,0] = -V[:,2,0]

        I = numpy.zeros((n_v, 3, 3))
        I[:] = numpy.eye(3)

        R = I + V + numpy.einsum('ijk,i->ijk', numpy.einsum(
            'hij,hjk->hik', V, V), (1 - c) / s**2)

    return R


def get_bend_energy(persistence_length, bead_unit_vectors, bead_distances):
    """
    this function finds the E/kT for N beads connected by N-1 rods
    lpl is the persistence length divided by the rod length: lp/l
    u contains the N-1 unit vectors representing each rod
    """

    vec_k0 = bead_unit_vectors[:-1, :]
    vec_k1 = bead_unit_vectors[1:, :]
    result = numpy.sum((1-numpy.inner(vec_k0, vec_k1).diagonal()
                        )/bead_distances[:-1])

    return result * persistence_length


def get_wca_energy(w, coor, wca0, trial_bead):
    """
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective of the chain
    coor contains the xyz-coordinates of the beads
    and does all this using FORTRAN
    """

    wca1 = numpy.copy(wca0)
    # calculate the distance between moved DNA beads and fixed DNA beads
    # use trial_bead+1 because FORTRAN indexes from 1 instead of 0
    wca1 = dna_overlap.wca_d(coor, trial_bead+1, w, wca1)

    res = 4.*numpy.sum(wca1)

    return res, wca1


def get_dna_bp_reference_frame(other_self, dna_ids, bp_mol, bp_masks=None,
                               dna_id_type='segname'):
    """
    The x-axis points in the direction of the major groove along what would
    be the pseudo-dyad axis of an ideal Watson-Crick base-pair, i.e. the
    perpendicular bisector of the C1'...C1' vector spanning the base-pair.
    The y-axis runs along the long axis of the idealizosed base-pair in the
    direction of the sequence strand, parallel with the C1'...C1' vector,
    and displaced so as to pass through the intersection on the
    (pseudo-dyad) x-axis of the vector connecting the pyrimidine Y(C6) and
    purine R(C8) atoms. The z-axis is defined by the right-handed rule,
    i.e. z = x cross y. (doi:10.1006/jmbi.2001.4987)
    """
    log = other_self.log

    dna_resnames1 = ['CYT', 'THY', 'GUA', 'ADE']
    dna_resnames2 = ['DC', 'DT', 'DG', 'DA']
    if any([name in bp_mol.resnames() for name in dna_resnames1]):
        dna_resnames = dna_resnames1
    elif any([name in bp_mol.resnames() for name in dna_resnames2]):
        dna_resnames = dna_resnames2
    else:
        log.error('ERROR, failed to find dna resnames:\n{}\n{}'.format(
            dna_resnames1, dna_resnames2))


    # c6c8_string = ('(((resname[i] == "CYT" or resname[i] == "THY") and '
                   # 'name[i] == "C6") or '
                   # '((resname[i] == "GUA" or resname[i] == "ADE") and '
                   # 'name[i] == "C8"))')

    c6c8_string = ('(((resname[i] == "{}" or resname[i] == "{}") and '
                   'name[i] == "C6") or '
                   '((resname[i] == "{}" or resname[i] == "{}") and '
                   'name[i] == "C8"))'.format(dna_resnames[0], dna_resnames[1],
                                              dna_resnames[2], dna_resnames[3]))

    dna1_c1p_filter  = '%s[i] == "%s" and name[i] == "C1\'" ' % (
        dna_id_type.lower(), dna_ids[0])
    dna2_c1p_filter  = '%s[i] == "%s" and name[i] == "C1\'" ' % (
        dna_id_type.lower(), dna_ids[1])
    dna1_c6c8_filter = '%s[i] == "%s" and %s' % (
        dna_id_type.lower(), dna_ids[0], c6c8_string)
    dna2_c6c8_filter = '%s[i] == "%s" and %s' % (
        dna_id_type.lower(), dna_ids[1], c6c8_string)

    if not bp_masks:
        e0, dna1_c1p_mask = bp_mol.get_subset_mask(dna1_c1p_filter)
        e1, dna2_c1p_mask = bp_mol.get_subset_mask(dna2_c1p_filter)
        e2, dna1_c6c8_mask = bp_mol.get_subset_mask(dna1_c6c8_filter)
        e3, dna2_c6c8_mask = bp_mol.get_subset_mask(dna2_c6c8_filter)
        assert (4 == dna1_c1p_mask.sum() + dna1_c6c8_mask.sum() +
                dna2_c1p_mask.sum() + dna2_c6c8_mask.sum()), ("ERROR: input "
                    "molecule lacks atoms needed for determining orientation")
        bp_masks = (dna1_c1p_mask, dna2_c1p_mask, dna1_c6c8_mask,
                    dna2_c6c8_mask)
    else:
        dna1_c1p_mask, dna2_c1p_mask, dna1_c6c8_mask, dna2_c6c8_mask = bp_masks

    dna1_c1p = numpy.dot(dna1_c1p_mask, bp_mol.coor()[0])
    dna2_c1p = numpy.dot(dna2_c1p_mask, bp_mol.coor()[0])
    dna1_c6c8 = numpy.dot(dna1_c6c8_mask, bp_mol.coor()[0])
    dna2_c6c8 = numpy.dot(dna2_c6c8_mask, bp_mol.coor()[0])

    y_vec = dna1_c1p - dna2_c1p
    y_mag = numpy.sqrt(numpy.dot(y_vec, y_vec))
    y_hat = y_vec/y_mag

    # following: http://geomalgorithms.com/a05-_intersect-1.html
    Q0 = (dna1_c1p + dna2_c1p)/2
    P0 = dna2_c6c8
    P1 = dna1_c6c8
    w  = P0 - Q0
    u  = P1 - P0
    s1 = numpy.dot(-y_hat, w)/numpy.dot(y_hat, u)
    assert 0 <= s1 <= 1, "ERROR: problem in calculating bead origin"
    bp_origin = P0 + s1 * u

    a = bp_origin - dna2_c1p
    x_vec = a - numpy.dot(a, y_hat) * y_hat
    x_hat = x_vec/numpy.sqrt(numpy.dot(x_vec, x_vec))

    z_hat = numpy.cross(x_hat, y_hat)

    bp_axes = numpy.array([x_hat, y_hat, z_hat])

    return bp_origin, bp_axes, bp_masks


def rotate_coarse_grained_beads(other_self, group_number):
    """
    method to choose a cg bead then rotate
    """

    log = other_self.log

    mvars = other_self.mvars
    nvars = other_self.mvars.nvars
    run_path = other_self.runpath
    mcvars = other_self.mcvars
    group_mol = other_self.group_molecules[group_number]
    group_mask = other_self.group_masks[group_number]
    flexible_mol = other_self.group_flexible_molecules[group_number]
    flexible_mask = other_self.group_flexible_masks[group_number]
    post_mol = other_self.group_post_molecules[group_number]
    post_mask = other_self.group_post_masks[group_number]

    n_flex_regions = mvars.number_of_flexible_regions
    theta_max = mvars.delta_theta_array[group_number]

    cg_flex_mol = nvars.coarse_grained_flexible_mol[group_number]
    xyz_vectors = nvars.xyz_vectors[group_number]
    all_beads = nvars.all_beads[group_number]
    all_overlap_mask = nvars.all_overlap_mask
    post_overlap_mask = nvars.post_overlap_mask[group_number]
    bead_masks = nvars.nucleic_acid_bead_masks[group_number]
    number_beads = nvars.number_beads[group_number]
    dna_segnames = nvars.dna_segnames[group_number]

    # ~~Refresh the group coordinates~~
    # error, group_coor = full_mol.get_coor_using_mask(0, group_mask)
    # group_mol.setCoor(group_coor)
    error, post_coor  = group_mol.get_coor_using_mask(0, post_mask)
    post_mol.setCoor(post_coor)
    error, flexible_coor = group_mol.get_coor_using_mask(0, flexible_mask)
    flexible_mol.setCoor(flexible_coor)

    # considered having this as an advanced input it depends on buffer salt
    # (users could input based depending on sample conditons)
    try:
        persistence_length = nvars.persistence_length
    except:
        # 10 mM Na^+ value (10.1126/science.8079175, 10.1126/science.1439819)
        persistence_length = nvars.persistence_length = 534.0 # 534.0 A

    try:
        trials = nvars.trials
    except:
        trials = nvars.trials = 1

    n_cg_beads = cg_flex_mol.natoms()

    # refresh the coordinates of the cg_beads and their orientations,
    # if not already done, save the bp masks used for calculating twist energy
    cg_coor = cg_flex_mol.coor()
    bp_masks = nvars.bp_masks[group_number]
    save_bp_masks = False
    if not bp_masks:
        save_bp_masks = True
        bp_masks = [None] * number_beads

    for j in xrange(number_beads):
        # update the bead coordinates from the all-atom flexible molecule
        bead = all_beads[j]
        error, bead_coor = flexible_mol.get_coor_using_mask(0, bead_masks[j])
        bead.setCoor(bead_coor)

        # get the bead origin and orientation vectors, store in cg_flex_mol
        (cg_coor[0,j,:], xyz_vectors[:,j,:], bp_masks[j]
         ) = get_dna_bp_reference_frame(other_self, dna_segnames, bead,
                                        bp_masks[j])

        # get the coordinates of the bp relative to the bp_origin and axes
        bead.setCoor(numpy.dot(bead.coor()-cg_coor[0,j,:],
                               xyz_vectors[:,j,:].transpose()))

        # store the bead coordinates for reverse coarse-graining
        all_beads.append(bead)

    if save_bp_masks:
        nvars.bp_masks[group_number] = bp_masks

    cg_flex_mol.setCoor(cg_coor)

    # initialize variables for each run
    xyz = numpy.copy(xyz_vectors)
    flex_coor = numpy.copy(cg_flex_mol.coor()[0]) # unique memory for each
    post_coor = numpy.copy(post_mol.coor()[0]) # unique memory for each
    # steps_from_0 = numpy.zeros(steps, dtype='int64')

    (bead_unit_vectors, bead_distances) = get_chain_parameters(flex_coor)

    rod_length = bead_distances.mean()
    dna_cutoff = 19.0 # dna atomic width measured using VMD: 19 A
    dna_bead_skip = int(numpy.ceil(dna_cutoff/rod_length))
    post_width = 1.0  # 1 A
    overlap_distance = (rod_length + post_width) / 2.0

    try:
        nvars.chain_width[group_number]
    except:
        nvars.chain_width = [None] * n_flex_regions
    chain_width = nvars.chain_width[group_number]
    if chain_width is None:
        chain_width = rod_length / 2.0**(1.0/6.0)
        nvars.chain_width[group_number] = chain_width

    # calculate the energy of the starting positions
    wca_energy_matrix_old = numpy.zeros((n_cg_beads, n_cg_beads))
    bend_energy_old = get_bend_energy(persistence_length, bead_unit_vectors,
                                      bead_distances)

    (wca_energy_old, wca_energy_matrix_old) = get_wca_energy(chain_width,
                                        flex_coor, wca_energy_matrix_old, 0)

    twist_energy_old = get_twist_energy(xyz)

    total_energy_old = bend_energy_old + wca_energy_old + twist_energy_old

    n_accept   = 0 # total times configuration was accepted
    n_reject   = 0 # total times configuration was rejected
    n_written  = 0 # total times dcd write has been called
    n_reload = [0] # list containing the i_goback values

    while n_accept < trials:
        # Choose a bead to rotate (excludes the first and last beads)
        # trial_bead = numpy.random.randint(1, n_cg_beads-1)
        trial_bead = numpy.random.randint(2, n_cg_beads-1)

        # Determine rotation to perform (angles in degrees)
        theta_max_xyz = numpy.array([theta_max]*3)
        theta_xyz = numpy.zeros(3)
        theta_xyz[numpy.random.randint(3)] = 1.0 # select rotation axis, w/ z
        # theta_xyz[numpy.random.randint(2)] = 1.0 # select rotation axis, w/o z
        # theta_xyz[2] = 1.0 # select rotation axis, always z

        # select angle between -theta_max_xyz and +theta_max_xyz
        theta_xyz *= theta_max_xyz * (2 * numpy.random.random() - 1)
        log.debug('[x, y, z] rotation angle: [%f, %f, %f]' %
                  (theta_xyz[0], theta_xyz[1], theta_xyz[2]) )

        # generate a newly rotated model
        (flex_coor[trial_bead:], xyz[:, trial_bead:], post_coor) = bead_rotate(
            other_self, flex_coor[trial_bead-1:], xyz[:, trial_bead-1:],
            theta_xyz, post_coor)

        # calculate the change in energy (dU) and the boltzman factor (p)
        (bead_unit_vectors, bead_distances) = get_chain_parameters(flex_coor)
        bend_energy_new = get_bend_energy(persistence_length,
                                          bead_unit_vectors, bead_distances)

        # ~~~~ DNA interaction energies  ~~~~~~#
        (wca_energy_new, wca_energy_matrix_new) = get_wca_energy(chain_width,
                                flex_coor, wca_energy_matrix_old, trial_bead)
        twist_energy_new = get_twist_energy(xyz)

        total_energy_new =  bend_energy_new + wca_energy_new + twist_energy_new
        energy_change = total_energy_new - total_energy_old

        probability = get_boltzman_factor(energy_change, log)

        test = numpy.random.random()
        overlap = 0
        if test > probability:  # this means test <= probability will pass
            ds_nucleic_acid_pass = False
            log.debug(' ~~~~ trial step failed: dsDNA bend/twist energy ~~~~')
        else:
            ds_nucleic_acid_pass = True

            #~~ check the group for overlap ~~#

            # verify distant strands of cgDNA are not too close
            # (skip closest cg_beads)
            if dna_bead_overlap(other_self, flex_coor, dna_bead_skip,
                                dna_cutoff):
                overlap += 1
                log.debug(' ~~~~ trial step failed: CG bead overlap ~~~~')

            # verify post atoms do not overlap stationary DNA beads
            # stationary NA:      flex_coor[:trial_bead]
            # post (heavy atoms): post_heavy_coor
            if post_mol.natoms() > 0 and not overlap:
                error, post_heavy_coor = post_mol.get_coor_using_mask(
                    0, post_overlap_mask)
                handle_error(other_self, error)
                post_heavy_coor = post_heavy_coor[0, :, :]
                if two_body_overlap(other_self, flex_coor[:trial_bead],
                                    post_heavy_coor, overlap_distance):
                    overlap += 1
                    log.debug(' ~~~~ trial step failed: post overlap ~~~~')

            else:
                post_heavy_coor = numpy.array((0,3))

        if ds_nucleic_acid_pass and not overlap:
            mcvars.trial_accepted = True
            n_accept += 1

            xyz_vectors = numpy.copy(xyz)    # update DNA bp orientations
            cg_flex_mol.setCoor(numpy.array([flex_coor])) # update flex coor
            post_mol.setCoor(numpy.array([post_coor]))    # update post coor

            # update dsDNA WCA energy and total energy
            wca_energy_matrix_old = numpy.copy(wca_energy_matrix_new)
            total_energy_old = total_energy_new

            # ~~recover the aa-DNA representation~~
            recover_all_atom_dsNA(other_self, cg_flex_mol, flexible_mol,
                                  xyz_vectors, all_beads, bead_masks)

            # ~~Combine aa Complete Structure~~
            group_mol.set_coor_using_mask(post_mol, 0, post_mask)
            group_mol.set_coor_using_mask(flexible_mol, 0, flexible_mask)

        else:
            n_reject += 1
            if (mcvars.fail_tally + n_reject == mvars.goback and
                mvars.goback > 1):
                mcvars.trial_accepted = False
                mcvars.number_overlap_failures += 1
                mcvars.fail_tally = mvars.goback - 1 # monte_carlo will add 1

                log.debug('too many steps rejected, signalling goback')
                break
            else:
                # reset group coordinates for another trial
                flex_coor = numpy.copy(cg_flex_mol.coor()[0])
                post_coor = numpy.copy(post_mol.coor()[0])
                xyz = numpy.copy(xyz_vectors)

    n_trials = n_accept + n_reject
    log.debug('exiting ds_nucleic, accepted %d/%d trials (%0.2f %%)' %
              (n_accept, n_trials, 100.0 * float(n_accept) / n_trials) )


def post_process(other_self, group_number):
    """
    track the number of accepted steps
    """
    mvars = other_self.mvars
    nvars = mvars.nvars
    mcvars = other_self.mcvars

    this_step = mcvars.accepted

    # increment following group counters (b/c some steps may not be dsDNA)
    nvars.steps_per_bead[this_step:, group_number] += (1.0 /
                nvars.number_beads[group_number])

    # # openMM destroys DNA
    # log = other_self.log
    # try:
        # min_frequency = nvars.min_frequency
    # except:
        # min_frequency = nvars.min_frequency = 0.1
    # psf_flag = mvars.psf_flag
    # pgui = other_self.run_utils.print_gui
    # # ~~perform minimization~~
    # if nvars.steps_per_bead[this_step, group_number] >= min_frequency:
        # if psf_flag:
            # log.info('minimizing group %d' % group_number)
            # flexible_mol = other_self.group_flexible_molecules[group_number]
            # flexible_mol.write_pdb('before_min.pdb', 0, 'w')
            # group_flexible_simulation = other_self.group_flexible_simulation[
                # group_number]

            # tic = time.time()
            # open_mm.energy_minimization(group_flexible_simulation,
                                        # flexible_mol, mvars.max_steps,
                                        # mvars.energy_convergence)
            # toc = time.time() - tic

            # flexible_mol.write_pdb('after_min.pdb', 0, 'w')

            # flexible_mask = other_self.group_flexible_masks[group_number]
            # group_mol = other_self.group_molecules[group_number]
            # group_mol.set_coor_using_mask(flexible_mol, 0, flexible_mask)

            # if DEBUG:
                # pgui('Minimization for group %d (%d atoms) took %0.2f s' %
                           # (group_number, flexible_mol.natoms(), toc))
            # log.debug('minimizing group %d took %0.2f seconds' % (group_number,
                                                                  # toc))
            # nvars.steps_per_bead[this_step:, group_number] = 0
        # else:
            # log.info('minimization not setup or psf_flag=False')


def epilogue(other_self):
    mvars = other_self.mvars
    steps_per_bead = mvars.nvars.steps_per_bead
    accepted = other_self.mcvars.accepted
    steps_per_bead_out = os.path.join(other_self.runpath,
                                         mvars.dcdfile.replace('dcd', 'spb'))

    index = numpy.arange(1, accepted+1).reshape(accepted, 1)
    spb = numpy.concatenate((index, numpy.array(steps_per_bead[1:accepted+1])),
                            axis=1)
    fmt = '%09d'
    for i in xrange(mvars.number_of_flexible_regions):
        fmt += ' %0.6f'
    numpy.savetxt(steps_per_bead_out, spb, fmt=fmt,
                  header='dcd frame, steps per bead for each flexible region')


def step_counter_goback(other_self, goback_frame):
    """
    track the number of accepted steps (used to determine when to minimize)
    """
    steps_per_bead = other_self.mvars.nvars.steps_per_bead
    next_step = other_self.mcvars.accepted + 1
    steps_per_bead[next_step:] = steps_per_bead[goback_frame]


def get_boltzman_factor(energy_change, log):
    """
    given a change in energy, output the related Boltzmann probablity:
    probability = exp(-energy_change)

    note:
        The try/except and warnings/errors were implemented because some large
        changes in energy once created a problem, likely from overflow (though
        the issue may not have come up since, not certain)

    inputs:
        energy_change - the change in energy
        log - run log for error reporting

    returns:
        probability - probability of changing to the new energy state
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('error') # need this for np warnings
        try:
            probability = numpy.exp(-energy_change)
            log.debug('probability of new structure: %f' % probability)
        except Warning:
            if energy_change > 99:
                probability =  0
                log.debug('large energy increase, probability set to 0')
            elif energy_change < 0:
                probability =  1
                log.debug('energy decreased, probability set to 1')
        except:
                log.error("ERROR: Unexpected error calculating Boltzmann "
                "probability:" + sys.exc_info()[0])

    return probability


def recover_all_atom_dsNA(other_self, cg_flex_mol, flexible_mol, xyz_vectors,
                          all_beads, bead_masks):
    """
    Recover the all-atom representation of coarse-grained DNA
    """

    n_beads = numpy.copy(cg_flex_mol.natoms())
    bead_origins = numpy.copy(cg_flex_mol.coor()[0])

    # split vecXYZ into three matrices
    x_vec = numpy.copy(xyz_vectors[0])
    y_vec = numpy.copy(xyz_vectors[1])
    z_vec = numpy.copy(xyz_vectors[2])

    for i in xrange(n_beads):
        # align_to_origin transforms DNA beads to their original position
        # align_to_origin.transpose() transforms to the final position
        align_to_origin = align_to_xyz_vectors(other_self, x_vec[i, :],
                                               y_vec[i, :], z_vec[i, :])
        # reference frame of the bead
        align_to_final = align_to_origin.transpose()
        # translate to bp's final position
        align_to_final[3, :3] = bead_origins[i, :]

        rows, columns = all_beads[i].coor()[0].shape
        bead_coor = numpy.ones((rows, 4))     # initialize the bead_coor matrix
        orig_coor = numpy.ones((1, rows, 3))
        orig_coor = all_beads[i].coor()
        bead_coor[:, :3] = numpy.copy(orig_coor[0])
        bead_coor = numpy.dot(bead_coor, align_to_final)

        tmp_coor = numpy.zeros((1, rows, 3))          # initialize the tmp_coor
        tmp_coor[0] = bead_coor[:, :3]
        all_beads[i].setCoor(tmp_coor)

        #recombine all the beads back into one pdb
        error = flexible_mol.set_coor_using_mask(all_beads[i], 0,
                                                 bead_masks[i])
        handle_error(other_self, error)

        #reset the bead coordinates for next iteration
        all_beads[i].setCoor(orig_coor)


def align_to_xyz_vectors(other_self, x_vec, y_vec, z_vec):

    log = other_self.log
    zero = 1E-5

    tmp_coor = numpy.zeros((2, 4))
    tmp_coor[:, 3] = 1
    tmp_coor[1, 0:3] = z_vec
    align_z = align_to_z(other_self, tmp_coor)
    new_x_vec = numpy.dot(x_vec, align_z[0:3, 0:3])
    if new_x_vec[2] >= zero:
        log.error('ERROR: z-component of new_x_vec is %f and it should be 0'
                  % new_x_vec[2])

    theta_z = -numpy.arctan2(new_x_vec[1], new_x_vec[0]) # angle rotate about z

    rotate_z = rotate4x4(other_self, 'z', theta_z)

    align_to_origin = numpy.dot(align_z, rotate_z)

    # test that this align
    new_y_vec = numpy.dot(y_vec, align_to_origin[0:3, 0:3])

    if new_y_vec[0]+new_y_vec[2] > zero:
        log.error('ERROR: new_y_vec is not aligned to the y-axis as intended')

    return align_to_origin


def dna_bead_overlap(other_self, dna_coor, n_skip, cutoff):
    """
    this function uses FORTRAN executable to check for overlap between
    coordinates, skipping n_skip nearest-neighbors
    """

    log = other_self.log

    # make sure cg_DNA beads are further than the cutoff value
    # (exclude n_skip nearest neighbors)
    check, a, b = dna_overlap.overlap_skip_n(dna_coor, n_skip, cutoff)
    if check:
        a_b = dna_coor[a-1] - dna_coor[b-1]
        dist = numpy.sqrt(a_b.dot(a_b))
        log.debug('trial rejected: overlap between DNA beads')
        log.debug('coordinates %d and %d were %0.2f apart' %
                  (a-1, b-1, dist))

    return check


def two_body_overlap(other_self, mol_a_coor, mol_b_coor,
                     cutoff):
    """
    this function uses FORTRAN executable to check for overlap between
    two sets of coordinates
    """

    log = other_self.log
    na_coor = len(mol_a_coor)
    nb_coor = len(mol_b_coor)
    if na_coor > 0 and  nb_coor > 0:
        # make sure none of the distances between mol_a and mol_b coordinates
        # are less than the cutoff value
        check, a, b = dna_overlap.overlap2(mol_a_coor, mol_b_coor, cutoff)
        if check:
            coor_a = mol_a_coor[a-1]
            coor_b = mol_b_coor[b-1]
            ab = coor_b - coor_a
            dist = numpy.sqrt(ab.dot(ab))
            log.debug('trial rejected: collision between post atoms and fixed '
                      'DNA bead')
            log.debug('DNA bead %d and post coordinate # %d were %0.2f apart' %
                      (a-1, b-1, dist))
    else:
        check = 0
        log.debug('no coordinates to check for overlap: \n'
                  'natoms_mol_a = %d ; natoms_mol_b = %d' % (na_coor, nb_coor))

    return check


def get_move_to_origin_matrix(other_self, coor4):

    log = other_self.log

    translation_matrix = numpy.eye(4)
    reverse_translation_matrix = numpy.copy(translation_matrix)

    try:
        translation_matrix[3, 0:3] = -coor4[0, 0:3]
        reverse_translation_matrix[3, 0:3] = coor4[0, 0:3]
    except:
        log.debug('no coordinates to translate:\n', coor4)

    return (translation_matrix, reverse_translation_matrix)


def bead_rotate(other_self, coor3, xyz_vectors, theta_xyz, post_coor3):
    """
    this function is designed to generate a modified version of the input
    coordinates (coor3)
    1) translate all coordinates so the first one is at the origin
    2) orient the second coordinate along the z-axis
    3) performs rotations using the input thetas averaged over n_soft beads
    4) reverse the orienting to the z-axis
    5) reverse transaltion of all coordinates so the first is where it started
    """

    (n_atoms, col) = coor3.shape

    # initialize vector arrays for coordinates and orientation vectors
    # changing them from 3 component vectors into 4 component vectors
    # this incorporates transaltions into the matrix math
    flex_coor4 = numpy.ones((n_atoms, 4), numpy.float)
    post_coor4 = numpy.ones((len(post_coor3), 4), numpy.float)
    X = numpy.copy(flex_coor4)
    Y = numpy.copy(flex_coor4)
    Z = numpy.copy(flex_coor4)

    # populate the arrays with the input values
    flex_coor4[:, 0:3] = coor3
    post_coor4[:, 0:3] = post_coor3
    X[:, 0:3] = numpy.copy(xyz_vectors[0])
    Y[:, 0:3] = numpy.copy(xyz_vectors[1])
    Z[:, 0:3] = numpy.copy(xyz_vectors[2])


    # create the translation-rotation matrix
    # This is intended to be multiplied from the right (unlike standard matrix
    # multiplication) so as not to require transposing the coordinate vectors.
    theta_xyz_rad = theta_xyz * numpy.pi / 180.0 # radians

    [cx, cy, cz] = numpy.cos(theta_xyz_rad)
    [sx, sy, sz] = numpy.sin(theta_xyz_rad)

    # initialize the rotation
    # consolidated method of defining the rotation matrices
    rotate = numpy.eye(4, dtype=numpy.float)
    rotate[0][0:3] = [ cy*cz,          cy*sz,          -sy   ]
    rotate[1][0:3] = [ sx*sy*cz-cx*sz, sx*sy*sz+cx*cz, sx*cy ]
    rotate[2][0:3] = [ sx*sz+cx*sy*cz, cx*sy*sz-sx*cz, cx*cy ]

    (move_to_origin_0, return_from_origin_0) = get_move_to_origin_matrix(
        other_self, flex_coor4)
    # note: Ti0 != T0, negative off diag elements

    flex_coor4 = numpy.dot(flex_coor4, move_to_origin_0)
    z_align_matrix = align_to_z(other_self, flex_coor4)
    test = numpy.dot(flex_coor4, z_align_matrix)
    zero = 1E-10
    assert numpy.sum(test[1,:2]) < zero, '2nd coordinate not aligned to z-axis'
    align_and_rotate = numpy.dot(z_align_matrix, rotate)

    # first of n_soft rotations
    flex_coor4 = numpy.dot(flex_coor4, align_and_rotate)

    # the coarse grained beads local coordinates should not be translated,
    # only rotated
    X = numpy.dot(X, align_and_rotate)
    Y = numpy.dot(Y, align_and_rotate)
    Z = numpy.dot(Z, align_and_rotate)

    post_coor4 = numpy.dot(numpy.dot(post_coor4, move_to_origin_0),
                           align_and_rotate)

    flex_coor4 = numpy.dot(flex_coor4, z_align_matrix.transpose())
    flex_coor4 = numpy.dot(flex_coor4, return_from_origin_0)
    X = numpy.dot(X, z_align_matrix.transpose())
    Y = numpy.dot(Y, z_align_matrix.transpose())
    Z = numpy.dot(Z, z_align_matrix.transpose())

    post_coor4 = numpy.dot(post_coor4, z_align_matrix.transpose())
    post_coor4 = numpy.dot(post_coor4, return_from_origin_0)

    xyz_vectors[0, 1:] = X[1:, 0:3]
    xyz_vectors[1, 1:] = Y[1:, 0:3]
    xyz_vectors[2, 1:] = Z[1:, 0:3]

    # this returns the modified positions and orientations for all but the
    # first bead or reference bead
    return (flex_coor4[1:, 0:3], xyz_vectors[:, 1:], post_coor4[:, 0:3])


def rotate4x4(other_self, rotation_axis, theta):
    """
    rotate a 4x4 matrix an angle "theta" about the "rotation_axis"
    """

    log = other_self.log

    rotate = numpy.eye(4)
    ct = numpy.cos(theta)
    st = numpy.sin(theta)
    if rotation_axis.lower()=='x':
        (rotate[1, 1], rotate[1, 2]) = (ct, st)
        (rotate[2, 1], rotate[2, 2]) = (-st, ct)
    elif rotation_axis.lower()=='y':
        (rotate[0, 0], rotate[0, 2]) = (ct, -st)
        (rotate[2, 0], rotate[2, 2]) = (st, ct)
    elif rotation_axis.lower()=='z':
        (rotate[0, 0], rotate[0, 1]) = ( ct, st)
        (rotate[1, 0], rotate[1, 1]) = (-st, ct)
    else:
        log.error('ERROR: did not recognize rotation axis')

    return rotate


def align_to_z(other_self, coor4):
    """
    align the axis connecting the first 2 coordinates of a 1x4 array of
    coodinate vectors to the z-axis
    """

    log = other_self.log

    align = numpy.eye(4, dtype=numpy.float)

    assert all(coor4[0] == [0., 0., 0., 1., ]), (
        "coordinates passed to align2z were not translated to the origin")

    if coor4.shape > (1, 4):
        ## old method
        # (u, v, w) = coor4[1, 0:3]
        # zero = 1E-5 # small limit to make sure not dividing by zero

        # d1 = numpy.sqrt(u**2+v**2)
        # d2 = numpy.sqrt(u**2+v**2+w**2)
        # if d1 > zero:
            # align[0][0:3] = [ u/d1*w/d2, -v/d1, u/d2 ]
            # align[1][0:3] = [ v/d1*w/d2,  u/d1, v/d2 ]
            # align[2][0:3] = [ -d1/d2,        0, w/d2 ]

        # improved method
        v1 = numpy.array([0,0,1])
        v2 = coor4[1,0:3]
        align[:3,:3] = rotate_v2_to_v1(v1, v2).transpose()
    else:
        log.debug('no beads to align')

    return align


def setup_overlap(other_self, group_number, nvars):
    """
    method to setup overlap mask
    currently only 'all', 'backbone', and 'heavy' option are implemented
    all these options will have the same cutoff value of 0.8
    this overlap check maybe polished by further numerical studies
    """

    log = other_self.log
    mvars = other_self.mvars
    mol = other_self.group_molecules[group_number]

    n_flex_regions = mvars.number_of_flexible_regions
    post_mol = other_self.group_post_molecules[group_number]

    log.debug('in setup_overlap')

    overlap_basis = mvars.overlap_basis

    log.debug('creating overlap basis for group using %s' % overlap_basis)

    if(overlap_basis.lower() == 'all'):
        overlap_basis_filter = "not name None"
        mvars.cutoff = 0.8

    elif(overlap_basis.lower() == 'backbone'):
        overlap_basis_filter = ("name O2P or name O1P or name P or "
                                "name O5' or name C5' or name C4' or "
                                "name C3' or name O3'")
        mvars.cutoff = 0.8

    elif(overlap_basis.lower() == 'heavy'):
        overlap_basis_filter = "not name H"
        mvars.cutoff = 0.8

    overlap_basis_filter_python = basis_to_python.parse_basis(
        overlap_basis_filter)

    log.info('overlap_basis_filter_python = %s' % overlap_basis_filter_python)

    error, all_overlap_mask = mol.get_subset_mask(overlap_basis_filter_python)
    handle_error(other_self, error)
    if post_mol.natoms() > 0:
        error, post_overlap_mask = post_mol.get_subset_mask(
            overlap_basis_filter_python)
        handle_error(other_self, error)
    else:
        post_overlap_mask = ''

    try:
        nvars.all_overlap_mask
    except:
        nvars.all_overlap_mask = [None] * n_flex_regions
        nvars.post_overlap_mask = [None] * n_flex_regions

    nvars.all_overlap_mask[group_number] = all_overlap_mask
    nvars.post_overlap_mask[group_number] = post_overlap_mask


def handle_error(other_self, error):
    """
    code to handle error output from sassie modules
    e.g. copy_molecule_using_mask, get_subset_mask
    """

    if error:
        line_number = inspect.currentframe().f_back.f_lineno - 1
        info_string = ''
        for (i, error_string) in enumerate(error):
            info_string += error_string + ' '
        other_self.log.error('ERROR: on line %d: %s' % (line_number,
                                                        error_string))


def setup_double_stranded_nucleic_parameters(other_self, group_number, nvars):
    """
    this method uses the basis_string arrays to set up the coarse-grained DNA
    for later Monte Carlo trials
    """

    log = other_self.log
    log.debug('in setup_double_stranded_nucleic_parameters')

    mvars = other_self.mvars
    mol = other_self.group_molecules[group_number]
    flexible_mol = other_self.group_flexible_molecules[group_number]

    flexible_basis = mvars.basis_string_array_python[group_number]

    base_pairs_per_bead = 1

    n_flex_regions = mvars.number_of_flexible_regions

    dna_chain1_mol = system.Molecule(0)
    dna_chain2_mol = system.Molecule(0)
    coarse_grained_flexible_mol = system.Molecule(0)

    # coarse-grain the flexible_basis_mol
    tic = time.time()
    try:
        dna_chain_basis = nvars.dna_chain_basis
    except:
        # initialize data structures
        nvars.dna_chain_basis = dna_chain_basis = [None] * n_flex_regions
        nvars.base_pairs_in_each_bead = [None] * n_flex_regions
        nvars.xyz_vectors = [None] * n_flex_regions
        nvars.all_beads = [None] * n_flex_regions
        nvars.dna_segnames = [None] * n_flex_regions
        nvars.resids1 = [None] * n_flex_regions
        nvars.resids2 = [None] * n_flex_regions
        nvars.coarse_grained_flexible_mol = [None] * n_flex_regions
        nvars.bp_masks = [None] * n_flex_regions
        nvars.nucleic_acid_bead_masks = [None] * n_flex_regions
        nvars.number_beads = [0] * n_flex_regions
        nvars.steps_per_bead = numpy.zeros((mvars.trial_steps + 1,
                                            n_flex_regions))

    nvars.all_beads[group_number] = []
    bp_in_each_bead = nvars.base_pairs_in_each_bead

    # create a separate sasmol object for each DNA chain
    dna_chain_basis[group_number] = flexible_basis.split(' or ')
    error, dna_chain1_mask = flexible_mol.get_subset_mask(
        dna_chain_basis[group_number][0])
    handle_error(other_self, error)
    error, dna_chain2_mask = flexible_mol.get_subset_mask(
        dna_chain_basis[group_number][1])
    handle_error(other_self, error)
    error = flexible_mol.copy_molecule_using_mask(dna_chain1_mol,
                                                  dna_chain1_mask, 0)
    handle_error(other_self, error)
    error = flexible_mol.copy_molecule_using_mask(dna_chain2_mol,
                                                  dna_chain2_mask, 0)
    handle_error(other_self, error)

    if nvars.dna_segnames[group_number] is None:
        log.info("WARNING: dna_segnames not in input, may guess wrong")
        nvars.dna_segnames[group_number] = [dna_chain1_mol.segnames()[0],
                                            dna_chain2_mol.segnames()[0]]

    # check that there are the same number of bases in the paired chains
    nvars.resids1[group_number] = resids1 = numpy.array(
        dna_chain1_mol.resids())
    nvars.resids2[group_number] = resids2 = numpy.array(
        dna_chain2_mol.resids())
    if resids1.size != resids2.size:
        log.error("ERROR: number of paired nucleic acid bases for flexible"
                  " group number %d are not the same" % group_number)

    # make the DNA beads
    nvars.number_beads[group_number] = number_beads = int(numpy.round(
        resids1.size))
    # every bp has a C4' atom, using it to create dummy properties
    error, make_beads_mask = dna_chain1_mol.get_subset_mask(
        'name[i] == "C4\'" ')

    if number_beads != sum(make_beads_mask):
        log.error("ERROR: failed to create the correct number of "
                  "coarse-grained beads, expected to make %d but made %d"
                  % (number_beads, sum(make_beads_mask)))

    # create sasmol object for the cg beads and initialize bead coordinates
    error = dna_chain1_mol.copy_molecule_using_mask(
        coarse_grained_flexible_mol, make_beads_mask, 0)
    handle_error(other_self, error)
    nvars.coarse_grained_flexible_mol[group_number] = (
        coarse_grained_flexible_mol)

    bp_in_each_bead[group_number] = numpy.ones(number_beads, int)

    # create the bead masks if they are not already created
    get_bead_masks(other_self, group_number, nvars)

    # store the molecules for each bead
    for j in xrange(number_beads):
        bead = system.Molecule(0)
        error = flexible_mol.copy_molecule_using_mask(
            bead, nvars.nucleic_acid_bead_masks[group_number][j], 0)
        handle_error(other_self, error)
        nvars.all_beads[group_number].append(bead) # contains dummy coordinates

    # setup bead coordinate reference frame tracking system
    xyz_vectors = numpy.zeros((3, number_beads, 3))
    nvars.xyz_vectors[group_number] = xyz_vectors

    toc = time.time() - tic
    log.info('dsDNA coarse-graining took %0.3f seconds' % toc)


def get_bead_masks(other_self, group_number, nvars):
    """
    method for generating bead-masks
    """

    log = other_self.log

    # Get input:
    pgui            = other_self.run_utils.print_gui
    flexible_mol    = other_self.group_flexible_molecules[group_number]
    bp_in_each_bead = nvars.base_pairs_in_each_bead[group_number]
    number_beads    = nvars.number_beads[group_number]
    chain1_res_a    = nvars.resids1[group_number][0]
    chain2_res_b    = nvars.resids2[group_number][-1]

    number_base_pairs = nvars.resids1[group_number].size
    [nucleic_segname1, nucleic_segname2] = nvars.dna_segnames[group_number]

    sasmol_cpp_failed = False
    if use_sasmol_cpp:
        try:
            log.info('using c++ sasmol to CG dna')
            cpp_object = create_sasmol_as_cpp(flexible_mol)
            bead_masks = get_bead_masks_in_cpp(cpp_object, bp_in_each_bead,
                number_beads, chain1_res_a, chain2_res_b, number_base_pairs,
                nucleic_segname1, nucleic_segname2)
            pgui('successfully created DNA masks using sasmol cpp')
        except:
            pgui('sasmol cpp failed to create DNA masks, using python')
            sasmol_cpp_failed = True
    if sasmol_cpp_failed or not use_sasmol_cpp:
        log.info('using python sasmol to CG dna')
        bead_masks = []
        for j in xrange(number_beads):
            if bp_in_each_bead[j] > 1:
                # create the basis filters fore ach bead
                if j+1 == number_beads:
                    # accomodate for a non-divisible number of residues
                    bp_in_each_bead[j] = (number_base_pairs -
                                          bp_in_each_bead[j] * j)

                # Get the atoms from DNA strand 1
                chain1_res_b = chain1_res_a + bp_in_each_bead[j]
                chain1_bead_filter = ("((resid[i] >= " + str(chain1_res_a) +
                                      " and resid[i] < " + str(chain1_res_b) +
                                      ") and (segname[i] == '" +
                                      nucleic_segname1 + "')) or ")
                chain2_res_a = chain2_res_b - bp_in_each_bead[j]
                chain2_bead_filter = ("((resid[i] > " + str(chain2_res_a) +
                                      " and resid[i] <= " + str(chain2_res_b) +
                                      ") and (segname[i] == '" +
                                      nucleic_segname2 + "'))")

                # setup for next iteration
                chain1_res_a = chain1_res_b
                chain2_res_b = chain2_res_a
            else:
                chain1_bead_filter = ("((resid[i] == " + str(chain1_res_a) +
                                      ") and (segname[i] == '" +
                                      nucleic_segname1 + "')) or ")
                chain2_bead_filter = ("((resid[i] == " + str(chain2_res_b) +
                                      ") and (segname[i] == '" +
                                      nucleic_segname2 + "'))")

                # setup for next iteration
                chain1_res_a += 1
                chain2_res_b -= 1

            bead_filter = chain1_bead_filter + chain2_bead_filter

            # create the bead masks to select the atoms from the aa mol
            error, bead_mask = flexible_mol.get_subset_mask(bead_filter)
            handle_error(other_self, error)

            # store the mask for the reverse coarse-graining
            bead_masks.append(bead_mask)

    nvars.nucleic_acid_bead_masks[group_number] = bead_masks


def setup(other_self, group_number):
    """
    this method prepares a group for a flexible double stranded nucleic region
    """

    mvars = other_self.mvars
    log = other_self.log

    log.debug('in setup')

    # get the double stranded nucleic acid variables
    try:
        nvars = mvars.nvars
    except:
        nvars = double_stranded_nucleic_variables()

    flexible_basis = mvars.basis_string_array[group_number]
    log.info('initializing group = %d: %s' % (group_number, flexible_basis))

    setup_double_stranded_nucleic_parameters(other_self, group_number, nvars)

    setup_overlap(other_self, group_number, nvars)

    mvars.nvars = nvars


if __name__=="__main__":

    import doctest
    doctest.testmod()
