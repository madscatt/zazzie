# $Id: convergence_test.py 3188 2016-06-13 16:58:12Z schowell $
import glob
import logging
import numpy
import os
import pandas
import time
import sasmol.sasmol as sasmol

try:
    dummy = os.environ["DISPLAY"]
except:
    # allows for creating plots without an xserver
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

class AlignInputs(object):

    def __init__(self, goal_pdb, move, ref_pdb, out_fname, **kwargs):
        self.goal_pdb = goal_pdb
        self.ref_pdb = ref_pdb
        self.move = move
        self.out_fname = out_fname
        self.path = kwargs.get('path', './')
        self.basis_atoms = kwargs.get('basis_atoms', 'CA')
        self.seg_or_chain = kwargs.get('seg_or_chain', 'segname')
        self.seg_chain = kwargs.get('seg_chain', 'GAG')
        self.min_resid = kwargs.get('min_resid', 20)
        self.max_resid = kwargs.get('max_resid', 30)
        default_filter = ('(({}[i] == "{}") and (name[i] == "{}") and '
                          '(resid[i] >= {}) and (resid[i] <= {}))'.format(
                              self.seg_or_chain, self.seg_chain,
                              self.basis_atoms, self.min_resid, self.max_resid))
        self.goal_filter = kwargs.get('goal_filter', default_filter)
        self.move_filter = kwargs.get('move_filter', default_filter)
        logging.debug('goal_pdb: {}'.format(self.goal_pdb))
        logging.debug('ref_pdb: {}'.format(self.ref_pdb))
        logging.debug('move: {}'.format(self.move))
        logging.debug('out_fname: {}'.format(self.out_fname))
        logging.debug('path: {}'.format(self.path))
        logging.debug('goal_filter: {}'.format(self.goal_filter))
        logging.debug('move_filter: {}'.format(self.move_filter))

def align(inputs):
    '''
    input:
    ------
        inputs: object should contain the following attributes
            goal:    goal pdb
            ref:     reference pdb containing molecule info for moving pdb/dcd
            move:    pdb/dcd to align
            out:     output dcd file
            path:    output path
            goal_filter:     goal basis filter
            move_filter:     move basis filter

    note: inputs.ref and inputs.move can ofter be the same pdb
    '''

    aa_goal_pdb = inputs.goal_pdb
    aa_move_pdb = inputs.ref_pdb
    aa_move_fname = inputs.move
    save_fname = inputs.out_fname
    path = inputs.path

    if save_fname == aa_move_fname:
        in_place = True
        save_fname = 'temp' + save_fname[-4:]

    try:
        goal_filter = inputs.goal_filter
    except:
        basis_atoms = inputs.basis_atoms
        goal_seg_or_ch = inputs.goal_seg_or_chain
        goal_segname = inputs.goal_seg_chain
        goal_res_max = inputs.goal_max
        goal_res_min = inputs.goal_min

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

    # check input
    assert os.path.exists(aa_move_fname), ('ERROR: no such file - %s' %
                                          aa_move_fname)
    assert os.path.exists(aa_move_pdb), ('ERROR: no such file - %s' %
                                         aa_move_pdb)
    assert os.path.exists(aa_goal_pdb), ('ERROR: no such file - %s' %
                                         aa_goal_pdb)

    # create the SasMol objects
    sub_goal = sasmol.SasMol(0)
    sub_move = sasmol.SasMol(0)
    aa_goal = sasmol.SasMol(0)
    aa_move = sasmol.SasMol(0)

    aa_goal.read_pdb(aa_goal_pdb)
    aa_move.read_pdb(aa_move_pdb)

    if aa_move_fname[-3:] == 'pdb':
        aa_move.read_pdb(aa_move_fname)
        n_frames = aa_move.number_of_frames()
        in_type = 'pdb'
    elif aa_move_fname[-3:] == 'dcd':
        dcd_file = aa_move.open_dcd_read(aa_move_fname)
        n_frames = dcd_file[2]
        in_type = 'dcd'
    else:
        message = "\n~~~ ERROR, unknown input type ~~~\n"
        print_failure(message, txtOutput)
        return

    out_type = save_fname[-3:].lower()
    if 'dcd' == out_type:
        dcd_out_file = aa_move.open_dcd_write(path + save_fname)
    elif 'pdb' == out_type:
        dcd_out_file = None

    error, goal_seg_mask = aa_goal.get_subset_mask(goal_filter)
    assert not error, error
    error, move_seg_mask = aa_move.get_subset_mask(move_filter)
    assert not error, error

    error = aa_goal.copy_molecule_using_mask(sub_goal, goal_seg_mask, 0)
    assert not error, error
    error = aa_move.copy_molecule_using_mask(sub_move, move_seg_mask, 0)
    assert not error, error

    # calculate the center of mass of the subset of m1
    com_sub_goal = sub_goal.calccom(0)
    sub_goal.center(0)                         # center the m1 coordinates
    # get the m1 centered coordinates
    coor_sub_goal = sub_goal.coor()[0]

    for i in xrange(n_frames):
        if in_type == 'dcd':
            aa_move.read_dcd_step(dcd_file, i)
            # move m2 to be centered at the origin
            aa_move.center(0)
            error, sub_move.coor = aa_move.get_coor_using_mask(
                0, move_seg_mask)
            sub_move.setCoor(sub_move.coor)
            # calculate the center of mass of the subset of m2
            com_sub_move = sub_move.calccom(0)
            # move the subset of m2 to be centered at the origin
            sub_move.center(0)
            # get the new coordinates of the subset of m2
            coor_sub_move = sub_move.coor[0]
            # align m2 using the transformation from sub_m2 to sub_m1
            aa_move.align(
                0, coor_sub_move, com_sub_move, coor_sub_goal, com_sub_goal)
        elif in_type == 'pdb':
            # move m2 to be centered at the origin
            aa_move.center(i)
            error, sub_move.coor = aa_move.get_coor_using_mask(
                i, move_seg_mask)
            sub_move.setCoor(sub_move.coor)
            # calculate the center of mass of the subset of m2
            com_sub_move = sub_move.calccom(0)
            # move the subset of m2 to be centered at the origin
            sub_move.center(0)
            # get the new coordinates of the subset of m2
            coor_sub_move = sub_move.coor[0]
            # align m2 using the transformation from sub_m2 to sub_m1
            aa_move.align(
                i, coor_sub_move, com_sub_move, coor_sub_goal, com_sub_goal)

        aa_move.write_dcd_step(dcd_out_file, 0, i + 1)

    if in_type == 'dcd':
        aa_move.close_dcd_read(dcd_file[0])
    if out_type == 'dcd':
        aa_move.close_dcd_write(dcd_out_file)

    if in_place:
        os.remove(aa_move_fname)
        os.rename(save_fname, aa_move_fname)

    logging.info('Alingment of {} complete. \m/ >.< \m/'.format(aa_move_fname))

def calc_sas_convergence_all(sas_folders, output_prefix=None,
                             granularity=int(1e3), show=False, sas_ext='iq'):

    assert len(sas_folders) == 1, ("ERROR: mode for examining multiple "
                                   "folders requires additional debugging")

    if not output_prefix:
        output_prefix = 'sas_convergence'

    # initialize data sets
    iq_all = []
    list_new_grids = []
    list_occupied_grids = []

    n_q, n_spec = load_iq(sas_folders, sas_ext, iq_all)

    count_sas_grids(sas_folders, iq_all, n_q, n_spec,
                    list_new_grids, list_occupied_grids, granularity)

    total_spec = n_spec.sum()
    new_grids = numpy.zeros((total_spec, len(sas_folders)+1))
    new_grids[:, 0] = numpy.arange(total_spec)
    occupied_grids = numpy.copy(new_grids)
    for i in xrange(len(sas_folders)):
        rows = list_new_grids[i][:, 0] - 1
        new_grids[rows, 1] = list_new_grids[i][:, 1]
        occupied_grids[rows, 1] = list_occupied_grids[i][:, 1]

    # create output text files
    fname_occupied_grids = output_prefix + '_occupied_grids.npy'
    fname_new_grids = output_prefix +'_new_grids.npy'
    numpy.savetxt(fname_occupied_grids, occupied_grids)
    numpy.savetxt(fname_new_grids, new_grids)
    print 'output text files: \n%s \n%s' % (fname_occupied_grids,
                                            fname_new_grids)

    plot_convergence(new_grids, sas_folders, occupied_grids,
                     output_prefix, show, spatial=False)


def calc_sas_convergence_by_run(sas_folders, output_prefix=None,
                                granularity=int(1e3), show=False, sas_ext='iq'):

    assert len(sas_folders) == 1, ("ERROR: mode for examining multiple "
                                   "folders requires additional debugging")

    if not output_prefix:
        output_prefix = 'sas_convergence'

    # initialize data sets
    iq_all = []
    list_new_grids = []
    list_occupied_grids = []

    n_q, n_spec = load_iq(sas_folders, sas_ext, iq_all)

    count_sas_grids(sas_folders, iq_all, n_q, n_spec,
                    list_new_grids, list_occupied_grids, granularity)

    total_spec = n_spec.sum()
    new_grids = numpy.zeros((total_spec, len(sas_folders)+1), dtype=int)
    new_grids[:, 0] = numpy.arange(total_spec)
    occupied_grids = numpy.copy(new_grids)
    for i in xrange(len(sas_folders)):
        rows = list_new_grids[i][:, 0] -1
        new_grids[rows, i+1] = list_new_grids[i][:, 1]
        occupied_grids[rows, i+1] = list_occupied_grids[i][:, 1]

    # create output text files
    fname_occupied_grids = output_prefix + '_occupied_grids_by_run.npy'
    fname_new_grids = output_prefix +'_new_grids_by_run.npy'
    numpy.savetxt(fname_occupied_grids, occupied_grids)
    numpy.savetxt(fname_new_grids, new_grids)
    print 'output text files: \n%s \n%s' % (fname_occupied_grids,
                                            fname_new_grids)

    plot_convergence(new_grids, sas_folders, occupied_grids,
                     output_prefix, show, spatial=False)


def calc_spatial_convergence_all(pdb_fname, dcd_fnames, output_prefix=None,
                                 show=False, **kwargs):

    assert len(sas_folders) == 1, ("ERROR: mode for examining multiple "
                                   "folders requires additional debugging")

    if not output_prefix:
        output_prefix = pdb_fname[:-4]

    # initialize data sets
    list_new_voxels = []
    list_occupied_voxels = []

    count_spatial_voxels(pdb_fname, dcd_fnames, list_new_voxels,
                         list_occupied_voxels, **kwargs)

    n_structures = sum([len(new_voxels) for new_voxels in list_new_voxels])
    new_voxels = numpy.zeros((n_structures, 2))
    occupied_voxels = numpy.zeros((n_structures, 2))
    new_voxels[:, 0] = numpy.arange(n_structures)
    occupied_voxels[:, 0] = numpy.arange(n_structures)
    for i in xrange(len(dcd_fnames)):
        rows = list_new_voxels[i][:, 0] - 1
        new_voxels[rows, 1] = list_new_voxels[i][:, 1]
        occupied_voxels[rows, 1] = list_occupied_voxels[i][:, 1]

    # create output text files
    fname_occupied_voxels = output_prefix + '_occupied_voxels.npy'
    fname_new_voxels = output_prefix +'_new_voxels.npy'
    numpy.savetxt(fname_occupied_voxels, occupied_voxels)
    numpy.savetxt(fname_new_voxels, new_voxels)
    print 'output text files: \n%s \n%s' % (fname_occupied_voxels,
                                            fname_new_voxels)


    plot_convergence(new_voxels, dcd_fnames, occupied_voxels,
                     output_prefix, show)


def calc_spatial_convergence_by_run(pdb_fname, dcd_fnames, output_prefix=None,
                                    show=False, **kwargs):

    assert len(sas_folders) == 1, ("ERROR: mode for examining multiple "
                                   "folders requires additional debugging")

    if not output_prefix:
        output_prefix = pdb_fname[:4]

    # initialize data sets
    list_new_voxels = []
    list_occupied_voxels = []

    count_spatial_voxels(pdb_fname, dcd_fnames, list_new_voxels,
                         list_occupied_voxels, **kwargs)

    n_structures = sum([len(new_voxels) for new_voxels in list_new_voxels])
    new_voxels = numpy.zeros((n_structures, len(dcd_fnames)+1))
    occupied_voxels = numpy.zeros((n_structures, len(dcd_fnames)+1))
    new_voxels[:, 0] = numpy.arange(n_structures)
    occupied_voxels[:, 0] = numpy.arange(n_structures)
    for i in xrange(len(dcd_fnames)):
        rows = list_new_voxels[i][:, 0] -1
        new_voxels[rows, i+1] = list_new_voxels[i][:, 1]
        occupied_voxels[rows, i+1] = list_occupied_voxels[i][:, 1]

    # create output text files
    fname_occupied_voxels = output_prefix + '_occupied_voxels_by_run.npy'
    fname_new_voxels = output_prefix +'_new_voxels_by_run.npy'
    numpy.savetxt(fname_occupied_voxels, occupied_voxels)
    numpy.savetxt(fname_new_voxels, new_voxels)
    print 'output text files: \n%s \n%s' % (fname_occupied_voxels,
                                            fname_new_voxels)

    plot_convergence(new_voxels, dcd_fnames, occupied_voxels,
                     output_prefix, show)


def count_new_spatial_voxels(coors, voxel_set, delta):
    number_new_voxels = 0
    for coor in coors:
        voxel_number = get_spatial_voxel_number(coor, delta)
        if voxel_number not in voxel_set:
            number_new_voxels += 1
            voxel_set.add(voxel_number)
    return number_new_voxels


def count_sas_grids(sas_folders, iq_all, n_q, n_spec, list_new_grids,
                    list_occupied_grids, granularity=int(1e3), iq_low=0, iq_high=2):

    den = float(iq_high - iq_low)
    # grid = numpy.zeros((n_q, granularity+1))
    delta_i = 1.0/granularity # using I(0) = 1 as the default
    number_of_occupied_grids = 0
    cwd = os.getcwd()

    tic = time.time()
    for (i_folder, this_folder) in enumerate(sas_folders):
        logging.info('processing spec files from: {}\n'.format(this_folder))

        output_prefix = os.path.join(cwd, this_folder, '{}_of_{}'.format(
            i_folder+1, len(sas_folders)) )
        output_new_grids = output_prefix + '_new_grids.npy'
        output_occupied_grids = output_prefix + '_occupied_grids.npy'

        try:
            # try loading output from previous run
            this_folder_new_grids = numpy.load(output_new_grids)
            this_folder_occupied_grids = numpy.load(output_occupied_grids)
            logging.info('Successfully loaded new voxels and occupied '
                         'voxels for {} from:\n{} \n{}'.format(
                             this_folder, output_new_grids,
                             output_occupied_grids))
        except:
            # calculate and create output
            logging.info('Calculating convergence. Did not find output '
                         'files from previous calculation. Storing the output '
                         'to:\n%s \n%s' % (output_new_grids,
                                           output_occupied_grids))

            this_folder_new_grids = numpy.zeros((n_spec[i_folder], 2), dtype=int)
            this_folder_new_grids[:, 0] = numpy.arange(n_spec[i_folder]) + 1
            this_folder_occupied_grids = numpy.copy(this_folder_new_grids)
            occupied_grids = {}

            # convert I(Q) to bin number
            binned_iqs = numpy.array((iq_all[i_folder] - 1.0)/delta_i, dtype=int)

            for i_spec in xrange(n_spec[i_folder]):
                number_of_new_grids = 0
                for q in xrange(n_q):
                    grids_this_q = occupied_grids.get(q, {})
                    if not grids_this_q.get(binned_iqs[i_spec, q], 0):
                        grids_this_q[binned_iqs[i_spec, q]] = 1
                        number_of_new_grids += 1
                        occupied_grids[q] = grids_this_q
                number_of_occupied_grids += number_of_new_grids
                this_folder_occupied_grids[i_spec, 1] = number_of_occupied_grids
                this_folder_new_grids[i_spec, 1] = number_of_new_grids

            # print "temporarily not saving output"
            numpy.save(output_new_grids, this_folder_new_grids)
            numpy.save(output_occupied_grids, this_folder_occupied_grids)

        list_new_grids.append(this_folder_new_grids)
        list_occupied_grids.append(this_folder_occupied_grids)

    toc = time.time() - tic
    logging.info("time used: {}".format(toc))


def old_count_sas_grids(sas_folders, iq_low, iq_high, iq_all, n_q, n_spec,
                        list_new_grids, list_occupied_grids, n_grids):

    iq_low = numpy.array(iq_low).min(axis=0)
    iq_high = numpy.array(iq_high).max(axis=0)
    grid = numpy.zeros((n_q, n_grids+1))
    number_of_occupied_grids = 0
    i_spec = 0

    cwd = os.getcwd()

    tic = time.time()
    for (i_folder, this_folder) in enumerate(sas_folders):
        print 'processing spec files from: %s\n' % this_folder
        output_prefix = os.path.join(cwd, this_folder, '%d_of_%d' %
                                     (i_folder+1, len(sas_folders)) )
        output_new_grids = output_prefix + '_new_grids.npy'
        output_occupied_grids = output_prefix + '_occupied_grids.npy'

        try:
            # try loading output from previous run
            this_folder_new_grids = numpy.load(output_new_grids)
            this_folder_occupied_grids = numpy.load(output_occupied_grids)
            print ('Successfully loaded new voxels and occupied voxels '
                   'for %s from:\n%s \n%s' % (this_folder,
                                              output_new_grids,
                                              output_occupied_grids))
        except:
            # calculate and create output
            print ('Calculating convergence. Did not find output files from '
                   'previous calculation. Storing the output to:\n%s \n%s' % (
                       output_new_grids,
                       output_occupied_grids))

            this_folder_new_grids = numpy.zeros((n_spec[i_folder], 2),
                                                dtype=int)
            this_folder_occupied_grids = numpy.zeros((n_spec[i_folder], 2),
                                                     dtype=int)

            for i_spec_folder in xrange(n_spec[i_folder]):
                number_of_new_grids = 0
                for q in xrange(n_q):
                    num = iq_all[i_folder][q, i_spec_folder] - iq_low[q]
                    den = iq_high[q] - iq_low[q]
                    try:
                        n = int(n_grids * (num / den))
                    except ValueError:
                        n = int(numpy.nan_to_num(n_grids * (num / den)))
                    if not grid[q, n]:
                        grid[q, n] = 1
                        number_of_new_grids += 1
                number_of_occupied_grids += number_of_new_grids
                this_folder_new_grids[i_spec_folder, :] = [
                    i_spec, number_of_new_grids]
                this_folder_occupied_grids[i_spec_folder, :] = [
                    i_spec, number_of_occupied_grids]
                i_spec += 1

            numpy.save(output_new_grids, this_folder_new_grids)
            numpy.save(output_occupied_grids, this_folder_occupied_grids)

        list_new_grids.append(this_folder_new_grids)
        list_occupied_grids.append(this_folder_occupied_grids)

    toc = time.time() - tic
    print "time used: ", toc


def count_spatial_voxels(pdb_fname, dcd_fnames, list_new_voxels,
                         list_occupied_voxels, voxel_size=5.0,
                         basis_filter=None, filter_label='', align_dcd=False,
                         **kwargs):

    # initialize molecule and mask
    mol=sasmol.SasMol(0)
    mol.read_pdb(pdb_fname)
    n_dcds = len(dcd_fnames)
    cap_filter = '(name[i]=="CA" or name[i]=="P")'
    if basis_filter:
        error, mask = mol.get_subset_mask('%s and %s' % (
            basis_filter, cap_filter))
    else:
        error, mask = mol.get_subset_mask(cap_filter)
    assert not error, error
    voxel_set = set([])
    number_occupied_voxels = 0

    tic = time.time()
    # i_frame = 0
    for (i_dcd, dcd_fname) in enumerate(dcd_fnames):
        print 'processing dcd: %s\n' % dcd_fname
        dcd_output_prefix = '%s_%d_of_%d' % (dcd_fname[:-4], i_dcd+1,
                                             n_dcds)
        output_new_voxels = '%s%s_new_voxels.npy' % (
            dcd_output_prefix, filter_label)
        output_occupied_voxels = '%s%s_occupied_voxels.npy' % (
            dcd_output_prefix, filter_label)

        try:
            # try loading output from previous run
            this_dcd_new_voxels = numpy.load(output_new_voxels)
            this_dcd_occupied_voxels = numpy.load(output_occupied_voxels)
            print ('Successfully loaded new voxels and occupied voxels '
                   'for %s from:\n%s \n%s' % (dcd_fname,
                                              output_new_voxels,
                                              output_occupied_voxels))
        except:
            # calculate and create output
            print ('Calculating convergence. Did not find output files from '
                   'previous calculation. Storing the output to:\n%s \n%s' % (
                       output_new_voxels,
                       output_occupied_voxels))

            if align_dcd:
                inputs = AlignInputs(pdb_fname, dcd_fname, pdb_fname,
                                     dcd_fname, **kwargs)
                align(inputs)

            dcd_file = mol.open_dcd_read(dcd_fname)
            number_of_frames = dcd_file[2]
            this_dcd_new_voxels = numpy.zeros((number_of_frames, 2),dtype=int)
            this_dcd_new_voxels[:, 0] = numpy.arange(n_spec[i_folder]) + 1
            this_dcd_occupied_voxels = numpy.copy(this_dcd_new_voxels)

            for nf in xrange(number_of_frames):
                mol.read_dcd_step(dcd_file, nf)
                error, coor = mol.get_coor_using_mask(0, mask)
                assert not error, error
                number_new_voxels = count_new_spatial_voxels(coor[0], voxel_set,
                                                             voxel_size)
                number_occupied_voxels += number_new_voxels
                this_dcd_occupied_voxels[nf, 1] = number_occupied_voxels
                this_dcd_new_voxels[nf, 1] = number_new_voxels
                # i_frame += 1

            numpy.save(output_new_voxels, this_dcd_new_voxels)
            numpy.save(output_occupied_voxels, this_dcd_occupied_voxels)

        list_new_voxels.append(this_dcd_new_voxels)
        list_occupied_voxels.append(this_dcd_occupied_voxels)

    toc = time.time() - tic
    logging.info("time used: {}".format(toc))


def get_spatial_voxel_number(coor, delta):
    idx = int(coor[0]/delta)
    idy = int(coor[1]/delta)
    idz = int(coor[2]/delta)
    return (idx, idy, idz)


def load_iq(sas_folders, sas_ext, iq_all):
    n_folders = len(sas_folders)

    n_q = numpy.zeros(n_folders, dtype=int)
    n_spec = numpy.zeros(n_folders, dtype=int)

    cwd = os.getcwd()

    for (i_folder, this_folder) in enumerate(sas_folders):
        logging.info('loading spec files from: {}'.format(this_folder))
        output_prefix = os.path.join(cwd, this_folder, '{}_of_{}'.format(
            i_folder+1, n_folders))
        output_iq = output_prefix + '_iq.h5'

        sas_search_path = os.path.join(cwd, this_folder, '*.' + sas_ext)
        file_list = glob.glob(sas_search_path)

        n_spec[i_folder] = len(file_list)

        if n_spec[i_folder] < 1:
            logging.info('No I(Q) files found in: {}'.format(sas_search_path))
        else:
            try:
                # try loading iq_array from previous run
                store = pandas.HDFStore(output_iq)
                these_iqs_df = store['iq']
                q_vals = store['q']
                n_q[i_folder] = len(q_vals)
                these_iqs = numpy.array(these_iqs_df)
                logging.info('Successfully loaded iq_array for {} from:\n{}'.format(
                    this_folder, output_iq))
            except:
                logging.info('Loading in iq data from {}. Output stored to:\n{}'.format(
                    this_folder, output_iq))

                file_list.sort()
                ref_iq = numpy.loadtxt(file_list[0])
                q_vals = pandas.Series(ref_iq[:, 0])
                n_q[i_folder] = len(q_vals)
                these_iqs = numpy.zeros((n_spec[i_folder], n_q[i_folder]))

                for (j, this_file) in enumerate(file_list):
                    this_iq = numpy.loadtxt(this_file)
                    if not numpy.all( 0.0 == (this_iq[:,0] - q_vals) ):
                        logging.error('Q values do not match for iq file: {0}'.format(
                            iq_file))

                    these_iqs[j] = this_iq[:, 1] / this_iq[0,1] # I(0) = 1

                these_iqs_df = pandas.DataFrame(these_iqs, columns=q_vals)
                store['iq'] = these_iqs_df
                store['q'] = q_vals

            store.close()

            iq_all.append(these_iqs)

            assert n_q[i_folder] == n_q[0], (
                'ERROR: inconsistent number of Q-grid points between spec '
                'files in %s and %s' % (sas_folders[0], this_folder)
            )
    n_q = n_q[0]

    return n_q, n_spec


def plot_convergence(new_voxels, dcd_fnames, occupied_voxels,
                     output_prefix, show=False, spatial=True):
    fig = plt.figure(figsize=(6, 10))
    gs = gridspec.GridSpec(2, 1, left=0.1, right=0.9, wspace=0, hspace=0)
    ax=[]
    ax.append(plt.subplot(gs[0]))
    ax.append(plt.subplot(gs[1]))
    n_plots = new_voxels.shape[1] - 1
    for i in xrange(n_plots):
        if 1 < n_plots < 100:
            label = dcd_fnames[i]
        else:
            label = ''
        if i > 0:
            # rows = (new_voxels[:, i+1] > 0)
            ax[0].plot(new_voxels[1:, 0], new_voxels[1:, i+1],
                       label=label)
        else:
            # rows = (new_voxels[:, i+1] > 0)[1:] # skip the initial frame
            ax[0].plot(new_voxels[1:, 0], new_voxels[1:, i+1],
                       label=label)
    ax[0].xaxis.set_ticklabels([])
    if n_plots > 1 :
        lg = ax[0].legend(bbox_to_anchor=(1, 1), loc=2)
        # lg.draw_frame(False)

    for i in xrange(n_plots):
        if i > 0:
            rows = (occupied_voxels[:, i+1] > 0) # only plot non-zero values
            ax[1].plot(occupied_voxels[rows, 0], occupied_voxels[rows, i+1])
        else:
            rows = (occupied_voxels[:, i+1] > 0)[1:] # skip the initial frame
            ax[1].plot(occupied_voxels[rows, 0], occupied_voxels[rows, i+1])
    ax[1].set_xlabel('Structures')
    ylim = ax[1].get_ylim()
    ax[1].set_ylim((ylim[0], ylim[1] * 1.1))

    if spatial:
        ax[1].set_ylabel('Number of Occupied Voxels')
        ax[0].set_ylabel('Number of New Voxels')
    else:
        ax[1].set_ylabel('Number of Occupied Grids')
        ax[0].set_ylabel('Number of New Grids')

    plot_name = output_prefix + '_convergence'
    plot_name = os.path.join(os.getcwd(), plot_name)
    plt.savefig(plot_name + '.eps', dpi=400, bbox_inches='tight')
    plt.savefig(plot_name + '.png', dpi=400, bbox_inches='tight')
    print 'Saving figure to: \nevince %s.eps &\neog %s.png &' % (plot_name,
                                                              plot_name)
    if show:
        plt.show()
    else:
        plt.close('all')

if __name__=='__main__':
    import sys

    mol=sasmol.SasMol(0)
    if len(sys.argv)<3:
        mol.read_pdb('min_dsDNA60.pdb')
        # mol.read_dcd('run3_100k_ngb/monte_carlo/min_dsDNA60.dcd')
        dcd_full_name = 'run3_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd'
    else:
        mol.read_pdb(sys.argv[1])
        dcd_full_name = sys.argv[2]

    voxel_set = set([])
    delta = 5.0
    list_number_new_voxels = []
    list_number_occupied_voxels = []
    number_occupied_voxels = 0
    error, mask = mol.get_subset_mask('name[i]=="CA" or name[i]=="P"')

    dcd_file = mol.open_dcd_read(dcd_full_name)
    number_of_frames = dcd_file[2]

    tic = time.time()
    output_file = "number_of_occupied_voxels.txt"
    fout = open(output_file, 'w')
    fout.write("#frame_number, number_of_occupied_voxels\n")

    for nf in xrange(number_of_frames):
        mol.read_dcd_step(dcd_file, nf+1)
        error, coors = mol.get_coor_using_mask(0, mask)
        assert not error, error
        number_new_voxels = count_new_spatial_voxels(coors[0], voxel_set, delta)
        number_occupied_voxels += number_new_voxels
        list_number_new_voxels.append(number_new_voxels)
        list_number_occupied_voxels.append(number_occupied_voxels)
        fout.write("%d %d\n"%(nf, number_occupied_voxels))
    fout.close()
    toc = time.time() - tic
    print "\ntime used: ", toc

    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(2, 1, left=0.2, right=0.95, wspace=0, hspace=0)
    ax=[]
    ax.append(plt.subplot(gs[0]))
    ax.append(plt.subplot(gs[1]))
    ax[0].plot(range(len(list_number_new_voxels)), list_number_new_voxels)
    ax[0].set_xlabel('Structure')
    ax[0].set_ylabel('number of new voxels')
    ax[0].set_yscale('log') #lim([0, max(list_number_new_voxels)*1.05])
    ax[0].xaxis.set_ticklabels([])

    ax[1].plot(range(len(list_number_occupied_voxels)), list_number_occupied_voxels)
    ax[1].set_xlabel('Structure')
    ax[1].set_ylabel('number of occupied voxels')
    ylim = ax[1].get_ylim()
    ax[1].set_ylim((ylim[0], ylim[1] * 1.1))
    plt.savefig('metric_convergence.eps', dpi=400, bbox_inches='tight')
    plt.savefig('metric_convergence.png', dpi=400, bbox_inches='tight')
    plt.show()

    print '\m/ >.< \m/'
