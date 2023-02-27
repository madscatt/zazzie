# $Id: cluster_ensemble.py 3078 2016-04-06 19:46:43Z schowell $
import numpy
import time
import os
import sasmol.sasmol as sasmol
import sassie.calculate.convergence_test as convergence_test
try:
    dummy = os.environ["DISPLAY"]
except:
    # allows for creating plots without an xserver
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

def com_clustering(pdb_file_name, dcd_file_names, output_prefix=None,
                   voxel_size=5.0, show=False, create=False):
    '''
    Clustering routine for identifying grouping structures within
    an ensemble. Uses the center of mass of the entire structure to
    differentiate one from another.

    This may be improved by using the center of mass of a certain region.
    '''
    if output_prefix:
        output_prefix += '_%dA' % int(voxel_size)
    else:
        output_prefix = '%s_cluster_%dA' % (pdb_file_name[:-4],
                                            int(voxel_size))

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file_name)

    com_coors = []
    new_voxels = []
    occupied_voxels = []
    number_occupied_voxels = 0
    voxel_set = set([])
    tic = time.time()
    if create:
        subset_mol = sasmol.SasMol(0)
        subset_mol.read_pdb(pdb_file_name)
        output_dcd = subset_mol.open_dcd_write(output_prefix + '.dcd')
        i = 0

    for (i_dcd, dcd_file_name) in enumerate(dcd_file_names):
        print 'processing dcd: %s\n' % dcd_file_name
        input_dcd = mol.open_dcd_read(dcd_file_name)
        number_of_frames = input_dcd[2]

        for frame in xrange(number_of_frames):
            mol.read_dcd_step(input_dcd, frame)
            com_coors.append(mol.calccom(0))
        mol.close_dcd_read(input_dcd[0])

        this_dcd_new_voxels = convergence_test.count_new_spatial_voxels(
            com_coors, voxel_set, voxel_size)
        number_occupied_voxels += this_dcd_new_voxels
        new_voxels.append(this_dcd_new_voxels)
        occupied_voxels.append(number_occupied_voxels)

        # select out the dcd_frames and save them to a new dcd file
        if create:
            voxel_list = list(voxel_set)
            voxel_number = numpy.array(com_coors)/voxel_size
            dcd_frames = []
            for voxel in voxel_list:
                v_diff = voxel_number-voxel
                i_min = numpy.argmin(numpy.sqrt(v_diff[:,0]**2 + v_diff[:,1]**2
                                                + v_diff[:,2]**2))
                dcd_frames.append(i_min)

            dcd_frames_file = dcd_file_name.replace('.dcd', '_cluster_%dA.dat'
                                                    % voxel_size)
            numpy.savetxt(dcd_frames_file, dcd_frames, fmt='%d')

            input_dcd = subset_mol.open_dcd_read(dcd_file_name)
            number_of_frames = input_dcd[2]
            for frame in xrange(number_of_frames):
                subset_mol.read_dcd_step(input_dcd, frame)
                if frame in dcd_frames:
                    i += 1
                    subset_mol.write_dcd_step(output_dcd, 0, i)

            subset_mol.close_dcd_read(input_dcd[0])
    if create:
        subset_mol.close_dcd_write(output_dcd)
    toc = time.time() - tic
    print "\ntime used: ", toc

    # convergence_test.plot_convergence(new_voxels, dcd_file_names,
                                      # occupied_voxels, output_prefix)

    return number_occupied_voxels

def voxel_scan(pdb_file_name, dcd_file_names, voxel_range,
                 output_prefix=None, show=False, create=False):

    if not output_prefix:
        output_prefix = '%s_cluster' % pdb_file_name[:-4]

    occupied_voxels = numpy.zeros((len(voxel_range), 2), dtype='int64')

    for (i, voxel_size) in enumerate(voxel_range):
        occupied_voxels[i, 0] = voxel_size
        occupied_voxels[i, 1] = com_clustering(pdb_file_name, dcd_file_names,
                                               output_prefix=output_prefix,
                                               voxel_size=voxel_size,
                                               show=show, create=create)

    output_prefix = os.path.join(os.getcwd(), output_prefix)
    out_name = output_prefix + '.dat'
    numpy.savetxt(out_name, occupied_voxels, fmt='%d')

    ax = plt.subplot(111)
    plt.plot(occupied_voxels[:,0], occupied_voxels[:,1])
    ax.set_yscale('log')
    plt.xlabel(r'voxel size ($\AA$)')
    plt.ylabel('occupied voxels')
    plot_name = output_prefix
    plt.savefig(plot_name + '.eps', dpi=400, bbox_inches='tight')
    plt.savefig(plot_name + '.png', dpi=400, bbox_inches='tight')
    if show:
        plt.show()

