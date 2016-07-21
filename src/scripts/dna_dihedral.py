#!/usr/bin/python
#
# Author:  Steven C. Howell
# Purpose: methods for calculating DNA dihedral angles and generating plots
# Created: 22 September 2014
#
# $Id: dna_dihedral.py 2907 2015-11-06 19:20:58Z curtisj $
#
# 000000011111111112222222222333333333344444444445555555555666666666677777777778
# 345678901234567890123456789012345678901234567890123456789012345678901234567890

import sys
import logging
import numpy
import pandas as pd
import os.path as op
import matplotlib
import matplotlib.pyplot as plt
import sasmol.sasmol as sasmol
# import sassie_2_na.dihedral_angles as dihedral_angles
import sassie.calculate.dihedral_angles as dihedral_angles

def plot_me(x_angle, y_angle, x_angle_label, y_angle_label):

    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    ax1 = plt.subplot(111)
    ax1.scatter(x_angle, y_angle)

    ax1.set_xlim([-180, 180])
    ax1.set_ylim([-180, 180])

    ax1.set_xlabel(x_angle_label)
    ax1.set_ylabel(y_angle_label)
    plt.title(y_angle_label + ' vs ' + x_angle_label)

    plt.show()

def positive(angle):
    if angle < 0:
        return angle + 360


def plot_dna_dihedral_grid(all_angles):

    plt.ion()
    (n_frames, n_bps, n_angles) = all_angles.shape

    angles = {r'$\alpha$': 1, r'$\beta$': 2, r'$\gamma$': 3, r'$\delta$': 4,
              r'$\epsilon$': 5, r'$\zeta$': 6, r'$\chi$': 7}

    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    for frame in xrange(n_frames):
        plt.figure()

        plt.subplot(331)
        x = r'$\zeta$'
        y = r'$\alpha$'

        zero_to_threesixty(all_angles[frame])
        plt.plot(all_angles[frame][0:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')

        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelbottom='off', ) #, labeltop='on')

        plt.subplot(332)
        x = r'$\zeta$'
        y = r'$\beta$'
        plt.plot(all_angles[frame][0:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off', ) #, labeltop='on')

        plt.subplot(333)
        x = r'$\zeta$'
        y = r'$\epsilon$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(
            labelleft='off', labelbottom='off', ) #, labeltop='on', labelright='on')

        plt.subplot(334)
        x = r'$\alpha$'
        y = r'$\gamma$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelbottom='off')

        plt.subplot(335)
        x = r'$\zeta$'
        y = r'$\chi$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off')

        plt.subplot(336)
        x = r'$\delta$'
        y = r'$\chi$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off') #, labelright='on')

        plt.subplot(337)
        x = r'$\zeta$'
        y = r'$\zeta$'
        plt.plot(all_angles[frame][:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()

        plt.subplot(338)
        x = r'$\epsilon$'
        y = r'$\epsilon$'
        lines = []
        labels = []
        plt.plot(all_angles[frame][:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off')

        plt.suptitle('Starting Structure')
        plt.savefig('DNA-dihedrals_0steps.pdf', bbox_ingches='tight')
        plt.savefig('DNA-dihedrals_0steps.png', bbox_ingches='tight')

        plt.draw()

    plt.show()


def plot_dna_compare_dihedral(n0_all_angles, nm_all_angles, m__all_angles):
    angles = {'alpha': 0, 'beta': 1, 'gamma': 2, 'delta': 3,
              'epsilon': 4, 'zeta': 5, 'chi': 6}
    zero_to_threesixty(n0_all_angles)
    zero_to_threesixty(nm_all_angles)
    zero_to_threesixty(m_all_angles)

    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    plt.title('Scatter plots of selected torsional angles.')

    ax1 = plt.subplot(331)
    x = 'zeta'
    y = 'alpha'
    ax1.scatter(n0_all_angles[angles[x]][
                0:-1], n0_all_angles[angles[y]][1:], c='b', marker='x', s=80)
    ax1.scatter(nm_all_angles[angles[x]][
                0:-1], nm_all_angles[angles[y]][1:], c='r', marker='o', s=25)
    ax1.scatter(m__all_angles[angles[x]][
                0:-1], m__all_angles[angles[y]][1:], c='g', marker='s', s=25)
    ax1.set_xlabel(x)
    ax1.set_ylabel(y + ' + 1')
    grid_format(ax1)

    ax2 = plt.subplot(332)
    x = 'zeta'
    y = 'beta'
    ax2.scatter(n0_all_angles[angles[x]][
                0:-1], n0_all_angles[angles[y]][1:], c='b', marker='x', s=80)
    ax2.scatter(nm_all_angles[angles[x]][
                0:-1], nm_all_angles[angles[y]][1:], c='r', marker='o', s=25)
    ax2.scatter(m__all_angles[angles[x]][
                0:-1], m__all_angles[angles[y]][1:], c='g', marker='s', s=25)
    ax2.set_xlabel(x)
    ax2.set_ylabel(y + ' + 1')
    grid_format(ax2)

    all_angles = n0_all_angles

    ax3 = plt.subplot(333)
    x = 'zeta'
    y = 'epsilon'
    ax3.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax3.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax3.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax3.set_xlabel(x)
    ax3.set_ylabel(y)
    grid_format(ax3)

    ax4 = plt.subplot(334)
    x = 'alpha'
    y = 'gamma'
    ax4.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax4.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax4.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax4.set_xlabel(x)
    ax4.set_ylabel(y)
    grid_format(ax4)

    ax5 = plt.subplot(335)
    x = 'zeta'
    y = 'chi'
    ax5.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax5.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax5.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax5.set_xlabel(x)
    ax5.set_ylabel(y)
    grid_format(ax5)

    ax6 = plt.subplot(336)
    x = 'delta'
    y = 'chi'
    ax6.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax6.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax6.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax6.set_xlabel(x)
    ax6.set_ylabel(y)
    grid_format(ax6)

    ax7 = plt.subplot(337)
    x = 'zeta'
    y = 'zeta'
    ax7.scatter(n0_all_angles[angles[x]][
                :-1], n0_all_angles[angles[y]][1:], c='b', marker='x', s=80)
    ax7.scatter(nm_all_angles[angles[x]][
                :-1], nm_all_angles[angles[y]][1:], c='r', marker='o', s=25)
    ax7.scatter(m__all_angles[angles[x]][
                :-1], m__all_angles[angles[y]][1:], c='g', marker='s', s=25)
    ax7.set_xlabel(x)
    ax7.set_ylabel(y + ' + 1')
    grid_format(ax7)

    ax8 = plt.subplot(338)
    x = 'epsilon'
    y = 'epsilon'
    ax8.scatter(n0_all_angles[angles[x]][:-1], n0_all_angles[angles[y]][1:],
                c='b', marker='x', s=80, label="starting")
    ax8.scatter(nm_all_angles[angles[x]][:-1], nm_all_angles[angles[y]][1:],
                c='r', marker='o', s=25, label="final")
    ax8.scatter(m__all_angles[angles[x]][:-1], m__all_angles[angles[y]][1:],
                c='g', marker='s', s=25, label="minimized")
    ax8.set_xlabel(x)
    ax8.set_ylabel(y + ' + 1')
    grid_format(ax8)
    ax8.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)

    ax9 = plt.subplot(339)
    x = 'zeta'
    y = 'P'
    ax9.scatter(n0_all_angles[angles[x]][:-1], n0_all_angles[angles[y]][1:],
                c='b', marker='x', s=80)
    ax9.scatter(nm_all_angles[angles[x]][:-1], nm_all_angles[angles[y]][1:],
                c='r', marker='o', s=25)
    ax9.scatter(m__all_angles[angles[x]][:-1], m__all_angles[angles[y]][1:],
                c='g', marker='s', s=25)
    ax9.set_xlim([0, 360])
    ax9.set_ylim([0, 360])
    ax9.set_xlabel(x)
    ax9.set_ylabel(y)

    plt.show()


def plot3d_dihedral(min_angles, raw_angles):
    from mayavi import mlab

    logging.debug('testing')
    (n_thetas, n_frames, n_bps, n_angles) = min_angles.shape
    (n_thetas_r, n_frames_r, n_bps_r, n_angles_r) = raw_angles.shape


def plot_dna_min_dihedral(min_angles, raw_angles):
    logging.debug('testing')
    plt.ion()
    (n_thetas, n_frames, n_bps, n_angles) = min_angles.shape
    (n_thetas_r, n_frames_r, n_bps_r, n_angles_r) = raw_angles.shape

    # assert that raw_angle parameters match the min parameters

    angles = {r'$\alpha$': 1, r'$\beta$': 2, r'$\gamma$': 3, r'$\delta$': 4,
              r'$\epsilon$': 5, r'$\zeta$': 6, r'$\chi$': 7}
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    step = 1
    frame_0 = step - 1

    for frame in xrange(frame_0, n_frames, step):
        # for frame in xrange(n_frames):
        plt.figure()

        plt.subplot(331)
        x = r'$\zeta$'
        y = r'$\alpha$'
        for tm in xrange(n_thetas):
            # only need to do this once
            zero_to_threesixty(min_angles[tm, frame])
            zero_to_threesixty(raw_angles[tm, frame])
            plt.plot(raw_angles[tm, frame][0:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')

        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelbottom='off') #, labeltop='on')

        plt.subplot(332)
        x = r'$\zeta$'
        y = r'$\beta$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][0:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off') #, labeltop='on')

        plt.subplot(333)
        x = r'$\zeta$'
        y = r'$\epsilon$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(
            labelleft='off', labelbottom='off') #, labeltop='on', labelright='on')

        plt.subplot(334)
        x = r'$\alpha$'
        y = r'$\gamma$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelbottom='off')

        plt.subplot(335)
        x = r'$\zeta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off')

        plt.subplot(336)
        x = r'$\delta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off') #, labelright='on')

        plt.subplot(337)
        x = r'$\zeta$'
        y = r'$\zeta$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()

        plt.subplot(338)
        x = r'$\epsilon$'
        y = r'$\epsilon$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off')
        # plt.legend(['raw', 'min'], loc='upper left',
        # bbox_to_anchor=[1.1475, 1.07], numpoints=1)

        plt.subplot(339)
        x = r'$\zeta$'
        y = r'$\delta$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off') #, labelright='on')
        # plt.xlim([-999, 999])
        # plt.ylim([-999, 999])
        plt.legend(['raw', 'min'], loc='upper left',
                   bbox_to_anchor=[1.1475, 1.07], numpoints=1)

        plt.suptitle('%d Steps: scatter plots of selected torsional angles' %
                     ((frame + 1) * 100))
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.pdf' %
                    ((frame + 1) * 100), bbox_inches='tight')
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.png' %
                    ((frame + 1) * 100), bbox_inches='tight')
        plt.draw()

    plt.show()


def plot_dna_dihedral(min_angles):
    logging.debug('testing')
    plt.ion()
    (n_thetas, n_frames, n_bps, n_angles) = min_angles.shape

    # assert that raw_angle parameters match the min parameters

    angles = {r'$\alpha$': 1, r'$\beta$': 2, r'$\gamma$': 3, r'$\delta$': 4,
              r'$\epsilon$': 5, r'$\zeta$': 6, r'$\chi$': 7}
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    step = 1
    frame_0 = step - 1

    for frame in xrange(frame_0, n_frames, step):
        # for frame in xrange(n_frames):
        plt.figure()

        plt.subplot(331)
        x = r'$\zeta$'
        y = r'$\alpha$'
        for tm in xrange(n_thetas):
            # only need to do this once
            zero_to_threesixty(min_angles[tm, frame])
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')

        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelbottom='off')  #, labeltop='on')

        plt.subplot(332)
        x = r'$\zeta$'
        y = r'$\beta$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off')  #, labeltop='on')

        plt.subplot(333)
        x = r'$\zeta$'
        y = r'$\epsilon$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(
            labelleft='off', labelbottom='off') #, labeltop='on', labelright='on')

        plt.subplot(334)
        x = r'$\alpha$'
        y = r'$\gamma$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelbottom='off')

        plt.subplot(335)
        x = r'$\zeta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off')

        plt.subplot(336)
        x = r'$\delta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off') #, labelright='on')

        plt.subplot(337)
        x = r'$\zeta$'
        y = r'$\zeta$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()

        plt.subplot(338)
        x = r'$\epsilon$'
        y = r'$\epsilon$'
        lines = []
        labels = []
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off')
        plt.legend(['min', 'raw'], loc='upper left',
                   bbox_to_anchor=[1.1475, 1.07], numpoints=1)

        if False:
            plt.subplot(339)
            x = r'$\zeta$'
            y = r'$\delta$'
            lines = []
            labels = []
            for tm in xrange(n_thetas):
                plt.plot(all_angles[tm, frame][:, angles[x]],
                         all_angles[tm, frame][:, angles[y]], 'gx')
                labels.append('run %d' % theta_max[tm])
            plt.tick_params(labelleft='off') #, labelright='on')
            plt.xlabel(x)
            plt.ylabel(y)
            plt_format()
            plt.xlim([-999, 999])
            plt.ylim([-999, 999])
            plt.legend(labels)

        plt.suptitle('%d Steps: scatter plots of selected torsional angles' %
                     ((frame + 1) * 10))
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.pdf' %
                    ((frame + 1) * 10), bbox_inches='tight')
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.png' %
                    ((frame + 1) * 10), bbox_inches='tight')
        plt.draw()

    plt.show()


def plt_format():

    plt.xlim([0, 360])
    plt.ylim([0, 360])
    plt.xticks([0, 60, 120, 180, 240, 300, 360])
    plt.yticks([0, 60, 120, 180, 240, 300, 360])


def grid_format(ax):

    ax.set_xlim([0, 360])
    ax.set_ylim([0, 360])
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_yticks([0, 60, 120, 180, 240, 300, 360])


def zero_to_threesixty(all_angles):
    for (i, angle_list) in enumerate(all_angles):
        zero_to_threesixty_list(angle_list)


def zero_to_threesixty_list(angle_list):
    for (j, angle) in enumerate(angle_list):
        if angle < 0:
            angle_list[j] = angle + 360


def load_angles(file_name):
    i_frame = 0
    all_frames = []
    with open(file_name, 'r') as in_file:
        for line in in_file.readlines():
            if line.find("#", 0) == 0:
                if line.find("frame") > 0:
                    # start a new frame
                    i_frame += 1
                    if i_frame > 1:
                        all_frames.append(frame)
                    frame = []
                    continue
                continue
            frame.append(numpy.array(line.split('\t')[2:], dtype='float64'))

    all_frames.append(frame)
    return numpy.array(all_frames)


def read_flex_resids(flex_file):
    '''
    Read flexible DNA resids from the file created using 'write_flex_resids'
    First line of file should have the DNA chain names
    Following lines should have the resids from each chain that are designated
    as flexible, for example:
    A B
    1 4
    2 3
    3 2
    4 1
    '''
    lines = [line.strip() for line in open(flex_file)]
    segnames = []
    segnames.append(lines[0].split()[0])
    segnames.append(lines[0].split()[1])
    flex_resids = numpy.genfromtxt(flex_file, dtype='int', delimiter=" ")[1:]
    return (segnames, flex_resids)



def get_angles_df_from_mask(dcd_file_name, pdb_file_name, first_last_resids,
                            flex_file, nf=0):
    """
    This is not complete or tested.  It would be nice to have
    because the ordering of DNA atoms in PDB files is not consistent.
    Additionally, the atoms included in the last base on each chain
    is not consistent.
    """
    assert op.exists(pdb_file_name), 'No such file: %s' % pdb_file_name
    assert op.exists(dcd_file_name), 'No such file: %s' % dcd_file_name
    txtOutput = []
    molecule_type = "dna"

    # read flex file
    # flex_file = pdb_file_name[:-3] + 'flex'
    dna_segnames, flex_resids = read_flex_resids(flex_file)
    assert len(dna_segnames) == 2, 'More than 2 DNA segnames in %s' % flex_file
    flex_resids = [flex_resids[:, 0], flex_resids[:, 1][::-1]]

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file_name)

    dna1 = sasmol.SasMol(0)
    dna1_filter = ' segname[i] == "%s" ' % dna_segnames[0]
    print 'DNA segname 1: %s' % dna_segnames[0]
    error, dna1_mask = mol.get_subset_mask(dna1_filter)
    assert not error, error
    error = mol.copy_molecule_using_mask(dna1, dna1_mask, 0)
    assert not error, error

    dna2 = sasmol.SasMol(0)
    dna2_filter = ' segname[i] == "%s" ' % dna_segnames[1]
    print 'DNA segname 2: %s' % dna_segnames[1]
    error, dna2_mask = mol.get_subset_mask(dna2_filter)
    assert not error, error
    error = mol.copy_molecule_using_mask(dna2, dna2_mask, 0)
    assert not error, error

    masks = [dna1_mask, dna2_mask]
    dna_mols = [dna1, dna2]
    segnames = dna_segnames

    dcdfile = mol.open_dcd_read(dcd_file_name)
    if nf == 0:
        nf = dcdfile[2]  # number of frames

    bp_mols = [[None] * nf, [None] * nf]
    bp_masks = [[None] * nf, [None] * nf]
    angles = []
    for i in xrange(nf):
        frame = []
        segname = []
        resid = []
        resname = []
        all_alpha = []
        all_beta = []
        all_gamma = []
        all_delta = []     # <--- significant
        all_epsilon = []
        all_zeta = []
        all_chi = []       # <--- significant

        mol.read_dcd_step(dcdfile, i)

        for (j, dna_strand) in enumerate(dna_mols):
            first_last_resid = first_last_resids[j]
            numranges = 1
            res_low = [flex_resids[j][0]]
            n_cont = [len(flex_resids[j][:])]
            res_high = [res_low[0] + n_cont[0]]
            assert numpy.alltrue(flex_resids[j] == numpy.arange(
                res_low[0], res_high[0])), (
                    'ERROR: update code to handle discontinuous flexible '
                    'regions or break input into continuous flexible regions')

            flexible_residues = dihedral_angles.get_flexible_residues(
                numranges, res_low, n_cont)

            base_indices, base_masks = (
                dihedral_angles.get_rotation_indices(
                    dna_strand, molecule_type, flexible_residues, txtOutput))

            error, coor = mol.get_coor_using_mask(0, masks[j])
            assert not error, error

            frame += [i + 1] * n_cont[0]
            resid += range(res_low[0], res_low[0] + n_cont[0])
            segname += [segnames[j]] * n_cont[0]

            for res in xrange(n_cont[0]):
                bp_mol = bp_mols[j][res]
                bp_mask = bp_masks[j][res]
                if not bp_mol:
                    bp_mol = bp_mols[j][res] = sasmol.SasMol(0)
                    bp_filter = 'resid %d' % resid
                    bp_filter_python = basis_to_python.parse_basis(bp_filter)
                    error, bp_mask = dna_strand.get_subset_mask(
                        bp_filter_python)
                    assert not error, error
                    bp_masks[j][res] = bp_mask
                    error = dna_strand.copy_mol_using_mask(bp_mol, bp_mask, 0)
                    assert not error, error

                q0 = res_low[0] + res
                resname.append(dna_strand.resname()[numpy.where(
                    dna_strand.resid() == q0)[0][0]][0])

                # get the indices and mask for this residue
                indices = base_indices[q0]
                this_mask = numpy.array(base_masks[q0])

                # alpha
                alpha_angle = dihedral_angles.measure(
                    coor, indices, "alpha", this_mask, q0,
                    first_last_resid, molecule_type)
                str_alpha_angle = '%.1f' % alpha_angle
                all_alpha.append(alpha_angle)
                logging.debug('alpha: ' + str_alpha_angle)

                # beta
                beta_angle = dihedral_angles.measure(
                    coor, indices, "beta", this_mask, q0,
                    first_last_resid, molecule_type)
                str_beta_angle = '%.1f' % (beta_angle)
                all_beta.append(beta_angle)
                logging.debug('beta: ' + str_beta_angle)

                # gamma
                gamma_angle = dihedral_angles.measure(
                    coor, indices, "gamma", this_mask, q0,
                    first_last_resid, molecule_type)
                str_gamma_angle = '%.1f' % (gamma_angle)
                all_gamma.append(gamma_angle)
                logging.debug('gamma: ' + str_gamma_angle)

                # delta
                delta_angle = dihedral_angles.measure(
                    coor, indices, "delta", this_mask, q0,
                    first_last_resid, molecule_type)
                str_delta_angle = '%.1f' % (delta_angle)
                all_delta.append(delta_angle)
                logging.debug('delta: ' + str_delta_angle)

                # epsilon
                epsilon_angle = dihedral_angles.measure(
                    coor, indices, "epsilon", this_mask, q0,
                    first_last_resid, molecule_type)
                str_epsilon_angle = '%.1f' % (epsilon_angle)
                all_epsilon.append(epsilon_angle)
                logging.debug('epsilon: ' + str_epsilon_angle)

                # zeta
                zeta_angle = dihedral_angles.measure(
                    coor, indices, "zeta", this_mask, q0,
                    first_last_resid, molecule_type)
                str_zeta_angle = '%.1f' % (zeta_angle)
                all_zeta.append(zeta_angle)
                logging.debug('zeta: ' + str_zeta_angle)

                # chi
                chi_angle = dihedral_angles.measure(
                    coor, indices, "chi", this_mask, q0,
                    first_last_resid, molecule_type)
                str_chi_angle = '%.1f' % (chi_angle)
                all_chi.append(chi_angle)
                logging.debug('alpha: ' + str_alpha_angle)

        all_angles = {'frame': frame, 'segname': segname, 'resid': resid,
                      'resname': resname, 'alpha': all_alpha,
                      'beta': all_beta, 'gamma': all_gamma,
                      'delta': all_delta, 'epsilon': all_epsilon,
                      'zeta': all_zeta, 'chi': all_chi}
        for key in ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                    'chi']:
            zero_to_threesixty_list(all_angles[key])

        column_labels = ('frame', 'segname', 'resid', 'resname', 'alpha',
                         'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                         'chi')
        angles.append(pd.DataFrame(all_angles, columns=column_labels))
        print ('\nfinished calculating %d dihedral angles from frame %d '
               'of %s\n' % (7 * n_cont[0], i + 1, dcd_file_name))
    df = pd.concat(angles)
    df.to_hdf(dcd_file_name[:-3] + 'hdf', 'angles')
    return df


def get_angles_df(dcd_file_name, pdb_file_name, first_last_resids, flex_file,
                  nf=0, drude=False):
    if drude:
        print 'Using drude DNA atom ordering'
    assert op.exists(pdb_file_name), 'No such file: %s' % pdb_file_name
    assert op.exists(dcd_file_name), 'No such file: %s' % dcd_file_name
    txtOutput = []
    molecule_type = "dna"

    # read flex file
    # flex_file = pdb_file_name[:-3] + 'flex'
    dna_segnames, flex_resids = read_flex_resids(flex_file)
    assert len(dna_segnames) == 2, 'More than 2 DNA segnames in %s' % flex_file
    flex_resids = [flex_resids[:, 0], flex_resids[:, 1][::-1]]

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file_name)

    dna1 = sasmol.SasMol(0)
    dna1_filter = ' segname[i] == "%s" ' % dna_segnames[0]
    print 'DNA segname 1: %s' % dna_segnames[0]
    error, dna1_mask = mol.get_subset_mask(dna1_filter)
    assert not error, error
    error = mol.copy_molecule_using_mask(dna1, dna1_mask, 0)
    assert not error, error

    dna2 = sasmol.SasMol(0)
    dna2_filter = ' segname[i] == "%s" ' % dna_segnames[1]
    print 'DNA segname 2: %s' % dna_segnames[1]
    error, dna2_mask = mol.get_subset_mask(dna2_filter)
    assert not error, error
    error = mol.copy_molecule_using_mask(dna2, dna2_mask, 0)
    assert not error, error

    masks = [dna1_mask, dna2_mask]
    dna_mols = [dna1, dna2]
    segnames = dna_segnames

    dcdfile = mol.open_dcd_read(dcd_file_name)
    if nf == 0:
        nf = dcdfile[2]  # number of frames

    angles = []
    for i in xrange(nf):
        frame = []
        segname = []
        resid = []
        resname = []
        all_alpha = []
        all_beta = []
        all_gamma = []
        all_delta = []     # <--- significant
        all_epsilon = []
        all_zeta = []
        all_chi = []       # <--- significant

        mol.read_dcd_step(dcdfile, i)

        for (j, dna) in enumerate(dna_mols):
            first_last_resid = first_last_resids[j]
            numranges = 1
            res_low = [flex_resids[j][0]]
            n_cont = [len(flex_resids[j][:])]
            res_high = [res_low[0] + n_cont[0]]
            assert numpy.alltrue(flex_resids[j] == numpy.arange(
                res_low[0], res_high[0])), (
                'ERROR: update code to handle discontinuous flexible '
                'regions or break input into continuous flexible regions')

            flexible_residues = dihedral_angles.get_flexible_residues(
                numranges, res_low, n_cont)

            base_indices, base_masks = (
                dihedral_angles.get_rotation_indices(
                    dna, molecule_type, flexible_residues, txtOutput))

            error, coor = mol.get_coor_using_mask(0, masks[j])
            assert len(error) < 1, error

            frame += [i + 1] * n_cont[0]
            resid += range(res_low[0], res_low[0] + n_cont[0])
            segname += [segnames[j]] * n_cont[0]

            for res in xrange(n_cont[0]):
                q0 = res_low[0] + res
                resname.append(dna.resname()[numpy.where(
                    dna.resid() == q0)[0][0]][0])

                # get the indices and mask for this residue
                indices = base_indices[q0]
                this_mask = numpy.array(base_masks[q0])

                # alpha
                alpha_angle = dihedral_angles.measure(
                    coor, indices, "alpha", this_mask, q0,
                    first_last_resid, molecule_type, drude)
                str_alpha_angle = '%.1f' % alpha_angle
                all_alpha.append(alpha_angle)
                logging.debug('alpha: ' + str_alpha_angle)

                # beta
                beta_angle = dihedral_angles.measure(
                    coor, indices, "beta", this_mask, q0,
                    first_last_resid, molecule_type, drude)
                str_beta_angle = '%.1f' % (beta_angle)
                all_beta.append(beta_angle)
                logging.debug('beta: ' + str_beta_angle)

                # gamma
                gamma_angle = dihedral_angles.measure(
                    coor, indices, "gamma", this_mask, q0,
                    first_last_resid, molecule_type, drude)
                str_gamma_angle = '%.1f' % (gamma_angle)
                all_gamma.append(gamma_angle)
                logging.debug('gamma: ' + str_gamma_angle)

                # delta
                delta_angle = dihedral_angles.measure(
                    coor, indices, "delta", this_mask, q0,
                    first_last_resid, molecule_type, drude)
                str_delta_angle = '%.1f' % (delta_angle)
                all_delta.append(delta_angle)
                logging.debug('delta: ' + str_delta_angle)

                # epsilon
                epsilon_angle = dihedral_angles.measure(
                    coor, indices, "epsilon", this_mask, q0,
                    first_last_resid, molecule_type, drude)
                str_epsilon_angle = '%.1f' % (epsilon_angle)
                all_epsilon.append(epsilon_angle)
                logging.debug('epsilon: ' + str_epsilon_angle)

                # zeta
                zeta_angle = dihedral_angles.measure(
                    coor, indices, "zeta", this_mask, q0,
                    first_last_resid, molecule_type, drude)
                str_zeta_angle = '%.1f' % (zeta_angle)
                all_zeta.append(zeta_angle)
                logging.debug('zeta: ' + str_zeta_angle)

                # chi
                chi_angle = dihedral_angles.measure(
                    coor, indices, "chi", this_mask, q0,
                    first_last_resid, molecule_type, drude)
                str_chi_angle = '%.1f' % (chi_angle)
                all_chi.append(chi_angle)
                logging.debug('alpha: ' + str_alpha_angle)

        all_angles = {'frame': frame, 'segname': segname, 'resid': resid,
                      'resname': resname, 'alpha': all_alpha,
                      'beta': all_beta, 'gamma': all_gamma,
                      'delta': all_delta, 'epsilon': all_epsilon,
                      'zeta': all_zeta, 'chi': all_chi}
        for key in ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                    'chi']:
            zero_to_threesixty_list(all_angles[key])

        column_labels = ('frame', 'segname', 'resid', 'resname', 'alpha',
                         'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                         'chi')
        angles.append(pd.DataFrame(all_angles, columns=column_labels))
        print ('\nfinished calculating %d dihedral angles from frame %d '
               'of %s\n' % (7 * n_cont[0], i + 1, dcd_file_name))
    df = pd.concat(angles)
    df.to_hdf(dcd_file_name[:-3] + 'hdf', 'angles')
    return df


def get_angles_file(dcd_file_name, pdb_file_name, first_last_resids, flex_file,
                    nf=0, drude=False):

    txtOutput = []
    all_alpha = []
    all_beta = []
    all_gamma = []
    all_delta = []     # <--- significant
    all_epsilon = []
    all_zeta = []
    all_chi = []       # <--- significant

    molecule_type = "dna"

    # read flex file
    # flex_file = pdb_file_name[:-3] + 'flex'
    dna_segnames, flex_resids = read_flex_resids(flex_file)
    assert len(dna_segnames) == 2, 'More than 2 DNA segnames in %s' % flex_file
    flex_resids = [flex_resids[:, 0], flex_resids[:, 1][::-1]]

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file_name)

    dna1 = sasmol.SasMol(0)
    dna1_filter = ' segname[i] == "%s" ' % dna_segnames[0]
    print 'DNA segname 1: %s' % dna_segnames[0]
    error, dna1_mask = mol.get_subset_mask(dna1_filter)
    assert not error, error
    error = mol.copy_molecule_using_mask(dna1, dna1_mask, 0)
    assert not error, error

    dna2 = sasmol.SasMol(0)
    dna2_filter = ' segname[i] == "%s" ' % dna_segnames[1]
    print 'DNA segname 2: %s' % dna_segnames[1]
    error, dna2_mask = mol.get_subset_mask(dna2_filter)
    assert not error, error
    error = mol.copy_molecule_using_mask(dna2, dna2_mask, 0)
    assert not error, error

    masks = [dna1_mask, dna2_mask]
    dna_mols = [dna1, dna2]
    segnames = dna_segnames
    # res_mol = sasmol.SasMol(0)

    st1 = "# segname\tbase\talpha\tbeta\tgamma\tdelta\tepsilon\tzeta\tchi\n"
    with open(dcd_file_name[:-3] + 'ddat', 'w') as out_file:
        out_file.write("%s" % (st1))

        dcdfile = mol.open_dcd_read(dcd_file_name)
        if nf == 0:
            nf = dcdfile[2]  # number of frames

        for i in xrange(nf):

            mol.read_dcd_step(dcdfile, i)
            out_file.write("# frame %d\n" % (i + 1))

            for (j, dna) in enumerate(dna_mols):
                first_last_resid = first_last_resids[j]
                numranges = 1
                res_low = [flex_resids[j][0]]
                n_cont = [len(flex_resids[j][:])]
                res_high = [res_low[0] + n_cont[0]]
                assert numpy.alltrue(flex_resids[j] == numpy.arange(
                    res_low[0], res_high[0])), (
                    'ERROR: update code to handle discontinuous flexible '
                    'regions or break input into continuous flexible regions')

                flexible_residues = dihedral_angles.get_flexible_residues(
                    numranges, res_low, n_cont)

                base_indices, base_masks = (
                    dihedral_angles.get_rotation_indices(
                        dna, molecule_type, flexible_residues, txtOutput))
                # print "residue_rotation_indices = ",residue_rotation_indices

                # reslist = range(first_last_resid[0],first_last_resid[1]+1)
                # print 'reslist = ',reslist

                error, coor = mol.get_coor_using_mask(0, masks[j])
                assert len(error) < 1, error
                # dna.setCoor(coor) # not necessary
                # coor = dna.coor()

                for res in xrange(n_cont[0]):
                    q0 = res_low[0] + res
                    resname = dna.resname()[numpy.where(
                        dna.resid() == q0)[0][0]][0]
                    st = '%s\t%s\t%d' % (segnames[j], resname, q0)

                    # get the indices and mask for this residue
                    indices = base_indices[q0]
                    this_mask = numpy.array(base_masks[q0])

                    # alpha
                    alpha_angle = dihedral_angles.measure(
                        coor, indices, "alpha", this_mask, q0,
                        first_last_resid, molecule_type, drude)
                    str_alpha_angle = '%.1f' % alpha_angle
                    all_alpha.append(alpha_angle)
                    st = st + '\t' + str_alpha_angle
                    logging.debug('alpha: ' + str_alpha_angle)

                    # beta
                    beta_angle = dihedral_angles.measure(
                        coor, indices, "beta", this_mask, q0,
                        first_last_resid, molecule_type, drude)
                    str_beta_angle = '%.1f' % (beta_angle)
                    all_beta.append(beta_angle)
                    st = st + '\t' + str_beta_angle
                    logging.debug('beta: ' + str_beta_angle)

                    # gamma
                    gamma_angle = dihedral_angles.measure(
                        coor, indices, "gamma", this_mask, q0,
                        first_last_resid, molecule_type, drude)
                    str_gamma_angle = '%.1f' % (gamma_angle)
                    all_gamma.append(gamma_angle)
                    st = st + '\t' + str_gamma_angle
                    logging.debug('gamma: ' + str_gamma_angle)

                    # delta
                    delta_angle = dihedral_angles.measure(
                        coor, indices, "delta", this_mask, q0,
                        first_last_resid, molecule_type, drude)
                    str_delta_angle = '%.1f' % (delta_angle)
                    all_delta.append(delta_angle)
                    st = st + '\t' + str_delta_angle
                    logging.debug('delta: ' + str_delta_angle)

                    # epsilon
                    epsilon_angle = dihedral_angles.measure(
                        coor, indices, "epsilon", this_mask, q0,
                        first_last_resid, molecule_type, drude)
                    str_epsilon_angle = '%.1f' % (epsilon_angle)
                    all_epsilon.append(epsilon_angle)
                    st = st + '\t' + str_epsilon_angle
                    logging.debug('epsilon: ' + str_epsilon_angle)

                    # zeta
                    zeta_angle = dihedral_angles.measure(
                        coor, indices, "zeta", this_mask, q0,
                        first_last_resid, molecule_type, drude)
                    str_zeta_angle = '%.1f' % (zeta_angle)
                    all_zeta.append(zeta_angle)
                    st = st + '\t' + str_zeta_angle
                    logging.debug('zeta: ' + str_zeta_angle)

                    # chi
                    chi_angle = dihedral_angles.measure(
                        coor, indices, "chi", this_mask, q0,
                        first_last_resid, molecule_type, drude)
                    str_chi_angle = '%.1f' % (chi_angle)
                    all_chi.append(chi_angle)
                    st = st + '\t' + str_chi_angle
                    logging.debug('alpha: ' + str_alpha_angle)

                    out_file.write("%s\n" % (st))

            all_angles = [all_alpha, all_beta, all_gamma, all_delta,
                          all_epsilon, all_zeta, all_chi]

    print ('\nfinished calculating %d dihedral angles from frame %d of %s\n' %
           (7 * n_cont[0], i + 1, dcd_file_name))

    return all_angles


def main():
    NotImplemented

if __name__ == '__main__':

    if '-v' in sys.argv:
        logging.basicConfig(level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

        main()
