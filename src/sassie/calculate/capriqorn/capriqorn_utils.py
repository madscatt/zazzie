#
# Hailiang Zhang
# NIST & UTK
#

import sys
import os
import locale
import glob
import shutil
import h5py
import random
import numpy as np
import re

import sasmol.sasmol as sasmol

def process_preprocessor_yaml(file_yaml, mvars, scvars):
    '''
    '''
    pdbfile = mvars.pdbfile
    dcdfile = mvars.dcdfile
    if scvars.intype == 'dcd':
        xstfile = mvars.xstfile
    aliasfile = mvars.aliasfile

    # print 'file_yaml = ',file_yaml
    # print 'core_path = ',mvars.core_path
    # print 'core_path plus file_yaml =
    # ',os.path.join(mvars.core_path,file_yaml)

    # file_yaml = os.path.join(mvars.core_path,file_yaml)

    fo = open('tmp.txt', 'w')
    lines = open(file_yaml).readlines()
    for line in lines[7:]:  # skip the first 7 lines since they are not used
        words = line.split()
        if words[0] == "pdbName:":
            line = '    pdbName: ' + pdbfile + '\n'
        if words[0] == "trajName:":
            line = '    trajName: ' + dcdfile + '\n'
            if scvars.intype == 'dcd':
                line += '    xstName: ' + xstfile + '\n'
        if words[0] == "aliasName:":
            line = '    aliasName: ' + aliasfile + '\n'
        if words[0] == "last:":
            line = '    last: %d\n' % scvars.number_of_frames
        fo.write(line)
    fo.close()
    shutil.move('tmp.txt', file_yaml)


def process_hist_yaml(file_yaml):
    '''
    '''
    # print 'file_yaml = ',file_yaml
    # file_yaml = os.path.join(mvars.core_path,file_yaml)

    # print 'file_yaml = ',file_yaml
    fo = open('tmp.txt', 'w')
    for line in open(file_yaml).readlines():
        if (line == '  gpus: 2\n'):
            line = '  gpus: 0\n'
        fo.write(line)
    fo.close()
    shutil.move('tmp.txt', file_yaml)


def process_postprocessor_yaml(file_yaml, mvars, scvars):
    '''
    '''
    pdbfile = mvars.pdbfile
    dcdfile = mvars.dcdfile
    if scvars.intype == 'dcd':
        xstfile = mvars.xstfile

    # print 'file_yaml = ',file_yaml
    # print 'core_path = ',mvars.core_path
    # print 'core_path plus file_yaml =
    # ',os.path.join(mvars.core_path,file_yaml)

    # file_yaml = os.path.join(mvars.core_path,file_yaml)

    fo = open('tmp.txt', 'w')
    lines = open(file_yaml).readlines()
    for line in lines:
        words = line.split()
        if words[0] == "nq:":
            line = '    nq: %d\n' % mvars.number_q_values
        if words[0] == "dq:":
            line = '    dq: %f\n' % (mvars.q_max / (mvars.number_q_values - 1))
        fo.write(line)
    fo.close()
    shutil.move('tmp.txt', file_yaml)


def generate_yaml(output_folder, mvars, scvars):
    '''
    generate_preprocessor_input
    generate_histograms_input
    generate_postprocessor_input,
    '''
    os.chdir(output_folder)

    exe_path = os.path.dirname(sys.executable)
    os.system(os.path.join(exe_path, 'generate_preprocessor_input'))

    process_preprocessor_yaml('preprocessor.yaml', mvars, scvars)
    os.system(os.path.join(exe_path, 'generate_histograms_input'))

    process_hist_yaml('histograms.yaml')

    os.system(os.path.join(exe_path, 'generate_postprocessor_input'))
    process_postprocessor_yaml('postprocessor.yaml', mvars, scvars)
    os.chdir(os.path.join("..", "..", ".."))

def prepare_capriqorn_inputs(mol, mvars, scvars, output_folder):
    '''
    Functionality: deuterate the exchangeable H's based on the HD exhcange ratio in various ways
    Input: mol--sasmol object
           mvars--module variables
    Return: None
    '''
    # folder preparation
    output_dir = os.path.join(output_folder, 'yaml')
    os.mkdir(output_dir)

    # generate yaml files
    generate_yaml(output_dir, mvars, scvars)

    return


def generate_ascii_intensity(h5file, asciifile):
    '''
    geneate scii intensity file from a h5 file
    '''
    f = h5py.File(h5file)
    # field = f['10/intensity/I_solv']
    field = f['10/intensity/dI']
    fo = open(asciifile, 'w')
    for i in field:
        fo.write('%8.6f %12.5f %12.5f\n' % (i[0], i[1], i[1] / 10.))
    fo.close()

def get_crdbox_number_of_frames(file):
    '''
    get number of frames from crdbox
    '''
    lines = open(file).readlines()
    natoms = locale.atoi(
        re.compile(r'TITLE : Created by VMD with (\d+) atoms\n').findall(lines[0])[0])
    n = 0
    for line in lines[1:]:
        n += len(line.split())
    nframes = n / (3 * (natoms + 1))  # NOTE to ZHL: is this robust?
    return nframes

def status_update(self, fraction_done):

    pgui = self.run_utils.print_gui   
    report_string = 'STATUS\t%f' % fraction_done
    pgui(report_string)

    return

def calculate(self, mvars, scvars, output_folder):

    output_dir = os.path.join(output_folder, 'yaml')
    os.system("cp %s %s" % (mvars.pdbfile, output_dir))
    os.system("cp %s %s" % (mvars.dcdfile, output_dir))
    if scvars.intype == 'dcd':
        os.system("cp %s %s" % (mvars.xstfile, output_dir))
    os.system("cp %s %s" % (mvars.aliasfile, output_dir))
    os.chdir(output_dir)
    exe_path = os.path.dirname(sys.executable)

    status_update(self, 0.2)
    os.system(os.path.join(exe_path, 'preprocessor'))
    status_update(self, 0.4)
    os.system(os.path.join(exe_path, 'histograms_par'))
    status_update(self, 0.7)
    os.system(os.path.join(exe_path, 'postprocessor'))
    status_update(self, 0.8)
    os.system('cp %s ..' %
              os.path.join('postprocessor_output', 'intensities.h5'))
    os.chdir(os.path.join(".."))
    generate_ascii_intensity("intensities.h5", "Iq.txt")
    os.chdir(os.path.join("..", ".."))
    status_update(self, 0.9)
