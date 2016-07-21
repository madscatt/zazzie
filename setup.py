'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os,sys
from distutils.core import setup
from distutils      import sysconfig
#from distutils.extension import Extension
from numpy.distutils.core import Extension, setup

#       SETUP
#
#       12/01/2009      --      initial coding              :       jc
#       11/22/2014      --      adapted for 2.0             :       jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
      setup.py is the script to install and/or update sassie.

	installation: use the installer.py script (which calls this script)

	update: using the python version that was used to install python type:

	> sudo python setup.py build
	> sudo python setup.py install

'''

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

import sys
sys.path.append(os.path.join('sassie','util'))
import sasconfig as sasconfig

print 'sasconfig.arch = ',sasconfig.__arch__
sasconfig.arch = sasconfig.__arch__

if sasconfig.__arch__ == "cluster":
    os.environ["CC"] = "/share/apps/local/bin/gcc"
    os.environ["CXX"] = "/share/apps/local/bin/g++"

core_libraries_include = sasconfig.__core_libraries_include__
core_libraries_lib = sasconfig.__core_libraries_lib__
cuda_flag = sasconfig.__cuda__

#if CUDA and cuda_driver:
#    include_dir_names.append(os.path.join(cuda_dir,'include'))
#    library_dir_names.append(os.path.join(cuda_dir,'lib64'))
#    library_names.append('cuda')
#    library_names.append('cudart')
#    macros.append(('CUDA_DRIVER','1'))
#if CUDA and cuda_lib:
#    include_dir_names.append(os.path.join(core_libraries,'include'))
#    library_dir_names.append(os.path.join(core_libraries,'lib'))
#    library_names.append(cuda_library_name)
#    macros.append(('CUDA_LIB','1'))

if not cuda_flag:
    overlap_macros = [('CPP_LIB','1')]
    overlap_libraries = ["overlap"]
else:
    cuda_dir = sasconfig.__cuda_path__
    overlap_macros = []
    overlap_libraries = []
    cuda_libraries = []
    
    overlap_libraries.append('cuda')
    overlap_libraries.append('cudart')
    overlap_libraries.append('cudaOverlap')
    #cuda_libraries.append('cuda')
    #cuda_libraries.append('cudart')

    if sasconfig.__arch__ == "cluster":
    #    core_libraries_include.append(os.path.join(cuda_dir,'targets/x86_64-linux/include'))
    #    core_libraries_lib.append(os.path.join(cuda_dir,'targets/x86_64-linux/lib'))
        core_libraries_include.append(os.path.join(cuda_dir,'include'))
        core_libraries_lib.append(os.path.join(cuda_dir,'lib64'))
    elif sasconfig.__arch__ == "linux":
        core_libraries_include.append(os.path.join(cuda_dir,'include'))
        core_libraries_lib.append(os.path.join(cuda_dir,'lib64'))
    elif sasconfig.__arch__ == "mac":
        core_libraries_include.append(os.path.join(cuda_dir,'include'))
        core_libraries_lib.append(os.path.join(cuda_dir,'lib'))
    else:
        assert False, ('ERROR: sasconfig.__arch__ == %s is not a valid option\n'
                       'ERROR: please select "cluster/mac/linux"' % sasconfig.__arch__)

    overlap_macros.append(('CUDA_DRIVER','1'))
    overlap_macros.append(('CUDA_LIB','1')) 

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='sassie_2',
	version='1.99_rev_1777',
	author='Joseph E. Curtis',
	author_email='joseph.curtis@nist.gov',
	license='GPL 3',
	url='www.smallangles.net/sassie',
	platforms='Linux, Mac OS X',
	description=("A suite of programs to generate atomistic models of biological systems, calculate scattering observables, and compare results to experimental data"),
	long_description=read('README'),	
	classifiers=["Development Status :: 2.0 Release",
		"License :: OSI Approved :: GNU Public License 3",
		"Intended Audience :: Science/Research",
		"Natural Language :: English",
		"Operating System :: Linux :: MacOS :: MacOS X",
		"Programming Language :: Python :: C :: Fortran",
		"Topic :: Scientific/Engineering :: Chemistry :: Physics"],

	package_dir={'sassie':''},

	packages=['sassie','sassie/util','sassie/analyze','sassie/build','sassie/build/pdbscan','sassie/build/pdbrx','sassie/calculate','sassie/calculate/capriqorn','sassie/calculate/capriqorn/source','sassie/interface','sassie/simulate','sassie/simulate/openmm','sassie/simulate/monte_carlo','sassie/simulate/monte_carlo/extensions','sassie/simulate/monte_carlo/extensions/overlap','sassie/simulate/monte_carlo/monte_carlo_utilities','sassie/simulate/monte_carlo/monte_carlo_utilities/isopeptide_bond_torsion','sassie/simulate/monte_carlo/monte_carlo_utilities/single_stranded_nucleic_backbone_torsion','sassie/simulate/monte_carlo/monte_carlo_utilities/protein_backbone_torsion','sassie/simulate/monte_carlo/monte_carlo_utilities/double_stranded_nucleic','sassie/simulate/monte_carlo/monte_carlo_utilities/tamc_utilities','sassie/simulate/openmm','sassie/simulate/energy','sassie/simulate/energy/extensions/','sassie/simulate/energy/extensions/non_bonding','sassie/simulate/constraints','sassie/simulate/prody','sassie/tools', 'sassie/calculate/sascalc_library/', 'sassie/calculate/sascalc_library/cpp_extension'],

	ext_modules=[
	Extension('sassie.simulate.energy.non_bonding_intermolecular',['sassie/simulate/energy/extensions/non_bonding/non_bonding_intermolecular.f']),
	Extension('sassie.simulate.energy.vdw',['sassie/simulate/energy/extensions/non_bonding/vdw.f']),
    Extension('sassie.simulate.monte_carlo.pairs',['sassie/simulate/monte_carlo/extensions/pairs/pairs.f'],include_dirs=[numpy_include]),
    Extension('sassie.simulate.monte_carlo.overlap',['sassie/simulate/monte_carlo/extensions/overlap/overlap_extension.cpp'],include_dirs=core_libraries_include,library_dirs=core_libraries_lib, libraries=overlap_libraries, define_macros=overlap_macros),
    Extension('sassie.simulate.monte_carlo.ooverlap',['sassie/simulate/monte_carlo/extensions/ooverlap/ooverlap.c'],include_dirs=[numpy_include]),
    Extension('sassie.simulate.monte_carlo.vdw_overlap',['sassie/simulate/monte_carlo/extensions/vdw_overlap/vdw_overlap.f'],include_dirs=[numpy_include]),
	Extension('sassie.simulate.monte_carlo.dna_overlap',['sassie/simulate/monte_carlo/extensions/dna_overlap/dna_overlap.f']),
	Extension('sassie.simulate.energy.electrostatics',['sassie/simulate/energy/extensions/non_bonding/electrostatics.f'])]
	)


