### for Python 3 compatibility
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

'''
    SASSIE  Copyright (C) 2011-2016 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os,sys
from distutils.core import setup
from distutils      import sysconfig
from numpy.distutils.core import Extension, setup

#

#       SETUP
#
#       12/01/2009      --      initial coding              :       jc
#       11/22/2014      --      adapted for 2.0             :       jc
#       07/21/2016      --      adapted for 2.0 git branch  :       jc
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

### begin user edit ###
### begin user edit ###
### begin user edit ###

all_packages = ['sassie', 'sassie.util','sassie.build', 
    'sassie.interface', 'sassie.interface.align', 
    'sassie.tools', 'sassie.tools.align',
    'sassie.interface', 'sassie.interface.chi_square_filter', 
    'sassie.analyze', 'sassie.analyze.chi_square_filter',
    'sassie.calculate', 'sassie.calculate.sascalc',
    'sassie.interface.sascalc',
    'sassie.calculate.sascalc.sascalc_library',
    'sassie.simulate', 
    'sassie.simulate.torsion_angle_monte_carlo',
    'sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities',
    'sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.tamc_utilities',
    'sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.protein_backbone_torsion',
    'sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.single_stranded_nucleic_backbone_torsion',
    'sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.isopeptide_bond_torsion',
    'sassie.simulate.torsion_angle_monte_carlo.monte_carlo_utilities.double_stranded_nucleic',
    'sassie.interface.torsion_angle_monte_carlo',
    'sassie.simulate.energy', 
    'sassie.simulate.constraints',
    'sassie.simulate.prody',
    'sassie.interface.prody',
    'sassie.build',
    'sassie.build.pdbrx',
    'sassie.build.pdbrx.pdbrx',
    'sassie.build.pdbscan',
    'sassie.build.pdbscan.pdbscan'
    ]

### end user edit ###
### end user edit ###
### end user edit ###

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# for compatibility with from __future__ import unicode_literals

python_2_flag = True

if int(sys.version[0]) > 2:
    python_2_flag = False

if python_2_flag:
    all_packages = [x.encode('UTF8') for x in all_packages] 
    # the next line does NOT work
    #all_data_files = [x.encode('UTF8') for x in all_data_files] 


setup(name='sassie',
	version='2.00_rev_0',
	author='Joseph E. Curtis',
	author_email='joseph.curtis@nist.gov',
	license='GPL 3',
	url='www.smallangles.net/sassie',
	platforms='Linux, Mac OS X',
	description=("A suite of programs to generate atomistic models of biological systems, calculate scattering observables, and compare results to experimental data"),
	long_description=read('README.md'),
	classifiers=["Development Status :: 2.0 Pre-Alpha",
		"License :: OSI Approved :: GNU Public License 3",
		"Intended Audience :: Science/Research",
		"Natural Language :: English",
		"Operating System :: Linux :: MacOS :: MacOS X",
		"Programming Language :: Python :: 2 :: Only",
		"Topic :: Scientific/Engineering :: Chemistry :: Physics"],

    package_dir={'sassie':os.path.join('src','sassie')},
	
    packages = all_packages,

    package_data = {'':['sascalc_api.so']},

    ext_modules=[
        Extension('sassie.simulate.torsion_angle_monte_carlo.ooverlap',['src/sassie/simulate/torsion_angle_monte_carlo/extensions/ooverlap/ooverlap.c'],include_dirs=[numpy_include]),
    Extension('sassie.simulate.torsion_angle_monte_carlo.dna_overlap',['src/sassie/simulate/torsion_angle_monte_carlo/extensions/dna_overlap/dna_overlap.f'])
        ]

    )

#    this data type must be parsed correctly if using string literals !!! ouch
#    and it needs to be in setup() above

#    data_files = all_data_files
