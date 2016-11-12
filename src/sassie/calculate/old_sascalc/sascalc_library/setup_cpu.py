from __future__ import print_function
# System imports
from distutils.core import *
import os, platform

os_type = platform.system()

#os.environ["CC"] = "g++"

if os_type == "Darwin":
    cuda_dir = os.path.join(os.path.sep,'usr','local','cuda')
elif os_type == "Linux":
    cuda_dir = os.path.join(os.path.sep,'share','apps','local','cuda')

# Third-party modules - we depend on numpy for everything
import numpy

from numpy.distutils.core import Extension, setup

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

include_dir_names = [numpy_include,'extensions/src/']
macros = []
macros.append(('CPP_LIB','1'))
macros.append(('USE_CPU','1'))

# extension module
sascalc_api = Extension(name="sascalc_api",sources=['sascalc_library_extension.cpp',
                                                    'extensions/src/SasCalc.cpp',
                                                    'extensions/src/Debye.cpp',
                                                    'extensions/src/GV.cpp',
                                                   ], 
                    include_dirs = include_dir_names,
                    define_macros = macros,
                   )

# setup
setup(  name        = "sascalc_api",
        description = "Module for sascalc_api",
        author      = "Hailiang Zhang",
        ext_modules = [sascalc_api]
        )



### post compilation file move

try:
    lib_file = os.path.join('build', 'lib*', 'sascalc_api.*')
    os.system('mv ' + lib_file + ' .')
except:
    print('sascalc_api.* not found')
    print('sascalc_api.* not found')
    print('sascalc_api.* not found')
    print('\nINSTALLATION FAILED\n')
    print('INSTALLATION FAILED\n')
    print('INSTALLATION FAILED\n')
