from __future__ import print_function
# System imports
from distutils.core import *
import os 
import platform

os_type = platform.system()

#os.environ["CC"] = "g++"

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
overlap_api = Extension(name="overlap_api",sources=['overlap_extension.cpp',
                                                    'extensions/src/overlap.cpp',
                                                   ], 
                    include_dirs = include_dir_names,
                    define_macros = macros,
                   )

# setup
setup(  name        = "overlap_api",
        description = "Module for overlap_api",
        author      = "Hailiang Zhang",
        ext_modules = [overlap_api]
        )



### post compilation file move

try:
    lib_file = os.path.join('build', 'lib*', 'overlap_api.*')
    os.system('mv ' + lib_file + ' .')
except:
    print('overlap_api.* not found')
    print('overlap_api.* not found')
    print('overlap_api.* not found')
    print('\nINSTALLATION FAILED\n')
    print('INSTALLATION FAILED\n')
    print('INSTALLATION FAILED\n')
