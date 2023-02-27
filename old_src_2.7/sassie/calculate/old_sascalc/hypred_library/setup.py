from __future__ import print_function
# System imports
from distutils.core import *
import os, platform

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

include_dir_names = [numpy_include,'extensions/include/']

# extension module
hypred_api = Extension(name="hypred_api",sources=['hypred_library_extension.cpp',
                                                    'extensions/src/util.cc',
                                                    'extensions/src/coor.cc',
                                                    'extensions/src/Atom.cc',
                                                    'extensions/src/box.cc',
                                                    'extensions/src/pRDF.cc',
                                                   ], 
		            extra_compile_args=['-std=c++0x','-g','-O0'],
                    include_dirs = include_dir_names,
                   )

# setup
setup(  name        = "hypred_api",
        description = "Module for hypred_api",
        author      = "Hailiang Zhang",
        ext_modules = [hypred_api]
        )



### post compilation file move

try:
    lib_file = os.path.join('build', 'lib*', 'hypred_api.*')
    os.system('mv ' + lib_file + ' .')
except:
    print('hypred_api.* not found')
    print('hypred_api.* not found')
    print('hypred_api.* not found')
    print('\nINSTALLATION FAILED\n')
    print('INSTALLATION FAILED\n')
    print('INSTALLATION FAILED\n')
