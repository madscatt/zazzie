'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
# System imports
from distutils.core import *
from distutils      import sysconfig
import os
os.environ["CC"] = "/share/apps/local/bin/gcc"
os.environ["CXX"] = "/share/apps/local/bin/g++"

# Third-party modules - we depend on numpy for everything
import numpy

from numpy.distutils.core import Extension, setup

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

cpp_library_name = 'sasmol'
core_library_path=os.path.join(os.path.sep,'share','apps','local','core_libraries')

if os.path.isfile(os.path.join(core_library_path,'lib','libsasmol.a')):
    cpp_lib = True
else:
    cpp_lib = False

include_dir_names = [numpy_include]
library_dir_names = ['/share/apps/local/boost_1_55_0/libs/']
library_names = ['boost_regex-mt']
macros = []

if cpp_lib:
    include_dir_names.append(os.path.join(core_library_path,'include'))
    library_dir_names.append(os.path.join(core_library_path,'lib'))
    library_names.append(cpp_library_name)
    macros.append(('CPP_LIB','1'))

# extension module
get_bead_masks = Extension(name="get_bead_masks",sources=['get_bead_masks_extension.cpp'],
                    include_dirs = include_dir_names,
                    library_dirs = library_dir_names,
                    libraries = library_names,
                    extra_compile_args=['-std=c++11'],
                    define_macros = macros,
                   )

# setup
setup(  name        = "get_bead_masks",
        description = "Module of get_bead_masks",
        author      = "Hailiang Zhang",
        ext_modules = [get_bead_masks]
        )

