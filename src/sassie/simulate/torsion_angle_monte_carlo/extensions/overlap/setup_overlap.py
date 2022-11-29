'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
# System imports
from distutils.core import *
from distutils      import sysconfig
import os, sys
#os.environ["CC"] = "g++"

# Third-party modules - we depend on numpy for everything
import numpy
#import sassie.util.sasconfig as sasconfig

#sys.path.append('./')
#import sasconfig as sasconfig

try:
    import sassie.util.sasconfig as sasconfig
except:
    sys.path.append('./')
    import _sasconfig as sasconfig

from numpy.distutils.core import Extension, setup

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

core_libraries = os.path.join(os.path.sep, 'usr','local','core_libraries')
core_libraries = os.path.join(os.path.sep, 'share','apps','local','core_libraries')

cpp_library_name = 'overlap'
cuda_library_name = 'cudaOverlap'

#### OPEN NEED INSTALLATION INFORMATION FOR THIS

if sasconfig.__cuda__:
    cuda_dir = sasconfig.__cuda_path__
    if os.path.isfile(os.path.join(core_libraries,'lib','libcudaOverlap.a')):
        cuda_lib = True
    else:
        cuda_lib = False

cuda_driver = True

if os.path.isfile(os.path.join(core_libraries,'lib','liboverlap.a')):
    cpp_lib = True
else:
    print os.path.join(core_libraries,'lib','liboverlap.a')

    cpp_lib = False


#if not cpp_lib and not cuda_lib:
#    print ("Either cpp or cuda lib needs to be pre-built")
#    exit(0)
#if cuda_lib and not cuda_driver:
#    print ("Cuda lib found but no cuda driver detected")
#    exit(0)

include_dir_names = [numpy_include]
library_dir_names = []
library_names = []
macros = []

CUDA = sasconfig.__cuda__

if cpp_lib:
    include_dir_names.append(os.path.join(core_libraries,'include'))
    library_dir_names.append(os.path.join(core_libraries,'lib'))
    library_names.append(cpp_library_name)
    macros.append(('CPP_LIB','1'))
if CUDA and cuda_driver:
    include_dir_names.append(os.path.join(cuda_dir,'include'))
    library_dir_names.append(os.path.join(cuda_dir,'lib64'))
    library_names.append('cuda')
    library_names.append('cudart')
    macros.append(('CUDA_DRIVER','1'))
if CUDA and cuda_lib:
    include_dir_names.append(os.path.join(core_libraries,'include'))
    library_dir_names.append(os.path.join(core_libraries,'lib'))
    library_names.append(cuda_library_name)
    macros.append(('CUDA_LIB','1'))

# extension module
overlap = Extension(name="overlap",sources=['overlap_extension.cpp'],
                    include_dirs = include_dir_names,
                    library_dirs = library_dir_names,
                    libraries = library_names,
                    define_macros = macros,
                   )

# setup
setup(  name        = "Overlap",
        description = "Module checks for atomic overlap",
        author      = "Hailiang Zhang",
        ext_modules = [overlap]
        )

