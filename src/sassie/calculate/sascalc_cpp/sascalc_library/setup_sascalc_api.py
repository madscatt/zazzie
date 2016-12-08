# System imports
from distutils.core import *
import os, platform

os_type = platform.system()

os.environ["CC"] = "h5cc"

###alias h5c++='/usr/bin/h5c++
### https://github.com/ContinuumIO/anaconda-issues/issues/953

### USER EDIT
cpp_buildingBlock_dir=os.path.join('extensions')

### END USER EDIT

cpp_library_name = 'sascalc'

# Third-party modules - we depend on numpy for everything
import numpy

from numpy.distutils.core import Extension, setup

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

if os.path.isfile(os.path.join(cpp_buildingBlock_dir,'lib','libsascalc.a')):
    cpp_lib = True
else:
    cpp_lib = False

cpp_lib = True

if not cpp_lib:
    print ("Cpp lib needs to be pre-built")
    exit(0)

include_dir_names = [numpy_include]
library_dir_names = []
library_names = []
rpath = []
macros = []


if cpp_lib:
    include_dir_names.append(os.path.join(cpp_buildingBlock_dir,'include'))
    include_dir_names.append('/usr/local/include')
    library_dir_names.append(os.path.join(cpp_buildingBlock_dir,'lib'))
    library_dir_names.append(os.path.join('/usr/local/lib'))
    library_names.append(cpp_library_name)
    library_names.append("hdf5_cpp")
    library_names.append("hdf5")
    macros.append(('CPP_LIB','1'))

# extension module
sascalc_api = Extension(name="sascalc_api",sources=['sascalc_api_extension.cpp'],
                    include_dirs = include_dir_names,
                    library_dirs = library_dir_names,
                    libraries = library_names,
                    define_macros = macros,
                   )

# setup
setup(  name        = "sascalc_api",
        description = "Module for sascalc_api",
        author      = "Hailiang Zhang",
        ext_modules = [sascalc_api]
        )

