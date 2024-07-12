'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''

from setuptools import setup, Extension, find_packages
import numpy

# Define the extension module
dna_overlap = Extension(
    'dna_overlap',  # Name of the module
    sources=['dna_overlap.c'],  # Source files
    include_dirs=[numpy.get_include()],  # Include NumPy headers
    extra_compile_args=['-std=c99'],  # Optional: specify C standard if needed
)

# Setup function to specify details of the module
setup(
    name='dna_overlap',  
    version='0.2',  
    description='Python interface for the dna_overlap C extension',
    ext_modules=[dna_overlap],  # Extensions to build
)

