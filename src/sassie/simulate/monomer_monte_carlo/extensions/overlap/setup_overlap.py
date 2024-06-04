'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
from setuptools import setup, Extension
import numpy

overlap_module = Extension('overlap',
                           sources=['overlap.c'],
                           include_dirs=[numpy.get_include()])

setup(
    name='Overlap',
    version='1.0',
    description='Python package with C extension for overlap calculation',
    ext_modules=[overlap_module],
)
