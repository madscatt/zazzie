from setuptools import setup, Extension
import numpy

module = Extension('pairs', sources=['pairs.c'], include_dirs=[numpy.get_include()])

setup(
    name='PairsModule',
    version='1.0',
    description='Python package with C extension',
    ext_modules=[module],
)
