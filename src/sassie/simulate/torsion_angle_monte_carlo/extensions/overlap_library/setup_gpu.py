"""cudaOverlap installer.

This file is part of the Cadishi package.  See README.rst,
LICENSE.txt, and the documentation for details.
"""

# ez_setup attempts to download setuptools in case it is not available
# import ez_setup
# ez_setup.use_setuptools()
from setuptools import setup, Extension, Command
from Cython.Distutils import build_ext
import numpy
import sys
import os
import subprocess as sub


if "--debug" in sys.argv:
    # flags mainly intended to speed up the builds during development
    print("Using debug compiler flags ...")
    debug_build = True
    sys.argv.remove("--debug")
else:
    debug_build = False


def find_in_path(filenames):
    """Locate executables in PATH."""
    # adapted from
    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52224
    from os.path import exists, join, abspath
    from os import pathsep, environ
    search_path = environ["PATH"]
    paths = search_path.split(pathsep)
    for path in paths:
        for filename in filenames:
            if exists(join(path, filename)):
                return abspath(join(path, filename))


def locate_cuda():
    """Locate the CUDA environment on the system Returns a dict with keys 'home',
    'nvcc', 'include', and 'lib64' and values giving the absolute path to each
    directory. Starts by looking for the CUDAHOME env variable. If not found,
    everything is based on finding 'nvcc' in the PATH.
    """
    # adapted from
    # https://stackoverflow.com/questions/10034325/can-python-distutils-compile-cuda-code
    nvcc = None
    envs = ['CUDA_HOME', 'CUDA_ROOT', 'CUDAHOME', 'CUDAROOT']
    for env in envs:
        if env in os.environ:
            nvcc = os.path.join(os.environ[env], 'bin', 'nvcc')
            break
    else:
        # otherwise, search PATH for NVCC
        nvcc = find_in_path(['nvcc'])
    if nvcc is None:
        raise EnvironmentError('The nvcc executable could not be found.  ' +
            'Add it to $PATH or set one of the environment variables ' +
            ', '.join(envs))
    home = os.path.dirname(os.path.dirname(nvcc))
    #home = os.path.join(os.path.sep,'share','apps','local','cuda')
    #nvcc =  '/share/apps/local/cuda/bin/nvcc'
    cudaconfig = {'home':home,
                  'nvcc':nvcc,
                  'include': os.path.join(home, 'include'),
                  'lib64': os.path.join(home, 'lib64')}
    for k, v in cudaconfig.items():
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be located in %s' % (k, v))
    print("Found CUDA: " + home)
    return cudaconfig
CUDA = locate_cuda()


def get_cuda_ver(nvcc="nvcc"):
    cmd = [nvcc, '--version']
    major = -1
    minor = -1
    patch = -1
    try:
        raw = sub.check_output(cmd).split(b'\n')
        for line in raw:
            if line.startswith(b'Cuda'):
                tokens = line.decode().split(',')
                # we obtain a version string such as "7.5.17"
                verstr = tokens[2].strip().strip('V')
                vertup = verstr.split('.')
                major = int(vertup[0])
                minor = int(vertup[1])
                patch = int(vertup[2])
    except:
        raise
    ver = major, minor, patch
    return ver
CUDAVER = get_cuda_ver(CUDA['nvcc'])
# print CUDAVER


def customize_compiler_for_nvcc(self):
    """inject deep into distutils to customize how the dispatch to gcc/nvcc
    works.  If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going the
    OO route, I have this. Note, it's kindof like a wierd functional subclassing
    going on."""
    # adapted from
    # https://stackoverflow.com/questions/10034325/can-python-distutils-compile-cuda-code

    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')
    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile
    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            postargs = extra_postargs['gcc']
        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so
    # inject our redefined _compile method into the class
    self._compile = _compile

# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


# --- GCC flags (tested with gcc 4.8)
if debug_build:
    gcc_flags = ['-O0', '-g']
else:
    gcc_flags = ['-O3']
    if (find_in_path(['gcc']) is not None):
        # --- let us assume that gcc will be used
        gcc_flags += ['-ffast-math']
        gcc_flags += ['-msse4.2']  # runs well on any reasonably modern system
        # gcc_flags += ['-mavx']  # AVX, found to be slower than SSE 4.2
        #gcc_flags += ['-march=native']  # non portable, not clear what gcc is doing
    else:
        pass
gcc_flags_string = " ".join(gcc_flags)


# --- CUDA flags
if debug_build:
    nvcc_flags = ['-O0', '-g', '-G']
else:
    nvcc_flags = ['-O3']
    nvcc_flags += ['-use_fast_math']
    #nvcc_flags += ['--generate-code', 'arch=compute_20,code=sm_20']
    #nvcc_flags += ['--generate-code', 'arch=compute_20,code=sm_21']
    #nvcc_flags += ['--generate-code', 'arch=compute_30,code=sm_30']
    #nvcc_flags += ['--generate-code', 'arch=compute_32,code=sm_32']
    #nvcc_flags += ['--generate-code', 'arch=compute_35,code=sm_35']
    nvcc_flags += ['--generate-code', 'arch=compute_75,code=sm_75']
    #nvcc_flags += ['--generate-code', 'arch=compute_37,code=sm_37']
    #if (CUDAVER[0] >= 6):
    #    nvcc_flags += ['--generate-code', 'arch=compute_50,code=sm_50']
    #if (CUDAVER[0] >= 7):
    #    nvcc_flags += ['--generate-code', 'arch=compute_52,code=sm_52']
    #    nvcc_flags += ['--generate-code', 'arch=compute_53,code=sm_53']
    #if (CUDAVER[0] >= 8):
    #    pass
        # nvcc_flags += ['--generate-code', 'arch=compute_52,code=sm_52']
        # nvcc_flags += ['--generate-code', 'arch=compute_53,code=sm_53']
nvcc_flags += ['--compiler-options=' + gcc_flags_string + ' -fPIC']


class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    # https://stackoverflow.com/questions/3779915/why-does-python-setup-py-sdist-create-unwanted-project-egg-info-in-project-r
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build')
        os.system('rm -vrf ./cudaOverlap.egg-info')
        os.system('rm -vrf ./cudaOverlap.so')

include_dir_names = [numpy_include,os.path.join(CUDA['home'],'include'),'./extensions/src/']
macros = []
macros.append(('CUDA_LIB','1'))
macros.append(('USE_CUDA','1'))

# extension module
overlap_api = Extension(name="overlap_api",
                     sources=['overlap_extension.cpp',
                              './extensions/src/cudaKernel_overlap.cu',
                              './extensions/src/wrapperCudaKernel.cu',
                     ],
                    include_dirs = include_dir_names,
                    libraries=['cudart', 'stdc++'],
                    library_dirs=[CUDA['lib64']],
                    define_macros = macros,
                    runtime_library_dirs=[CUDA['lib64']],
                    extra_compile_args={'gcc':gcc_flags, 'nvcc':nvcc_flags}
                   )
# setup
setup(  name        = "overlap_api",
        description = "Module for overlap_api",
        author      = "Hailiang Zhang",
        packages=[os.path.join(os.getcwd(),'extensions','src')],
        ext_modules = [overlap_api],
        install_requires=['numpy'],
        cmdclass={'build_ext': custom_build_ext, 'clean': CleanCommand},
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
