import os
from setuptools import setup, Extension
import numpy
try:
    from pathlib import Path
    home = str(Path.home())
except ImportError:
    home = os.path.expanduser("~")

# This is the switch for ALM developement. Always True for general cases.
compile_with_sources = True

# Configurations to pass to extention.
# The following directory structure and use of conda are supposed.
#
# $HOME
# |-- ALM
# |   `-- ALM
# |       |-- include/
# |       |-- lib/
# |       |-- python/setup.py
# |       |-- src/
# |       |-- _build/
# |       |-- CMakeLists.txt
# |       `-- ...
#
# |-- $CONDA_PREFIX/include
# |-- $CONDA_PREFIX/include/eigen3
# `-- ...

if 'CONDA_PREFIX' in os.environ:
    conda_prefix = os.environ['CONDA_PREFIX']
else:
    conda_prefix = os.path.join(home, "miniconda", "envs", "alm")

library_dirs = []

# For linux
extra_link_args = ['-lgomp', '-llapack']
# For macOS
#extra_link_args = ['-lomp']

spglib_dir = os.path.join(conda_prefix, "lib")

include_dirs = []
include_dirs.append(numpy.get_include())
include_dirs.append(os.path.join(conda_prefix, "include"))
include_dirs.append(os.path.join(conda_prefix, "include", "eigen3"))

if compile_with_sources:
    cpp_files = ['alm.cpp',
                 'alm_cui.cpp',
                 'cluster.cpp',
                 'constraint.cpp',
                 'fcs.cpp',
                 'files.cpp',
                 'input_parser.cpp',
                 'input_setter.cpp',
                 'main.cpp',
                 'optimize.cpp',
                 'patterndisp.cpp',
                 'rref.cpp',
                 'symmetry.cpp',
                 'system.cpp',
                 'timer.cpp',
                 'writer.cpp']
    if os.path.exists('src'):
        source_dir = "src"
    else:
        source_dir = os.path.join("..", "src")

    include_dirs += [source_dir, ]
    sources = [os.path.join(source_dir, s) for s in cpp_files]
    sources += ['_alm.c', 'alm_wrapper.cpp']
    extra_link_args.append(os.path.join(spglib_dir, "libsymspg.a"))
else:  # compile with library
    sources = ['_alm.c', 'alm_wrapper.cpp']
    # 'libalmcxx.a' static link library has to come before depending
    # dynamic link libraries '-lgomp', '-llapack'.
    extra_link_args.insert(0, os.path.join("..", "lib", "libalmcxx.a"))
    extra_link_args.append(os.path.join(spglib_dir, "libsymspg.a"))

extension = Extension('alm._alm',
                      include_dirs=include_dirs,
                      library_dirs=library_dirs,
                      extra_compile_args = ['-fopenmp', '-std=c++11',
                                            '-DWITH_SPARSE_SOLVER'],
                      extra_link_args=extra_link_args,
                      sources=sources)

setup(name='alm',
      version='1.0.2',
      description='Force constants generator',
      setup_requires=['numpy', 'setuptools>=18.0'],
      author='Terumasa Tadano',
      author_email='terumasa.tadano@gmail.com',
      url='https://github.com/ttadano/ALM',
      packages=['alm'],
      install_requires=['numpy'],
      provides=['alm'],
      platforms=['all'],
      ext_modules=[extension])
