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
# |   |-- ALM
# |   |   |-- include/
# |   |   |-- lib/
# |   |   |-- python/setup.py
# |   |   |-- src/
# |   |   |-- _build/
# |   |   |-- CMakeLists.txt
# |   |   `-- ...
# |   `-- spglib
# |       |-- include/
# |       |-- lib/
# |       |-- _build/
# |       |-- CMakeLists.txt
# |       `-- ...
# |-- miniconda/envs/alm/include
# |-- miniconda/envs/alm/include/eigen3
# `-- ...

library_dirs = []
extra_link_args = []
include_dirs = []

spglib_dir = os.path.join(home, "ALM", "spglib", "lib")
#spglib_dir = os.path.join(home, "src", "spglib", "lib")
include_dirs_numpy = [numpy.get_include()]
include_dirs += include_dirs_numpy

# The following setting can work in place of "export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include:$CONDA_PREFIX/include/eigen3:$HOME/ALM/spglib/include".
# include_dir_spglib = os.path.join(home, "ALM", "spglib", "include")
# include_dirs.append(include_dir_spglib)
# include_dir_boost = os.path.join(home, "miniconda", "envs", "alm", "include")
# include_dirs.append(include_dir_boost)
# include_dir_eigen = os.path.join(include_dir_boost, "eigen3")
# include_dirs.append(include_dir_eigen)

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
    # static link library
    extra_link_args.append(os.path.join("..", "lib", "libalmcxx.a"))
    extra_link_args.append(os.path.join(spglib_dir, "libsymspg.a"))
    # dynamic link library
    # extra_link_args += ['-lalmcxx']
    # library_dirs.append(os.path.join("..", "lib"))

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
