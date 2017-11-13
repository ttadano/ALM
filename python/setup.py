import os
from setuptools import setup, Extension
# from distutils.core import setup, Extension
import numpy
include_dirs_numpy = [numpy.get_include()]

compile_with_sources = True

extra_link_args = ['-lstdc++', '-lgomp', '-llapack']
library_dirs = []
if compile_with_sources:
    sources = ['alm.cpp',
               'alm_core.cpp',
               'alm_cui.cpp',
               'constraint.cpp',
               'error.cpp',
               'fcs.cpp',
               'files.cpp',
               'fitting.cpp',
               'input_parser.cpp',
               'input_setter.cpp',
               'interaction.cpp',
               'main.cpp',
               'patterndisp.cpp',
               'symmetry.cpp',
               'system.cpp',
               'timer.cpp',
               'writer.cpp']
    if os.path.exists('src'):
        source_dir = "src"
    else:
        source_dir = "../src"
    
    include_dirs = [source_dir,]
    include_dirs += include_dirs_numpy
    for i, s in enumerate(sources):
        sources[i] = "%s/%s" % (source_dir, s) 
    sources += ['_alm.c', 'alm_wrapper.cpp']
else: # compile with library
    sources = ['_alm.c', 'alm_wrapper.cpp']
    # static link library
    extra_link_args += ['../lib/libalmcxx.a']
    # dynamic link library
    # extra_link_args += ['-lalmcxx']
    # library_dirs = ['../lib']

extension = Extension('alm._alm',
                      include_dirs=include_dirs_numpy,
                      library_dirs=library_dirs,
                      extra_compile_args = ['-fopenmp', '-std=c++11'],
                      extra_link_args=extra_link_args,
                      sources=sources)

setup(name='alm',
      version='0.9.8',
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
