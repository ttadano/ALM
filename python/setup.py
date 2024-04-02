"""Setuptools script for building alm package."""

import os
import warnings

import numpy
from setuptools import Extension, setup

try:
    from pathlib import Path

    home = str(Path.home())
except ImportError:
    home = os.path.expanduser("~")

# This is the switch for ALM development. Always True for general cases.
compile_with_sources = True

# Configurations to pass to extension.
# The following directory structure and use of conda are supposed.
#
# $HOME
# |-- ALM
# |   `-- ALM
# |       |-- bin/
# |       |-- include/
# |       |-- lib/
# |       |-- python/setup.py
# |       |-- src/
# |       |-- _build/
# |       |-- CMakeLists.txt
# |       `-- ...
# |
# |-- $CONDA_PREFIX/include
# |-- $CONDA_PREFIX/include/eigen3
# |-- $CONDA_PREFIX/lib
# `-- ...

if "CONDA_PREFIX" in os.environ:
    conda_prefix = os.environ["CONDA_PREFIX"]
else:
    conda_prefix = os.path.join(home, "miniconda", "envs", "alm")

library_dirs = []
extra_link_args = []

# OpenMP library
if "CC" in os.environ:
    if "clang" in os.environ["CC"]:
        if "clang++" not in os.environ["CC"]:
            warnings.warn("clang++ is used instead of clang.")
            os.environ["CC"] = os.environ["CC"].replace("clang", "clang++")
        extra_link_args.append("-lomp")
    elif "gcc" in os.environ["CC"]:
        if "g++" not in os.environ["CC"]:
            warnings.warn("g++ is used instead of gcc.")
            os.environ["CC"] = os.environ["CC"].replace("gcc", "g++")
        extra_link_args.append("-lgomp")
    elif "gnu-cc" in os.environ["CC"]:
        if "gnu-c++" not in os.environ["CC"]:
            warnings.warn("gnu-g++ is used instead of gnu-cc.")
            os.environ["CC"] = os.environ["CC"].replace("gnu-cc", "gnu-c++")
        extra_link_args.append("-lgomp")

if not extra_link_args:  # Default libgomp
    extra_link_args.append("-lgomp")

# Lapack library
extra_link_args.append("-llapack")

# spglib
extra_link_args.append("-lsymspg")

include_dirs = []
include_dirs.append(numpy.get_include())
include_dirs.append(os.path.join(conda_prefix, "include"))
include_dirs.append(os.path.join(conda_prefix, "include", "eigen3"))


if compile_with_sources:
    cpp_files = [
        "alm.cpp",
        "cluster.cpp",
        "constraint.cpp",
        "fcs.cpp",
        "files.cpp",
        "optimize.cpp",
        "patterndisp.cpp",
        "rref.cpp",
        "symmetry.cpp",
        "system.cpp",
        "timer.cpp",
        "writer.cpp",
    ]
    if os.path.exists("src"):
        source_dir = "src"
    else:
        source_dir = os.path.join("..", "src")

    include_dirs += [
        source_dir,
    ]
    sources = [os.path.join(source_dir, s) for s in cpp_files]
    sources += ["_alm.c", "alm_wrapper.cpp"]
else:  # compile with library
    sources = ["_alm.c", "alm_wrapper.cpp"]
    # 'libalmcxx.a' static link library has to come before depending
    # dynamic link libraries '-lgomp', '-llapack'.
    extra_link_args.insert(0, os.path.join("..", "lib", "libalmcxx.a"))

extension = Extension(
    "alm._alm",
    include_dirs=include_dirs,
    library_dirs=library_dirs,
    extra_compile_args=["-fopenmp", "-std=c++11"],
    extra_link_args=extra_link_args,
    sources=sources,
)

setup(
    platforms=["all"],
    ext_modules=[extension],
)
