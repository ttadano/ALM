.. _compile_with_conda_packages:

Building ALM using conda
=========================

ALM is written in C++. To build it, a set of build tools is
needed. Currently using conda gives a relatively simple and uniform
way to perform building the ALM python module and the ALM library for
C++. In this documentation, it is presented a step-by-step procedure
to build them using conda.

.. contents::
   :depth: 2
   :local:

Creating a conda environment
-----------------------------

At first, it is recommended `to prepare a conda environment
<https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_,
which is done by::

   % conda create --name alm

Here the name of the conda environment is chosen ``alm``. The
detailed instruction about the conda environment is found at

- https://conda.io/docs/user-guide/tasks/manage-environments.html

Preparing build tools by conda
-------------------------------

Compilers prepared by conda for linux and macOS are different. See:

- https://conda.io/docs/user-guide/tasks/build-packages/compiler-tools.html

The necessary conda packages to build ALM on linux and macOS are
respectively installed as shown below.

For linux
~~~~~~~~~~

::

   % conda install gcc_linux-64 gxx_linux-64 cmake boost eigen numpy ipython
   % conda install -c conda-forge openblas h5py

For macOS
~~~~~~~~~~

::

   % conda install clang_osx-64 clangxx_osx-64 llvm-openmp cmake boost eigen numpy ipython
   % conda install -c conda-forge openblas h5py

Along with this installation of compilers, conda activation scripts of
environment variables are installed under
``$CONDA_PREFIX/etc/conda/activate.d``. On macOS Mojave, unfortunately
this script doesn't work well at this moment (15/Oct/2018). So the
comment out the line below in the script
``$CONDA_PREFIX/etc/conda/activate.d/activate_clang_osx-64.sh`` is
required::

   # "CONDA_BUILD_SYSROOT,${CONDA_BUILD_SYSROOT:-$(xcrun --show-sdk-path)}"


.. _build_ALMlib:

Building ALM python module and/or ALM library for C++
------------------------------------------------------

Now the directory structure supposed in this document is shown as below::

   $HOME
   |-- ALM
   |   |-- ALM
   |   |   |-- include/
   |   |   |-- lib/
   |   |   |-- python/setup.py
   |   |   |-- src/
   |   |   |-- _build/
   |   |   |-- CMakeLists.txt
   |   |   `-- ...
   |   `-- spglib
   |       |-- include/
   |       |-- lib/
   |       |-- _build/
   |       |-- CMakeLists.txt
   |       `-- ...
   |-- miniconda/envs/alm/include
   |-- miniconda/envs/alm/include/eigen3
   `-- ...

In this directory structure, ``$CONDA_PREFIX`` is equivalent to
``$HOME/miniconda/envs/alm``. The location of ``miniconda`` directory
is chosen at the installation time of miniconda.

ALM and spglib are downloaded from github. ``ALM`` and ``spglib``
directories are created running the following commands::

   % git clone https://github.com/ttadano/ALM.git
   % git clone https://github.com/atztogo/spglib.git

When this is done on ``$HOME/ALM``, the above directory structure is
made. If git command doesn't exist in your system, it is also obtained
from conda by ``conda install git``.

Preparing spglib
~~~~~~~~~~~~~~~~

`spglib <https://github.com/atztogo/spglib>`_ is necessary for to
build both of the ALM python module and/or ALM library for C++. It is
assumed that we are in ``$HOME/ALM/spglib``, where spglib is built as
follows::

   % mkdir _build && cd _build
   % cmake -DCMAKE_INSTALL_PREFIX="" ..
   % make
   % make DESTDIR=.. install

Detailed build configuration is controlled to modify
``CMakeLists.txt`` if necessary.

Building ALM python module by Modifying ``setup.py``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ALM python module is built on the directory
``$HOME/ALM/ALM/python``. So first you have to move to this directory.

To link the spglib static link library (``libsymspg.a``),
``spglib_dir`` in ``setup.py`` has to be modified appropriately. If we
assume it exists in the directory ``$HOME/ALM/spglib/lib`` following
the :ref:`above <build_ALMlib>` directory structure, the modification
is::

   spglib_dir = os.path.join(home, "ALM", "spglib", "lib")

To build ALM, it is necessary to tell compiler where the libraries and
their header files exist. If the conda environment is used, the
library path would be set correctly. However the include paths
probably have to be set for boost and eigen by::

   % export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include:$CONDA_PREFIX/include/eigen3:$HOME/ALM/spglib/include

The build of the ALM python module is performed by a C++ compiler but
not a C compiler. To let the python `setuptools
<https://setuptools.readthedocs.io/en/latest/>`_ choose the C++
compiler installed using conda, the environment variables ``CC`` is
overwritten by ``CXX`` by

::

   % export CC=$CXX

Finally the build and installation in the user directory is done by

::

   % python setup.py build
   % pip install -e .

Building ALM library for C++ by modifying ``CMakeLists.txt``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you need only ALM python module, this section can be skipped.

Let's assume we are in the directory ``$HOME/ALM/ALM`` (see above
:ref:`directory structure <build_ALMlib>`). The ALM
library for C++ is built using cmake. The cmake's configuration file
has to have the filename ``CMakeLists.txt``. So its example of
``CMakeLists.txt.conda`` is renamed to ``CMakeLists.txt``, i.e.,

::

   % mv CMakeLists.txt.conda CMakeLists.txt

Then this ``CMakeLists.txt`` is to be modified appropriately.  At
least, the following lines for spglib library setting would be
modified depending on your location of the spglib library,

::

   include_directories("$ENV{HOME}/ALM/spglib/include")
   set(spglib "-L$ENV{HOME}/ALM/spglib/lib -lsymspg")

These lines are an example made along with the directory
structure shown :ref:`above <build_ALMlib>`. Using this
``CMakeLists.txt``, the ALM library for c++ is built by

::

   % mkdir _build && cd _build
   % cmake ..
   % make -j4
   % make install


The dynamic and static link libraries and the head file are installed
at

- ``$HOME/ALM/ALM/lib/libalmcxx.dylib`` or ``$HOME/ALM/ALM/lib/libalmcxx.so``
- ``$HOME/ALM/ALM/lib/libalmcxx.a``
- ``$HOME/ALM/ALM/include/alm.h``

These libraries are linked to spglib, openblas, and boost
dynamically. Therefore to use the ALM library for C++,
``LD_LIBRARY_PATH`` has to be set properly, e.g.,

::

   export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$HOME/ALM/spglib/lib:$LD_LIBRARY_PATH
