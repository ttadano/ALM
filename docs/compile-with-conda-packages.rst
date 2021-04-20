.. _compile_with_conda_packages:

Building ALM using conda
=========================

ALM is written in C++. To build it, a set of build tools is
needed. Currently using conda gives a relatively simple and uniform
way to perform building the ALM python module, ALM library for
C++, and ALMM command executable. In this documentation, it is
presented a step-by-step procedure to build them using conda.

.. contents::
   :depth: 2
   :local:

Preparing build tools by conda
-------------------------------

At first, it is recommended `to prepare a conda environment
<https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_ by::

   % conda create --name alm -c conda-forge python=3.8

Here the name of the conda environment is chosen ``alm``. The detailed
instruction about the conda environment is found `here
<https://conda.io/docs/user-guide/tasks/manage-environments.html>`_.
To build ALM on linux or macOS, the following conda packages are
installed by

::

   % conda install -c conda-forge numpy scipy h5py compilers "libblas=*=*mkl" spglib boost eigen cmake ipython mkl-include


.. _build_ALMlib:

Building ALM
-------------

Now the directory structure supposed in this document is shown as below::

   $HOME
   |-- alm
   |   `-- ALM
   |       |-- bin/
   |       |-- include/
   |       |-- lib/
   |       |-- python/setup.py
   |       |-- src/
   |       |-- _build/
   |       |-- CMakeLists.txt
   |       |-- CMakeLists.txt.conda
   |       `-- ...
   |
   |-- $CONDA_PREFIX/include
   |-- $CONDA_PREFIX/include/eigen3
   |-- $CONDA_PREFIX/lib
   `-- ...

``alm`` directory is created by us and we move to this directory. Now
we are in ``$HOME/alm``. In this direcotry, ALM is downloaded from
github. ``ALM`` directorie is created running the following commands::

   % git clone https://github.com/ttadano/ALM.git

Make sure that you are using develop branch by

::

   % cd ALM
   % git branch
   * develop

When this is done on ``$HOME/ALM``, the above directory structure is
made. If git command doesn't exist in your system, it is also obtained
from conda by ``conda install git``.

Building ALM python module
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ALM python module is built on the directory
``$HOME/alm/ALM/python``. So first you have to move to this directory.
The build and installation in the user directory is done by

::

   % python setup.py build
   % pip install -e .

..
   For macOS, we use clang instead of gcc in this documentation. In this
   case, ALM python module must be compiled by clang++
   command but not clang command. To let python `setuptools
   <https://setuptools.readthedocs.io/en/latest/>`_ choose the C++
   compiler installed using conda, the environment variables ``CC`` is
   overwritten by ``CXX`` by

   ::

      % export CC=$CXX

   and libomp is used as the OpenMP library, which is set in ``setup.py``

      extra_link_args = ['-lomp']


Building ALM executable and C++ library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you need only ALM python module, this section can be skipped.

Let's assume we are in the directory ``$HOME/alm/ALM`` (see above
:ref:`directory structure <build_ALMlib>`). The ALM
library for C++ is built using cmake. The cmake's configuration file
has to have the filename ``CMakeLists.txt``. So its example of
``CMakeLists.txt.conda`` is renamed to ``CMakeLists.txt``, i.e.,

::

   % cp CMakeLists.txt.conda CMakeLists.txt

Then this ``CMakeLists.txt`` may be modified appropriately when the
following compilation fails.
Using this ``CMakeLists.txt``, the ALM library for c++ is built for Linux by

::

   % mkdir _build && cd _build
   % cmake ..
   % make -j4
   % make install

and for macOS

::

   % mkdir _build && cd _build
   % cmake -DCMAKE_C_COMPILER='clang' -DCMAKE_CXX_COMPILER='clang++' ..
   % make -j4
   % make install

To see detailed make log, ``-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON`` option
for ``cmake`` should be added.

The dynamic and static link libraries and the head file are installed
at

- ``$HOME/alm/ALM/lib/libalmcxx.dylib`` or ``$HOME/alm/ALM/lib/libalmcxx.so``
- ``$HOME/alm/ALM/lib/libalmcxx.a``
- ``$HOME/alm/ALM/include/alm.h``

The executable is found at

- ``$HOME/alm/ALM/bin/alm``

To use the dynamic link library, it may be necessary to set
``$LD_LIBRARY_PATH`` to

::

   % export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$HOME/alm/ALM/lib:$LD_LIBRARY_PATH

and to use the executable

::

   % export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
