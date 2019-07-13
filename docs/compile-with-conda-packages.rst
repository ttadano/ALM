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

   % conda install -c conda-forge gcc_linux-64 gxx_linux-64 h5py scipy numpy boost eigen cmake spglib ipython

For macOS
~~~~~~~~~~

::

   % conda install
   % conda install -c conda-forge clang_osx-64 clangxx_osx-64 llvm-openmp cmake boost eigen numpy scipy h5py spglib ipython

..
   Along with this installation of compilers, conda activation and
   deactivation scripts of environment variables are installed under
   ``$CONDA_PREFIX/etc/conda/activate.d`` and
   ``$CONDA_PREFIX/etc/conda/deactivate.d``, respectively.

.. _build_ALMlib:

Building ALM python module and/or ALM library for C++
------------------------------------------------------

Now the directory structure supposed in this document is shown as below::

   $HOME
   |-- ALM
   |   `-- ALM
   |       |-- include/
   |       |-- lib/
   |       |-- python/setup.py
   |       |-- src/
   |       |-- _build/
   |       |-- CMakeLists.txt
   |       `-- ...
   |
   |-- $CONDA_PREFIX/include
   |-- $CONDA_PREFIX/include/eigen3
   `-- ...

ALM is downloaded from github. ``ALM`` directorie is created running
the following commands::

   % git clone https://github.com/ttadano/ALM.git

When this is done on ``$HOME/ALM``, the above directory structure is
made. If git command doesn't exist in your system, it is also obtained
from conda by ``conda install git``.

Building ALM python module by Modifying ``setup.py``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ALM python module is built on the directory
``$HOME/ALM/ALM/python``. So first you have to move to this directory.

For linux, the following line in ``setup.py`` has to be set::

   extra_link_args = ['-lgomp', '-llapack']

For macOS, the build of the ALM python module is performed by a C++
compiler but not a C compiler. To let the python `setuptools
<https://setuptools.readthedocs.io/en/latest/>`_ choose the C++
compiler installed using conda, the environment variables ``CC`` is
overwritten by ``CXX`` by

::

   % export CC=$CXX

In addition for macOS, the following line in ``setup.py`` has to be set::

   extra_link_args = ['-lomp']

Finally the build and installation in the user directory is done by

::

   % python setup.py build
   % pip install -e .

Building ALM library for C++ by modifying ``CMakeLists.txt``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you need only ALM python module, this section can be skipped. Note
that the ``CMakeLists.txt.conda`` doesn't work well for macOS.

Let's assume we are in the directory ``$HOME/ALM/ALM`` (see above
:ref:`directory structure <build_ALMlib>`). The ALM
library for C++ is built using cmake. The cmake's configuration file
has to have the filename ``CMakeLists.txt``. So its example of
``CMakeLists.txt.conda`` is renamed to ``CMakeLists.txt``, i.e.,

::

   % mv CMakeLists.txt.conda CMakeLists.txt

Then this ``CMakeLists.txt`` is to be modified appropriately.
Using this ``CMakeLists.txt``, the ALM library for c++ is built by

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
