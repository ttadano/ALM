Installation
============

Requirement
-----------

* C++ compiler (Intel compiler is recommended.)
* LAPACK library
* `Boost C++ library <http://www.boost.org>`_
* `Eigen3 library <http://eigen.tuxfamily.org/>`_
* `spglib library <https://atztogo.github.io/spglib/>`_

How to install
--------------

.. highlight:: bash

Install using CMake
~~~~~~~~~~~~~~~~~~~


Install using Makefile
~~~~~~~~~~~~~~~~~~~~~~

To build the binary by Makefile, please edit the Makefile in the src/ directory and issue make as::

    $ cd src/
    $ cp Makefile.linux Makefile
    (Edit Makefile)
    $ make -j

.. highlight:: makefile

Here's a typical setting for Linux with Intel compiler::

    CXX = icpc 
    CXXFLAGS = -O2 -xHOST -qopenmp -std=c++11 -DWITH_SPARSE_SOLVER
    INCLUDE = -I../include -I$(HOME)/include -I$(HOME)/src/spglib/include

    CXXL = ${CXX}
    LDFLAGS = -mkl -L$(HOME)/src/spglib/lib/ -lsymspg

    LAPACK = 
    LIBS = ${LAPACK}

Here's a setting for GCC on MacOS (installed via Homebrew)::

    CXX = g++-8
    CXXFLAGS = -O2 -fopenmp -DWITH_SPARSE_SOLVER
    INCLUDE = -I../include -I$(HOME)/include -I$(HOME)/src/spglib/include

    CXXL = ${CXX}
    LDFLAGS = -lgomp $(HOME)/src/spglib/lib/libsymspg.a

    LAPACK = -llapack -lblas
    LIBS = ${LAPACK}

To enable OpenMP parallelization, please add the ``-qopenmp`` (Intel) or ``-fopenmp`` (gcc) option in ``CXXFLAGS``.
In addition, the directory containing the boost, Eigen, and spglib header files must be given in ``INCLUDE``. 

