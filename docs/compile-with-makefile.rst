
Building ALM using Makefile
===========================

ALM can also be built with the Makefile in the ``src`` subdirectory.
This approach only generates the command line version of ALM (binary ``alm``).
Therefore, if you want to use ALM from python as well, please 
:ref:`build ALM using conda <compile_with_conda_packages>`.

Requirement
-----------

* C++ compiler (C++11 standard or newer)
* LAPACK 
* `Boost C++ library <http://www.boost.org>`_
* `Eigen3 library <http://eigen.tuxfamily.org/>`_
* `spglib <https://atztogo.github.io/spglib/>`_

How to install
--------------

.. highlight:: bash

1. Install the LAPACK, Boost C++, and Eigen3, and spglib.

   To install the Boost C++ library, please download a source file from the `website <http://www.boost.org>`_ and
   unpack the file. Then, copy the 'boost' subdirectory to the include folder in the home directory (or anywhere you like).
   This can be done as follows::
    
    $ cd
    $ mkdir etc; cd etc
    (Download a source file and mv it to ~/etc)
    $ tar xvf boost_x_yy_z.tar.bz2
    $ cd ../
    $ mkdir include; cd include
    $ ln -s ../etc/boost_x_yy_z/boost .

  In this example, we made a symbolic link to the 'boost' subdirectory in ``$HOME/include``.
  Instead of installing from source, you can install the Boost library with `Homebrew <http://brew.sh>`_ on Mac.

  In the same way, please install the Eigen3 include files as follows::

    $ cd
    $ mkdir etc; cd etc
    (Download a source file and mv it to ~/etc)
    $ tar xvf eigen-eigen-*.tar.bz2 (* is an array of letters and digits)
    $ cd ../
    $ cd include
    $ ln -s ../etc/eigen-eigen-*/Eigen .  

2. Clone ALM from the git repository and edit Makefile::

    $ git clone https://github.com/ttadano/ALM
    $ cd ALM/src/
    (Edit Makefile.linux or Makefile.osx)
    $ make -f Makefile.linux -j (or make -j Makefile.osx -j)

   In the ``src`` directory, we provide sample Makefiles for linux (Intel compiler) and MacOS (GCC installed via homebrew) as shown below.

.. literalinclude:: ../src/Makefile.linux
    :caption: **Makefile.linux**
    :language: makefile
    :linenos:
    :lines: 7-17

.. literalinclude:: ../src/Makefile.osx
    :caption: **Makefile.osx**
    :language: makefile
    :linenos:
    :lines: 7-18

3. Modify ``LD_LIBRARY_PATH`` as follows::

    bash, zsh
    $ export LD_LIBRARY_PATH=$(HOME)/src/spglib/lib:$LD_LIBRARY_PATH

    csh, tcsh
    $ setenv LD_LIBRARY_PATH $(HOME)/src/spglib/lib:$LD_LIBRARY_PATH
