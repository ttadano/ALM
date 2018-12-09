How to make a DFSET file
------------------------

.. _label_format_DFILE:

Format of ``DFSET``
~~~~~~~~~~~~~~~~~~~

The displacement-force data sets obtained by first-principles (or classical force-field) calculations
have to be saved to a file, say *DFSET*. Then, the force constants are estimated by setting ``DFSET = `` *DFSET* and with ``MODE = optimize``.

The *DFSET* file must contain the atomic displacements and corresponding forces in Cartesian coordinate for at least ``NDATA`` structures (displacement patterns)
in the following format: 

.. math::
    :nowrap:
     
    
    \text{# Structure number 1 (this is just a comment line)}
    \begin{eqnarray*}
     \text{# Structure number 1 (this is just a comment line)} & & & & & \\
     u_{x}(1) & u_{y}(1) & u_{z}(1) & f_{x}(1) & f_{y}(1) & f_{z}(1) \\
     u_{x}(2) & u_{y}(2) & u_{z}(2) & f_{x}(2) & f_{y}(2) & f_{z}(2) \\
              & \vdots   &          &          & \vdots   &          \\
     u_{x}(\mathrm{NAT}) & u_{y}(\mathrm{NAT}) & u_{z}(\mathrm{NAT}) & f_{x}(\mathrm{NAT}) & f_{y}(\mathrm{NAT}) & f_{z}(\mathrm{NAT}) \\
     \text{# Structure number 2} & & & & &  \\
     u_{x}(1) & u_{y}(1) & u_{z}(1) & f_{x}(1) & f_{y}(1) & f_{z}(1) \\
    \end{eqnarray*}

Here, ``NAT`` is the number of atoms in the supercell. 
The unit of displacements and forces must be **Bohr** and **Ryd/Bohr**, respectively.


Generation of ``DFSET`` by extract.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The script ``extract.py`` in the tools directory of ALM is useful to generate ``DFSET`` files from output files of some popular DFT codes.
Let us assume that we have calculated atomic forces for 10 different structures by VASP and saved the results as vasprun_01.xml ... vasprun_10.xml.
Then, a ``DFSET`` file can be generated as::

    **VASP**
    ::

    $ python extract.py --VASP=SPOSCAR --offset=vasprun0.xml vasprun??.xml > DFSET

Here, ``SPOSCAR`` is the supercell structure without atomic displacements, and vasprun0.xml is the result of DFT calculation for ``SPOSCAR``.
The ``--offset`` option subtract the offset (residual) components of atomic forces from the data in vasprun??.xml. 

Important::

    The ``--offset`` is optional, but we strongly recommend to use it when the fractional coordinates of atoms have degrees of freedom.

The ``extract.py`` can also parse the data from the output files of QE, OpenMX, xTAPP, and LAMMPS::

    **QE**
    ::

    $ python extract.py --QE=supercell.pw.in --offset=supercell.pw.out disp??.pw.out > DFSET

    **OpenMX**
    ::

    $ python extract.py --OpenMX=supercell.dat --offset=supercell.out disp??.out > DFSET

    **xTAPP**
    ::

    $ python extract.py --xTAPP=supercell.cg --offset=supercell.str disp??.str > DFSET


The LAMMPS case requires a special treatment. We first need to add the *dump* option in the LAMMPS input file as

    ::
    dump            1 all custom 1 XFSET id xu xy xz fx fy fz
    dump_modify     1 format float "%20.15f"

This option will generate the file *XFSET* which contains atomic coordinates and forces. 
After generating *XFSET* files for 10 structures and save them as *XFSET.01* ... *XFSET.10* , we can create ``DFSET`` as::

   **LAMMPS**
   :: 

   $ python extract.py --LAMMPS=supercell.lammps --offset=supercell.XFSET XFSET.?? > DFSET