.. _python_module:

Using ALM from python
=====================

ALM's python module can be made following :ref:`this page
<compile_with_conda_packages>`.

The job of ALM is to select force constants elements among all
possible elements for atoms in a supercell and then ALM fits relations
between displacements and forces to those force constants
elements. The selection of the force constants elements is done by the
order of force constants, force constants symmetry, user-input cutoff
distances, and maybe LASSO regression analysis. In the following, how
to use ALM is presented step by step.

.. contents::
   :depth: 2
   :local:

Initialization
--------------

ALM class instance is made by context manager as follows::

   from alm import ALM

   with ALM(lavec, xcoord, numbers) as alm:
       ...

``lavec``, ``xcoord``, and ``numbers`` are the essential parameters
and are the basis vectors, fractional coordinates of atoms, and atomic
numbers, respectively. ``lavec`` is the :math:`3 \times 3` matrix and
:math:`a`, :math:`b`, :math:`c` basis vectors are given as the row
vectors, i.e.,

.. math::

   \begin{pmatrix}
   a_x & a_y & a_z \\
   b_x & b_y & b_z \\
   c_x & c_y & c_z
   \end{pmatrix}.

``xcoord`` is the :math:`n\times 3` matrix. Each row gives the point
coordinates of the atom. ``numbers`` is a list of integer numbers.

Dataset: displacements and forces
---------------------------------

ALM requires a dataset composed of sets of pairs of atomic
displacements from the input crystal structure and forces due to the
displacements.

Because there should be many supercells with different
configurations of displacements, the array shape of the displacement
data is ``(number of supercells, number of atoms, 3)``. The array
shape of forces is the same as that of the displadements.

For example, the file format of datasets for ALM commandline interface
(``DFSET``) can be easily tranformed to the displacements and forces by

::

   dfset = np.loadtxt("DFSET").reshape((-1, number_of_atoms, 6))
   displacements = dfset[:, :, :3]
   forces = dfset[:, :, 3:]

These data are set to ALM by

::

   from alm import ALM

   with ALM(lavec, xcoord, numbers) as alm:
       ...
       alm.displacements = displacements
       alm.forces = forces

Selection of force constants elements
-------------------------------------

The basis force constants selection is performed by

::

    alm.define(maxorder, cutoff_radii, nbody)

``maxorder = 1`` for only harmonic (2nd-order) force constants, and
``maxorder = 2`` for 2nd- and 3rd-order force constants, and so on,
i.e. up to (n+1)th order force constants are included in the
consideration with ``maxorder=n``. ``cutoff_radii`` controls how far
the atomic pairs are searched to select force constants
elements. ``nbody`` is used to limit interaction of atoms.  The
details meanings are found in :ref:`interaction_field` for
``maxorder`` and ``nbody`` and in :ref:`cutoff_field` for
``cutoff_radii``.

Force constants calculation
---------------------------

Force constants are calculated by fitting dataset of displacements and
forces to the force constants elements selected by

::

    alm.optimize()

There are two solvers, which are chosen either by ``solver='dense'``
(default) or ``solver='SimplicialLDLT'`` as the keyword argument.

Extraction of force constants values
------------------------------------

The calculated force constants elements are stored in ALM intance in
the ALM's manner. They are extracted by

::

   alm.get_fc(fc_order, mode)

By ``fc_order=n``, (n+1)th order force constants are
extracted. ``mode`` chooses the format of force constants. There are
three modes, but two of them are important, which are ``all`` and
``origin``. By ``mode=all``, all elements of force constants are
returned except for the elements whose values are 0. With
``mode=origin``, the first atomic indices of force constants are
limited for only those in the primitive cell.

LASSO and elastic net regression
--------------------------------

As shown in :ref:`optimize_field`, ALM has a functionality to compute
force constants using LASSO and elestic net regression. This feature
is accessed from ALM python module such as by

::

   optcontrol = {'linear_model': 2,
                 'cross_validation': 4,
                 'num_l1_alpha': 50}
   alm.optimizer_control = optcontrol

The controllable parameters and their variable types are listed as
follows::

   optimizer_control_data_types = OrderedDict([
       ('linear_model', int),                         # LMODEL
       ('use_sparse_solver', int),                    # '=1' equivalent to SimplicialLDLT
       ('maxnum_iteration', int),                     # MAXITER
       ('tolerance_iteration', float),                # CONV_TOL
       ('output_frequency', int),                     # NWRITE
       ('standardize', int),                          # STANDARDIZE
       ('displacement_normalization_factor', float),  # ENET_DNORM
       ('debiase_after_l1opt', int),                  # DEBIAS_OLS
       ('cross_validation', int),                     # CV
       ('l1_alpha', float),                           # L1_ALPHA
       ('l1_alpha_min', float),                       # CV_MINALPHA
       ('l1_alpha_max', float),                       # CV_MAXALPHA
       ('num_l1_alpha', int),                         # CV_NALPHA
       ('l1_ratio', float),                           # L1_RATIO
       ('save_solution_path', int)])                  # SOLUTION_PATH


Wrap-up and example
-------------------

Some examples are found in ``example`` directory.

The following python script is an example to computer 2nd and 3rd
order force constants by ordinary least square fitting. Let's assume
the crystal is wurtzite-type AlN as found in the example
directory. For the 2nd order, no cutoff radii are used, but for the
3rd order, cutoff radii of 4 Angstrom is chosen between all pairs of
atoms (Al-Al, Al-N, N-N). The fitting is achieved by using lapack SVD
solver. Be sure that the memory usage is ~1.6GB and the whole
calculation takes a few minutes, depending on computers though.

::

   import h5py
   import numpy as np
   from alm import ALM

   with open("POSCAR_AlN") as f:
       lines = f.readlines()
   [lines.pop(0) for i in range(2)]
   [lines.pop(3) for i in range(3)]
   vals = np.array([np.fromstring(l, dtype='double', sep=' ') for l in lines])
   lavec = vals[:3]
   xcoord = vals[3:]
   numbers = [13, ] * 36 + [7, ] * 36
   natom = len(numbers)

   dfset = np.loadtxt("DFSET_AlN").reshape((-1, natom, 6))
   displacements = dfset[:, :, :3]
   forces = dfset[:, :, 3:]

   cutoff_radii = [np.ones((2, 2)) * -1, np.ones((2, 2)) * 4]

   with ALM(lavec, xcoord, numbers, verbosity=1) as alm:
       alm.displacements = displacements
       alm.forces = forces
       alm.define(2, cutoff_radii=cutoff_radii)
       alm.optimize()

       fc2 = np.zeros((natom, natom, 3, 3), dtype='double', order='C')
       fc3 = np.zeros((natom, natom, natom, 3, 3, 3), dtype='double', order='C')
       for fc, indices in zip(*alm.get_fc(1, mode='all')):
           v1, v2 = indices // 3
           c1, c2 = indices % 3
           fc2[v1, v2, c1, c2] = fc
       for fc, indices in zip(*alm.get_fc(2, mode='all')):
           v1, v2, v3 = indices // 3
           c1, c2, c3 = indices % 3
           fc3[v1, v2, v3, c1, c2, c3] = fc
       with h5py.File('fc.hdf5', 'w') as w:
           w.create_dataset('fc2', data=fc2, compression='gzip')
           w.create_dataset('fc3', data=fc3, compression='gzip')
