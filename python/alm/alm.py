"""The main python file for ALM package."""

import warnings
from collections import OrderedDict

import numpy as np

from . import _alm as alm

atom_names = (
    "X",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Uut",
    "Uuq",
    "Uup",
    "Uuh",
    "Uus",
    "Uuo",
)

# From src/optimize.h
# {sparsefolver: str} is omitted because this is set at ALM.optimize.
# This order is not allowed to change because it is explicitly used in _alm.c.
optimizer_control_data_types = OrderedDict(
    [
        ("linear_model", int),
        ("use_sparse_solver", int),
        ("maxnum_iteration", int),
        ("tolerance_iteration", float),
        ("output_frequency", int),
        ("standardize", int),
        ("displacement_normalization_factor", float),
        ("debiase_after_l1opt", int),
        ("cross_validation", int),
        ("l1_alpha", float),
        ("l1_alpha_min", float),
        ("l1_alpha_max", float),
        ("num_l1_alpha", int),
        ("l1_ratio", float),
        ("save_solution_path", int),
        ("mirror_image_conv", int),
    ]
)


class ALM(object):
    """Calculate harmonic and anharmonic interatomic force constants.

    Attributes
    ----------
    lavec : ndarray
        Basis vectors. a, b, c are given as row vectors.
        shape=(3, 3), dtype='double'
    xcoord : ndarray
        Fractional coordinates of atomic points.
        shape=(num_atoms, 3), dtype='double'
    numbers : ndarray
        Atomic numbers.
        shape=(num_atoms,), dtype='intc'
    kind_names : OrderedDict
        Pairs of (atomic number, element name). Since the atomic number is the
        key of OrderedDict, only unique atomic numbers are stored and the
        order of ``numbers`` is preserved in the keys of this OrderedDict.
    displacements : ndarray
        Displacements of atoms in supercells used as training data.
        shape=(supercells, num_atoms, 3), dtype='double', order='C'
    forces : ndarray
        Forces of atoms in supercells used as training data.
        shape=(supercells, num_atoms, 3), dtype='double', order='C'
    verbosity : int
        Level of the output frequency either 0 (no output) or
        1 (normal output). Default is 0.
    output_filename_prefix : str
        More detailed logs are stored in files when this is given. This string
        is used to the prefix of filenames of logs.
    optimizer_control : dict
        Parameters to use elastic net regression.
    cv_l1_alpha : float (read-only)
        Alpha value to minimize fitting error of elastic net regression
        obtained by cross validation.

    """

    def __init__(self, lavec, xcoord, numbers, verbosity=0):
        """Initialize the alm object with the input crystal structure.

        Parameters
        ----------
        lavec : array_like
            Basis vectors. a, b, c are given as column vectors.
            shape=(3, 3), dtype='double'
        xcoord : array_like
            Fractional coordinates of atomic points.
            shape=(num_atoms, 3), dtype='double'
        numbers : array_like
            Atomic numbers.
            shape=(num_atoms,), dtype='intc'
        verbosity : int
            Level of the output frequency either 0 (no output) or
            1 (normal output). Default is 0.

        """
        self._id = None
        self._lavec = None
        self._xcoord = None
        self._numbers = None
        self._verbosity = False
        self._kind_names = None
        self._iconst = 11
        self._maxorder = 1

        self.lavec = lavec
        self.xcoord = xcoord
        self.numbers = numbers
        self._verbosity = verbosity

        self._output_filename_prefix = None

        # Whether python parameters are needed to be copied to C++ instance
        # or not.
        self._need_transfer = True

        # self.define() has been done or not.
        self._defined = False

    @property
    def lavec(self):
        """Getter of basis vectors.

        Returns
        -------
        lavec : ndarray
        Copy of basis vectors. a, b, c are given as row vectors.
        shape=(3, 3), dtype='double', order='C'

        """
        return np.array(self._lavec, dtype="double", order="C")

    @lavec.setter
    def lavec(self, lavec):
        """Setter of basis vectors.

        Parameters
        ----------
        lavec : array_like
        Basis vectors. a, b, c are given as row vectors.
        shape=(3, 3), dtype='double', order='C'

        """
        self._need_transfer = True
        self._lavec = np.array(lavec, dtype="double", order="C")

    @property
    def xcoord(self):
        """Getter of atomic point coordinates.

        Returns
        -------
        xcoord : ndarray
        Atomic point coordinates.
        shape=(num_atom, 3), dtype='double', order='C'

        """
        return np.array(self._xcoord, dtype="double", order="C")

    @xcoord.setter
    def xcoord(self, xcoord):
        """Setter of atomic point coordinates.

        Returns
        -------
        xcoord : ndarray
        Atomic point coordinates.
        shape=(num_atom, 3), dtype='double', order='C'

        """
        self._need_transfer = True
        self._xcoord = np.array(xcoord, dtype="double", order="C")

    @property
    def numbers(self):
        """Getter of atomic numbers.

        Returns
        -------
        numbers : ndarray
        Atomic numbers.
        shape=(num_atom,), dtype='intc', order='C'

        """
        return np.array(self._numbers, dtype="intc")

    @numbers.setter
    def numbers(self, numbers):
        """Setter of atomic numbers.

        Parameters
        ----------
        numbers : ndarray
        Atomic numbers.
        shape=(num_atom,), dtype='intc', order='C'

        """
        self._need_transfer = True
        self._numbers = np.array(numbers, dtype="intc")
        self._kind_names = OrderedDict.fromkeys(self._numbers)
        for key in self._kind_names:
            self._kind_names[key] = atom_names[key % 118]

    @property
    def kind_names(self):
        """Return the kind names."""
        return self._kind_names

    @property
    def verbosity(self):
        """Return the output verbosity."""
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity):
        """Set verbosity of output.

        Parameters
        ----------
        verbosity : int
            Choose the level of the output frequency from
            0 (no output) or 1 (normal output).

        """
        self._verbosity = verbosity
        self._set_verbosity()

    def set_verbosity(self, verbosity):
        """Set the output verbosity."""
        self.verbosity = verbosity

    def __enter__(self):
        """Run when entered in a with block."""
        self.alm_new()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Delete the pointer when exit from a with block."""
        self.alm_delete()

    def alm_new(self):
        """Create ALM instance in C++.

        This is also called by context manager when entering the block.

        ex.::

           with ALM(lavec, xcoord, numbers) as alm:

        Note
        ----
        When an ALM instance is created by ``alm_new``, it must be deleted
        by ``alm_delete`` to avoid memory leak.

        """
        if self._id is None:
            self._id = alm.alm_new()
            if self._id < 0:
                raise RuntimeError("Too many ALM objects")
            if self._verbosity is not None:
                self.verbosity = self._verbosity
        else:
            raise ("This ALM object is already initialized.")

    def alm_delete(self):
        """Delete ALM instance in C++.

        This is also called by context manager when exiting the block.

        ex.::

           with ALM(lavec, xcoord, numbers) as alm:

        """
        if self._id is None:
            self._show_error_not_initizalied()

        alm.alm_delete(self._id)
        self._id = None

    def suggest(self):
        """Compute displacement patterns to obtain force constants."""
        if self._id is None:
            self._show_error_not_initizalied()

        if self._defined:
            alm.suggest(self._id)
        else:
            self._show_error_not_defined()

    def optimize(self, solver="dense"):
        """Fit force constants to forces.

        Parameters
        ----------
        solver : str, default='dense'
            Solver choice for fitting either 'dense' or 'SimplicialLDLT'.

            - When solver='dense', the fitting is performed with the
              singular value decomposition implemented in LAPACK.
            - When solver='SimplicialLDLT', the fitting is performed with
              the sparse solver class SimplicialLDLT implemented in
              Eigen3 library.

        Returns
        -------
        info : int
            This tells condition how fitting went.
            0 if the fitting is successful, 1 otherwise.

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if not self._defined:
            self._show_error_not_defined()

        solvers = {"dense": "dense", "simplicialldlt": "SimplicialLDLT"}

        if solver.lower() not in solvers:
            msgs = [
                "The given solver option is not supported.",
                "Available options are 'dense' and 'SimplicialLDLT'.",
            ]
            raise ValueError("\n".join(msgs))

        info = alm.optimize(self._id, solvers[solver.lower()])

        return info

    @property
    def output_filename_prefix(self):
        """Return the prefix used for the output filenames."""
        return self._output_filename_prefix

    @output_filename_prefix.setter
    def output_filename_prefix(self, prefix):
        """Set output prefix of output filename."""
        if self._id is None:
            self._show_error_not_initizalied()

        if type(prefix) is str:
            self._output_filename_prefix = prefix
            self._set_output_filename_prefix(prefix)

    def _set_output_filename_prefix(self, prefix):
        self.output_filename_prefix = prefix

    @property
    def optimizer_control(self):
        """Return the outimizer_control that containts the options for optimization."""
        if self._id is None:
            self._show_error_not_initizalied()

        optctrl = alm.get_optimizer_control(self._id)
        keys = optimizer_control_data_types.keys()
        optcontrol = dict(zip(keys, optctrl))
        return optcontrol

    @optimizer_control.setter
    def optimizer_control(self, optcontrol):
        if self._id is None:
            self._show_error_not_initizalied()

        supported_keys = optimizer_control_data_types.keys()
        optctrl = []
        optcontrol_l = {key.lower(): optcontrol[key] for key in optcontrol}

        for i, key in enumerate(optcontrol):
            if key.lower() not in supported_keys:
                msg = "%s is not a valide key for optimizer control." % key
                raise KeyError(msg)

        for i, key in enumerate(supported_keys):
            if key in optcontrol_l:
                optctrl.append(optcontrol_l[key])
            else:
                optctrl.append(None)

        alm.set_optimizer_control(self._id, optctrl)

    def set_optimizer_control(self, optcontrol):
        """Set the options for force constant optimization."""
        self.optimizer_control = optcontrol

    @property
    def displacements(self):
        """Get displacements.

        Returns
        --------
        u : ndarray
            Atomic displacement patterns in supercells in Cartesian.
            shape=(supercells, num_atoms, 3), dtype='double', order='C'

        """
        if self._id is None:
            self._show_error_not_initizalied()

        ndata = alm.get_number_of_data(self._id)
        u = np.zeros((ndata, len(self._xcoord), 3), dtype="double", order="C")
        succeeded = alm.get_u_train(self._id, u)
        if succeeded:
            return u
        else:
            return None

    @displacements.setter
    def displacements(self, u):
        """Set displacements.

        Parameters
        ----------
        u : array_like
            Atomic displacement patterns in supercells in Cartesian.
            shape=(supercells, num_atoms, 3), dtype='double'

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if u.ndim != 3:
            msg = "Displacement array has to be three dimensions."
            raise RuntimeError(msg)

        alm.set_u_train(self._id, np.array(u, dtype="double", order="C"))

    @property
    def forces(self):
        """Get forces.

        Returns
        --------
        f : ndarray
            Forces in supercells.
            shape=(supercells, num_atoms, 3), dtype='double', order='C'

        """
        if self._id is None:
            self._show_error_not_initizalied()

        ndata = alm.get_number_of_data(self._id)
        f = np.zeros((ndata, len(self._xcoord), 3), dtype="double", order="C")
        succeeded = alm.get_f_train(self._id, f)
        if succeeded:
            return f
        else:
            return None

    @forces.setter
    def forces(self, f):
        """Set forces.

        Parameters
        ----------
        f : array_like
            Forces in supercells.
            shape=(supercells, num_atoms, 3), dtype='double'

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if f.ndim != 3:
            msg = "Force array has to be three dimensions."
            raise RuntimeError(msg)

        alm.set_f_train(self._id, np.array(f, dtype="double", order="C"))

    def set_training_data(self, u, f):
        """Set displacements and respective forces in supercell.

        Parameters
        ----------
        u : array_like
            Atomic displacement patterns in supercells in Cartesian.
            dtype='double'
            shape=(supercells, num_atoms, 3)
        f : array_like
            Forces in supercells.
            dtype='double'
            shape=(supercells, num_atoms, 3)

        """
        self.displacements = u
        self.forces = f

    def set_displacement_and_force(self, u, f):
        """Set displacements and forces. Deprecated version."""
        warnings.warn(
            "set_displacement_and_force is deprecated. " "Use set_training_data.",
            DeprecationWarning,
        )

        self.set_training_data(u, f)

    def define(
            self, maxorder, cutoff_radii=None,
            nbody=None, symmetrization_basis="Lattice"
    ):
        """Define the Taylor expansion potential.

        Parameters
        ----------
        maxorder : int
            Maximum order of the Taylor expansion potential.

            - If ``maxorder = 1``, only harmonic (2nd-order) terms are
              considered.
            - If ``maxorder = 2``, both harmonic and cubic terms are
              considered.

        cutoff_radii : array_like, default = None
            Cutoff radii defined for each order.
            When a negative value is provided, the cutoff radius is not used.
            dtype='double'
            shape=(maxorder, num_elems, num_elems)

        nbody : array_like, default = None
            Option to neglect multi-body interactions.
            dtype='intc'
            shape=(maxorder,)

        symmetrization_basis : str, default='Lattice'
            Either 'Cartesian' or 'Lattice'. Symmetrization of force constants
            is done either in the matrix based on crystal coordinates
            ('Lattice') or Cartesian coordinates ('Cartesian').

        """
        if self._id is None:
            self._show_error_not_initizalied()

        self._transfer_parameters()

        if nbody is None:
            nbody = []
            for i in range(maxorder):
                nbody.append(i + 2)

        else:
            if len(nbody) != maxorder:
                msg = "The size of nbody must be equal to maxorder."
                raise RuntimeError(msg)

        if cutoff_radii is None:
            _cutoff_radii = None
        else:
            _cutoff_radii = np.array(cutoff_radii, dtype="double", order="C")
            nelem = len(_cutoff_radii.ravel())
            if (nelem // maxorder) * maxorder != nelem:
                msg = "The array shape of cutoff_radii is wrong."
                raise RuntimeError(msg)
            nkd = int(round(np.sqrt(nelem // maxorder)))
            if nkd ** 2 - nelem // maxorder != 0:
                msg = "The array shape of cutoff_radii is wrong."
                raise RuntimeError(msg)
            _cutoff_radii = np.reshape(_cutoff_radii, (maxorder, nkd, nkd), order="C")

        self._maxorder = maxorder

        if symmetrization_basis.lower() in ["lattice", "cartesian"]:
            fc_basis = symmetrization_basis.capitalize()
        else:
            fc_basis = "Lattice"

        alm.define(
            self._id, maxorder, np.array(nbody, dtype="intc"), _cutoff_radii, fc_basis
        )

        alm.init_fc_table(self._id)

        self._defined = True

    def set_constraint(self, translation=True, rotation=False):
        """Set constraints for the translational and rotational invariances.

        Parameters
        ----------
        translation : bool, optional (default = True)
            When set to ``True``, the translational invariance
            (aka acoustic sum rule) is imposed between force constants.

        rotation : bool, optional (default = False)
            When set to ``True``, the rotational invariance is imposed between
            force constants. This function is not implemented.

        """
        if rotation is True:
            raise ("Rotational invariance is not supported in python API.")

        if translation is True:
            iconst = 11
        else:
            iconst = 10

        self._iconst = iconst
        alm.set_constraint_type(self._id, self._iconst)

    def freeze_fc(self, fc_values, fc_indices):
        """Freeze force to the given values upon optimization.

        Parameters
        ----------
        fc_values: array_like, dtype='float', shape=(num_fc)
            Array of force constant values.

        fc_indices: array_like, dtype='int', shape=(num_fc, fc_order + 1)
            Array of flattened indices 3 * index_atom + index_xyz of the
            force constants.

        Note
        ----
        The order of force constant is judged using the dimension of
        the input `fc_indices`.

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if fc_indices.ndim != 2:
            msg = "fc_indices array has to be two dimensions."
            raise RuntimeError(msg)

        if len(fc_indices) != len(fc_values):
            msg = "length of fc_indices and fc_values must be the same."
            raise RuntimeError(msg)

        fc_indices = np.array(fc_indices, dtype="intc", order="C")
        fc_values = np.array(fc_values, dtype="float")
        alm.set_fcs_freeze(self._id, fc_values, fc_indices)

    def getmap_primitive_to_supercell(self):
        """Return the mapping information from the primitive cell to the supercell.

        Returns
        -------
        map_p2s : array_like
            The mapping information of atoms from the primitive cell to the
            supercell.
            dtype='intc'
            shape = (num_trans, num_atoms_primitive)

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if not self._defined:
            self._show_error_not_defined()

        map_p2s = np.zeros(len(self._xcoord), dtype="intc")
        ntrans = alm.get_atom_mapping_by_pure_translations(self._id, map_p2s)
        return map_p2s.reshape((ntrans, -1))

    def get_displacement_patterns(self, fc_order):
        """Return the displacement patterns to obtain force constants.

        Parameters
        ----------
        fc_order : int
            The order of force constants to get the displacement patterns.

            - If ``fc_order = 1``, returns patterns for harmonic force
              constants.
            - If ``fc_order = 2``, returns patterns for cubic force constants.
            - If ``fc_order = 3``, returns patterns for quartic force
              constants.
            - ...

        Returns
        -------
        all_disps : array_like, shape = (n_patterns,)
            The array of tuples (``atom_index``, ``direction``, ``basis``),
            where ``direction`` is the numpy.ndarray of size = (3,)
            representing the direction of the displacement,
            and ``basis`` is a string either "Cartesian" or "Fractional".

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if fc_order > self._maxorder:
            msg = (
                "The fc_order must not be larger than the maximum order " "(maxorder)."
            )
            raise ValueError(msg)

        numbers = self._get_number_of_displaced_atoms(fc_order)
        tot_num = np.sum(numbers)
        atom_indices = np.zeros(tot_num, dtype="intc")
        disp_patterns = np.zeros((tot_num, 3), dtype="double", order="C")
        nbasis = alm.get_displacement_patterns(
            self._id, atom_indices, disp_patterns, fc_order
        )
        basis = ["Cartesian", "Fractional"][nbasis]
        all_disps = []
        pos = 0
        for num in numbers:
            disp = []
            for i in range(num):
                disp.append((atom_indices[pos], disp_patterns[pos], basis))
                pos += 1
            all_disps.append(disp)
        return all_disps

    def get_fc(self, fc_order, mode="origin", permutation=True):
        """Return the force constant values.

        Parameters
        ----------
        fc_order : int
            The order of force constants to get.

            - If ``fc_order = 1``, returns harmonic force constants.
            - If ``fc_order = 2``, returns cubic force constants.
            - If ``fc_order = 3``, returns quartic force constants.
            - ...

        mode : str, optional (default="origin")
            The choice of the force constant list to be returned.

            - If "origin", returns the reducible set of force constants,
              whose first element corresponds to an atom in the
              primitive cell at the origin.
            - If "all", returns the all non-zero elements of force constants
              in the supercell.
            - If "irreducible" or "irred", returns the irreducible set of
              force constants.

        permutation : bool (default=True)
            The flag for printing out elements with permutation symmetry.
            Effective only when ``mode = origin`` or ``mode = all``.

            - If True, returns force constants after replicating elements
              by the permutation of indices.
            - If False, returns force constants without replicating elements
              by the permutation of indices. For "origin" and "all", all
              indices except the first index participate to the permutation
              of indices to reduce the number of the output values.

        Returns
        -------
        fc_values : array_like, dtype='double', shape=(num_fc,)
            Force constant values.

        elem_indices : array_like, dtype='int', shape=(num_fc, fc_order + 1)
            Array of flattened indices 3 * index_atom + index_xyz.

        Note
        ----
        This method returns force constants in Cartesian basis
        when ``mode = origin`` and ``mode = all`.
        When ``mode = irred``, it returns the irreducible set of
        force constants in the basis defined via "symmetrization_basis"
        of the alm.define method.

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if fc_order > self._maxorder:
            msg = (
                "The fc_order must not be larger than the maximum order " "(maxorder)."
            )
            raise ValueError(msg)

        perm_int = permutation * 1

        if mode == "origin":

            fc_length = self._get_number_of_fc_origin(fc_order, perm_int)
            fc_values = np.zeros(fc_length, dtype="double")
            elem_indices = np.zeros((fc_length, fc_order + 1), dtype="intc", order="C")

            alm.get_fc_origin(self._id, fc_values, elem_indices, perm_int)

            return fc_values, elem_indices

        elif mode == "irreducible" or mode == "irred":

            fc_length = self._get_number_of_irred_fc_elements(fc_order)
            fc_values = np.zeros(fc_length, dtype="double")
            elem_indices = np.zeros((fc_length, fc_order + 1), dtype="intc", order="C")

            alm.get_fc_irreducible(self._id, fc_values, elem_indices)

            return fc_values, elem_indices

        elif mode == "all":

            map_p2s = np.zeros(len(self._xcoord), dtype="intc")
            ntrans = alm.get_atom_mapping_by_pure_translations(self._id, map_p2s)
            fc_length = self._get_number_of_fc_origin(fc_order, perm_int) * ntrans
            fc_values = np.zeros(fc_length, dtype="double")
            elem_indices = np.zeros((fc_length, fc_order + 1), dtype="intc", order="C")

            alm.get_fc_all(self._id, fc_values, elem_indices, perm_int)

            return fc_values, elem_indices

        else:
            raise ValueError("Invalid mode in get_fc.")

    def set_fc(self, fc_in):
        """Copy force constant obtained by an external optimizer to the ALM instance.

        Parameters
        ----------
        fc_in : array_like
            The irreducible set of force constants.
            dtype='double'
            shape=(num_fc,)

        Note
        ----
        When an external optimizer, such as numpy.linalg.lstsq, is used to fit
        force constants, the force constants need to be passed to
        the ALM instance by ``set_fc`` to use the ``get_fc`` method.

        """
        if self._id is None:
            self._show_error_not_initizalied()

        maxorder = self._maxorder
        fc_length_irred = 0
        for i in range(maxorder):
            fc_length_irred += self._get_number_of_irred_fc_elements(i + 1)

        if fc_length_irred != len(fc_in):
            msg = "The size of the given force constant array is incorrect."
            raise RuntimeError(msg)

        alm.set_fc(self._id, np.array(fc_in, dtype="double", order="C"))

    def get_matrix_elements(self):
        """Return the sensing matrix A and force vector b.

        Returns
        -------
        amat : ndarray, dtype='double'
            shape=(3 * num_atoms * ndata_training, num_fc_irred), order='F'.
            The sensing matrix A calculated from the displacements.

        bvec : ndarray, dtype='double'
            shape=(3 * num_atoms * ndata_training,)
            The vector b calculated from the atomic forces.


        Note
        ----
        From the amat (``A``) and bvec (``b``), the force constant vector ``x``
        can be obtained by solving the least-square problem:
        x = argmin_{x} | Ax-b|^{2}.

        """
        if self._id is None:
            self._show_error_not_initizalied()

        maxorder = self._maxorder
        nrows = self._get_nrows_amat()

        fc_length = 0
        for i in range(maxorder):
            fc_length += self._get_number_of_irred_fc_elements(i + 1)

        amat = np.zeros(nrows * fc_length, dtype="double", order="C")
        bvec = np.zeros(nrows, dtype="double")
        alm.get_matrix_elements(self._id, amat, bvec)

        return (np.reshape(amat, (nrows, fc_length), order="F"), bvec)

    def save_fc(self, filename, format="alamode", maxorder_to_save=3):
        """Save the fitted force constants to a file with a given format and filename.

        Parameters
        ----------
        filename : str
            The filename where the force constants are saved.

        format : str, optional (default="alamode")
            The file format to save the force constants

            - If "alamode", save force constants in the ALAMODE XML format.
              This file can be read by the ANPHON code of ALAMODE.
            - If "shengbte" or "thirdorderpy", save third-order force constants
              in the thirdorder.py format only when the model includes the
              third-order terms. This file can be read by the ShengBTE code.
            - If "qe-fc", save harmonic force constants in the Quantum-ESPRESSO (QE)
              fc format. The file can be processed by QE matdyn.x or other third-party
              softwares such as EPW.

        maxorder_to_save : int, optional (default=3)
            Defines the maximum order of force constants to save in the file.

            Effective only when format="alamode".
            The default value is maxorder_to_save=3, which means that the
            force constants up to the fourth-order (order = 3) are saved
            and higher-order terms are omitted.
            The default setting should be adequate for most
            cases (thermal conductivity, self-consistent phonon calculations).

        """
        if self._id is None:
            self._show_error_not_initizalied()

        if filename is None:
            raise ValueError("filename is invalid.")

        if format == "alamode":
            alm.save_fc(self._id, filename, format, maxorder_to_save)

        elif format == "shengbte" or format == "thirdorderpy":
            if self._maxorder < 2:
                msg = (
                    "The current force constant model does not have"
                    " third-order terms"
                )
                raise ValueError(msg)
            else:
                alm.save_fc(self._id, filename, "shengbte", maxorder_to_save)

        elif format == "qe-fc":
            alm.save_fc(self._id, filename, "qefc", maxorder_to_save)

        else:
            raise ValueError("Invalid format in save_fc.")

    def get_fc_dependency(self, fc_order):
        """Return the force constant dependency information as dictionary.

        Parameters
        ----------
        fc_order : int
            The order of force constants to get.

            - If ``fc_order = 1``, returns harmonic force constants.
            - If ``fc_order = 2``, returns cubic force constants.
            - If ``fc_order = 3``, returns quartic force constants.
            - ...

        Returns
        -------
        dependency : dictionary
            Force constant dependency information

        Note
        ----
        The returned dictionary contains the information of equivalent
        force constants as well as the equality constraints between
        force constants that needs to be hold to satisfy the
        translational invariance (ASR, acoustic sum rule).
        If one want to exclude the constraints due to the ASR,
        please set alm.set_constraint(translation=False)
        before calling this method.
        """
        n_fc_irred = self._get_number_of_irred_fc_elements(fc_order)
        n_fc_origin = self._get_number_of_fc_origin(fc_order, 0)
        elem_indices_irred = np.zeros((n_fc_irred, fc_order + 1),
                                      dtype="intc", order="C")
        elem_indices_origin = np.zeros((n_fc_origin, fc_order + 1),
                                       dtype="intc", order="C")
        dependency_mat = np.zeros(n_fc_irred * n_fc_origin, dtype="double", order="C")

        alm.get_fc_dependency(self._id, elem_indices_irred,
                              elem_indices_origin, dependency_mat)
        xyz = "xyz"
        dependency_mat = np.reshape(dependency_mat, (n_fc_origin, n_fc_irred))
        nonzeros = np.where(np.abs(dependency_mat) > 1.0e-5)

        dependency_dic = {}
        i_prev = -1
        dic_tmp = {}
        for i, j in zip(nonzeros[0], nonzeros[1]):
            if i == i_prev:
                key = ''
                for entry in elem_indices_irred[j]:
                    key += '{}{}'.format(entry // 3 + 1, xyz[entry % 3])
                dic_tmp[key] = dependency_mat[i, j]
            else:
                if len(dic_tmp) > 0:
                    key = ''
                    for entry in elem_indices_origin[i_prev]:
                        key += '{}{}'.format(entry // 3 + 1, xyz[entry % 3])
                    dependency_dic[key] = dic_tmp
                dic_tmp = {}
                key = ''
                for entry in elem_indices_irred[j]:
                    key += '{}{}'.format(entry // 3 + 1, xyz[entry % 3])
                dic_tmp[key] = dependency_mat[i, j]
                i_prev = i

        if len(dic_tmp) > 0:
            key = ''
            for entry in elem_indices_origin[i_prev]:
                key += '{}{}'.format(entry // 3 + 1, xyz[entry % 3])
            dependency_dic[key] = dic_tmp

        return dependency_dic

    def get_cv_l1_alpha(self):
        """Return L1 alpha at minimum CV."""
        if self._id is None:
            self._show_error_not_initizalied()

        return alm.get_cv_l1_alpha(self._id)

    def _transfer_parameters(self):
        if self._need_transfer:
            self._set_cell()
            self._need_transfer = False
            self._defined = False

    def _set_cell(self):
        """Inject crystal structure in C++ instance."""
        if self._id is None:
            self._show_error_not_initizalied()

        if self._lavec is None:
            msg = "Basis vectors are not set."
            raise RuntimeError(msg)

        if self._xcoord is None:
            msg = "Atomic point coordinates (positions) are not set."
            raise RuntimeError(msg)

        if self._numbers is None:
            msg = "Atomic numbers are not set."
            raise RuntimeError(msg)

        if len(self._xcoord) != len(self._numbers):
            msg = "Numbers of atomic points and atomic numbers don't agree."
            raise RuntimeError(msg)

        kind_numbers = np.array(list(self._kind_names.keys()), dtype="intc")
        alm.set_cell(self._id, self._lavec, self._xcoord, self._numbers, kind_numbers)

    def _set_verbosity(self):
        """Inject verbosity in C++ instance."""
        if self._id is None:
            self._show_error_not_initizalied()

        alm.set_verbosity(self._id, self._verbosity)

    def _get_nrows_amat(self):
        """Private method to return the number of training data sets."""
        if self._id is None:
            self._show_error_not_initizalied()
        nrows_amat = alm.get_nrows_amat(self._id)
        return nrows_amat

    def _get_id(self):
        """Private method to return the instance ID."""
        return self._id

    def _get_number_of_displacement_patterns(self, fc_order):
        """Private method to return the number of displacement patterns.

        Parameters
        ----------
        fc_order : int
            The order of force constants.
            fc_order = 1 for harmonic, fc_order = 2 for cubic ...

        """
        return alm.get_number_of_displacement_patterns(self._id, fc_order)

    def _get_number_of_displaced_atoms(self, fc_order):
        """Private method to return the number of displaced atoms.

        Parameters
        ----------
        fc_order : int
            The order of force constants.
            fc_order = 1 for harmonic, fc_order = 2 for cubic ...

        """
        num_disp_patterns = self._get_number_of_displacement_patterns(fc_order)
        numbers = np.zeros(num_disp_patterns, dtype="intc")
        alm.get_number_of_displaced_atoms(self._id, numbers, fc_order)
        return numbers

    def _get_number_of_fc_elements(self, fc_order):
        """Private method to get the number of force constants.

        Parameters
        ----------
        fc_order : int
            The order of force constants.
            fc_order = 1 for harmonic, fc_order = 2 for cubic ...

        """
        return alm.get_number_of_fc_elements(self._id, fc_order)

    def _get_number_of_fc_origin(self, fc_order, permutation):
        """Private method to get the number of force constants for fc_origin.

        Parameters
        ----------
        fc_order : int
            The order of force constants.
            fc_order = 1 for harmonic, fc_order = 2 for cubic ...

        permutation: int
            Flag to include permuted elements
            permutation = 0 for skipping permuted elements,
            permutation = 1 for including them

        """
        return alm.get_number_of_fc_origin(self._id, fc_order, permutation)

    def _get_number_of_irred_fc_elements(self, fc_order):
        """Private method to get the number of irreducible set of force constants.

        Parameters
        ----------
        fc_order : int
            The order of force constants.
            fc_order = 1 for harmonic, fc_order = 2 for cubic ...

        """
        return alm.get_number_of_irred_fc_elements(self._id, fc_order)

    def _show_error_not_initizalied(self):
        """Private method to raise an error."""
        msg = (
            "This ALM instance has to be initialized by ALM.alm_new() or "
            "context manager."
        )
        raise RuntimeError(msg)

    def _show_error_not_defined(self):
        msg = "This ALM.define() has to be done beforehand."
        raise RuntimeError(msg)
