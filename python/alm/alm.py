import numpy as np
from . import _alm as alm
import warnings


class ALM:
    """Calculate harmonic and anharmonic interatomic force constants

    Attributes
    ----------

    """

    def __init__(self, lavec, xcoord, kd):
        """

        Parameters
        ----------
        lavec : array_like
            Basis vectors. a, b, c are given as column vectors.
            shape=(3, 3)
            dtype='double'
        xcoord : array_like
            Fractional coordinates of atomic points.
            shape=(num_atoms, 3)
            dtype='double'
        kd : array_like
            Atomic numbers.
            shape=(num_atoms,)
            dtype='intc'

        """

        self._id = None
        self._lavec = np.array(lavec, dtype='double', order='C')
        self._xcoord = np.array(xcoord, dtype='double', order='C')
        self._kd = np.array(kd, dtype='intc', order='C')
        self._iconst = 11
        self._verbosity = 0
        self._norder = 1

    def __enter__(self):
        self.alm_new()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.alm_delete()

    def alm_new(self):
        """Create ALM instance in C++.

        This is also called by context manager when entering the block.

        ex.::

           with ALM(lavec, xcoord, kd) as alm:

        """

        if self._id is None:
            self._id = alm.alm_new()
            if self._id < 0:
                print("Too many ALM objects")
                raise RuntimeError
            self._set_cell()
            self._set_verbosity()
        else:
            print("This ALM object is already initialized.")
            raise

    def alm_delete(self):
        """Delete ALM instance in C++.

        This is also called by context manager when exiting the block.

        ex.::

           with ALM(lavec, xcoord, kd) as alm:

        """

        if self._id is None:
            self._show_error_message()

        alm.alm_delete(self._id)
        self._id = None

    def suggest(self):
        """Suggest displacement patterns to obtain force constants."""

        if self._id is None:
            self._show_error_message()

        alm.run_suggest(self._id)

    def optimize(self, solver='dense'):
        """Fit force constants to forces.

        Parameters
        ----------
        solver : str, default='dense'
            Solver choice for fitting either 'dense' or 'SimplicialLDLT'.
            When solver='dense', the fitting is performed with the
            singular value decomposition implemented in LAPACK.
            When solver='SimplicialLDLT', the fitting is performed with
            the sparse solver class SimplicialLDLT implemented in
            Eigen3 library.

        Returns
        -------
        info : int
            This tells condition how fitting went.
            0 if the fitting is successful, 1 otherwise.

        """

        if self._id is None:
            self._show_error_message()

        if solver not in ['dense', 'SimplicialLDLT']:
            print("The given solver option is not supported.")
            print("Available options are 'dense' and 'SimplicialLDLT'.")
            raise ValueError

        info = alm.optimize(self._id, solver)

        return info

    def set_displacement_and_force(self, u, f):
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

        if self._id is None:
            self._show_error_message()

        alm.set_displacement_and_force(
            self._id,
            np.array(u, dtype='double', order='C'),
            np.array(f, dtype='double', order='C'))

    def define(self, norder, rcs, nbody=None):
        """Define the Taylor expansion potential.

        Parameters
        ----------
        norder : int
            Maximum order of the Taylor expansion potential.
            When norder = 1, only harmonic (2nd-order) terms are considered and
            when norder = 2, both harmonic and cubic terms are considered.

        rcs : array_like
            Cutoff radii defined for each order.
            When a negative value is provided, the cutoff radius is not used.
            dtype='double'
            shape=(norder, num_elems, num_elems)

        nbody : array_like, default = None
            Option to neglect multi-body interactions.
            dtype='int'
            shape=(norder,)

        """

        if self._id is None:
            self._show_error_message()

        if nbody is None:
            nbody = []
            for i in range(norder):
                nbody.append(i + 2)

        else:
            if len(nbody) != norder:
                print("The size of nbody must be equal to norder.")
                raise RuntimeError

        self._norder = norder
        self._set_norder()
        alm.set_cutoff_radii(self._id,
                             np.array(rcs, dtype='double', order='C'))
        alm.set_nbody_rule(self._id, np.array(nbody, dtype='intc'))

        alm.generate_force_constant(self._id)

    def set_constraint(self, translation=True, rotation=False):
        """Set constraints for the translational and rotational invariances

        Parameters
        ----------
        translation : bool, optional (default = True)
            When set to ``True``, the translational invariance (aka acoustic sum rule)
            is imposed between force constants.

        rotation : bool, optional (default = False)
            When set to ``True``, the rotational invariance is imposed between
            force constants. This function is not implemented.

        """

        if rotation is True:
            print("Rotational invariance is not supported in API.")
            raise 

        if translation is True:
            iconst = 11
        else:
            iconst = 10

        self._iconst = iconst
        alm.set_constraint_type(self._id, self._iconst)

    def set_verbosity(self, verbosity):
        """Set verbosity of output.

        Parameters
        ----------
        verbosity : int
            Choose the level of the output frequency from 0 (no output) or 1 (normal output).

        """

        if self._id is None:
            self._show_error_message()

        self._verbosity = verbosity
        alm.set_verbosity(self._id, self._verbosity)

    def _set_verbosity(self):
        """Set verbosity of output. Do the same job as the public function set_verbosity."""

        if self._id is None:
            self._show_error_message()

        alm.set_verbosity(self._id, self._verbosity)

    def get_atom_mapping_by_pure_translations(self):
        """Returns the mapping information from the primitive cell to the supercell.

        Returns
        -------
        map_p2s : array_like
            The mapping information of atoms from the primitive cell to the supercell.
            dtype='int'
            shape = (num_trans, num_atoms_primitive)

        """
        
        if self._id is None:
            self._show_error_message()

        map_p2s = np.zeros(len(self._xcoord), dtype='intc')
        ntrans = alm.get_atom_mapping_by_pure_translations(self._id, map_p2s)
        return map_p2s.reshape((ntrans, -1))

    def get_displacement_patterns(self, fc_order):
        """Returns the displacement patterns to obtain force constants.

        Parameters
        ----------
        fc_order : int
            The order of force constants to get the displacement patterns.
            fc_order = 1 means harmonic, and fc_order = 2 means cubic, etc.

        Returns
        -------
        all_disps : array_like, shape = (n_patterns,)
            The array of tuples (``atom_index``, ``direction``, ``basis``),
            where ``direction`` is the numpy.ndarray of size = (3,) representing
            the direction of the displacement, and ``basis`` is a string
            either "Cartesian" or "Fractional".

        """

        if self._id is None:
            self._show_error_message()

        numbers = self._get_numbers_of_displacements(fc_order)
        tot_num = np.sum(numbers)
        atom_indices = np.zeros(tot_num, dtype='intc')
        disp_patterns = np.zeros((tot_num, 3), dtype='double', order='C')
        nbasis = alm.get_displacement_patterns(self._id,
                                               atom_indices,
                                               disp_patterns,
                                               fc_order)
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

    def get_fc(self, fc_order, mode="origin"):
        """Returns the force constant values

        Parameters
        ----------
        fc_order : int
            The order of force constants to get the force constants.
            fc_order = 1 means harmonic, and fc_order = 2 means cubic, etc.

        mode : str, optional (default="origin")
            The choice of the force constant list to be returned.

            - If "origin", returns the reducible set of force constants,
              whose first element corresponds to an atom in the 
              primitive cell at the origin.
            - If "irreducible" or "irred", returns the irreducible set of
              force constants.
            - If "all", returns the all non-zero elements of force constants
              in the supercell.

        Returns
        -------
        fc_values : array_like, dtype='double', shape=(num_fc,)
            Force constant values.

        elem_indices : array_like, dtype='int', shape=(num_fc, fc_order + 1)
            Array of flattened indices 3 * index_atom + index_xyz.

        Note
        ----
        This method does not replicate force constants elements that
        can be generated from the set by the permutation of indices.

        """

        if self._id is None:
            self._show_error_message()

        if mode == "origin":

            fc_length = self._get_number_of_fc_elements(fc_order)
            fc_values = np.zeros(fc_length, dtype='double')
            elem_indices = np.zeros((fc_length, fc_order + 1),
                                    dtype='intc', order='C')

            alm.get_fc_origin(self._id, fc_values, elem_indices)

            return fc_values, elem_indices

        elif mode == "irreducible" or mode == "irred":

            fc_length = self._get_number_of_irred_fc_elements(fc_order)
            fc_values = np.zeros(fc_length, dtype='double')
            elem_indices = np.zeros((fc_length, fc_order + 1),
                                    dtype='intc', order='C')

            alm.get_fc_irreducible(self._id, fc_values, elem_indices)

            return fc_values, elem_indices

        elif mode == "all":

            map_p2s = np.zeros(len(self._xcoord), dtype='intc')
            ntrans = alm.get_atom_mapping_by_pure_translations(self._id,
                                                               map_p2s)
            fc_length = self._get_number_of_fc_elements(fc_order) * ntrans
            fc_values = np.zeros(fc_length, dtype='double')
            elem_indices = np.zeros((fc_length, fc_order + 1),
                                    dtype='intc', order='C')

            alm.get_fc_all(self._id, fc_values, elem_indices)

            return fc_values, elem_indices

        else:
            print("Invalid mode in get_fc.")
            raise ValueError


    def set_fc(self, fc_in):
        """Copy force constant obtained by an external optimizer to the ALM instance.

        Parameters
        ----------
        fc_in : array_like
            The irreducible set of force constants.
            dtype='double'
            shape=(num_fc,)

        """

        if self._id is None:
            self._show_error_message()

        norder = self._norder
        fc_length_irred = 0
        for i in range(norder):
            fc_length_irred += self._get_number_of_irred_fc_elements(i + 1)

        if fc_length_irred != len(fc_in):
            print("The size of the given force constant array is incorrect.")
            raise RuntimeError

        alm.set_fc(self._id, np.array(fc_in, dtype='double', order='C'))

    def get_matrix_elements(self):
        """Returns the sensing matrix A and force vector b

        Returns
        -------
        amat : array_like, dtype='double'
            shape=(3 * num_atoms * ndata_training, num_fc_irred)
            The sensing matrix A calculated from the displacements.

        bvec : array_like, dtype='double'
            shape=(3 * num_atoms * ndata_training,)
            The vector b calculated from the atomic forces.


        Note
        ----
        From the amat (``A``) and bvec (``b``), the force constant vector ``x``
        can be obtained by solving the least-square problem:
        x = argmin_{x} | Ax-b|^{2}.
        """

        if self._id is None:
            self._show_error_message()

        norder = self._norder
        nat = len(self._xcoord)
        ndata_used = self._get_ndata_used()

        fc_length = 0
        for i in range(norder):
            fc_length += self._get_number_of_irred_fc_elements(i + 1)

        amat = np.zeros(3 * nat * ndata_used * fc_length, dtype='double')
        bvec = np.zeros(3 * nat * ndata_used)
        alm.get_matrix_elements(self._id, nat, ndata_used, amat, bvec)

        return np.reshape(amat, (3 * nat * ndata_used, fc_length), order='F'), bvec

    def _set_cell(self):
        if self._id is None:
            self._show_error_message()

        alm.set_cell(self._id, self._lavec, self._xcoord, self._kd)

    def _set_norder(self):
        if self._id is None:
            self._show_error_message()

        alm.set_norder(self._id, self._norder)

    def _get_ndata_used(self):
        if self._id is None:
            self._show_error_message()

        ndata_used = alm.get_ndata_used(self._id)
        return ndata_used

    def _get_id(self):
        return self._id

    def _get_number_of_displacement_patterns(self, fc_order):
        return alm.get_number_of_displacement_patterns(self._id, fc_order)

    def _get_numbers_of_displacements(self, fc_order):
        num_disp_patterns = self._get_number_of_displacement_patterns(
            fc_order)
        numbers = np.zeros(num_disp_patterns, dtype='intc')
        alm.get_numbers_of_displacements(self._id, numbers, fc_order)
        return numbers

    def _get_number_of_fc_elements(self, fc_order):  # harmonic: fc_order=1
        return alm.get_number_of_fc_elements(self._id, fc_order)

    def _get_number_of_irred_fc_elements(self, fc_order):  # harmonic: fc_order=1
        return alm.get_number_of_irred_fc_elements(self._id, fc_order)

    def _show_error_message(self):
        print("This ALM object has to be initialized by ALM::alm_new()")
        raise RuntimeError
