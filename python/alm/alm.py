import numpy as np
from . import _alm as alm

class ALM:
    def __init__(self, lavec, xcoord, kd):
        self._id = None
        self._lavec = np.array(lavec, dtype='double', order='C')
        self._xcoord = np.array(xcoord, dtype='double', order='C')
        self._kd = np.array(kd, dtype='intc', order='C')
        self._iconst = 11

    def __enter__(self):
        self.alm_new()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.alm_delete()

    def alm_new(self):
        if self._id is None:
            self._id = alm.alm_new()
            if self._id < 0:
                print("Too many ALM objects")
                raise
            self._set_cell()
        else:
            print("This ALM object is already initialized.")
            raise
    
    def alm_delete(self):
        if self._id is None:
            self._show_error_message()

        alm.alm_delete(self._id)
        self._id = None


    def run_suggest(self):
        if self._id is None:
            self._show_error_message()

        alm.run_suggest(self._id)
    
    def optimize(self):
        if self._id is None:
            self._show_error_message()

        info = alm.optimize(self._id)
        return info
    
    def set_displacement_and_force(self, u, f):
        if self._id is None:
            self._show_error_message()

        alm.set_displacement_and_force(
            self._id,
            np.array(u, dtype='double', order='C'),
            np.array(f, dtype='double', order='C'))

    def find_force_constant(self, norder, rcs, nbody = []):
        # TODO: support nbody option
        if self._id is None:
            self._show_error_message()

        self._norder = norder
        self._set_norder()
        alm.set_cutoff_radii(self._id,
                             np.array(rcs, dtype='double', order='C'))

        alm.generate_force_constant(self._id)


    def set_constraint(self, translation=True, rotation=False):
        if rotation is True:
            print("Rotational invariance is not supported in API.")
            raise

        if translation is True:
            iconst = 11
        else:
            iconst = 0

        self._iconst = iconst
        alm.set_constraint_type(self._id, self._iconst)


    def set_cutoff_radii(self, rcs):
        if self._id is None:
            self._show_error_message()

        alm.set_cutoff_radii(self._id,
                             np.array(rcs, dtype='double', order='C'))


    def get_atom_mapping_by_pure_translations(self):
        if self._id is None:
            self._show_error_message()

        map_p2s = np.zeros(len(self._xcoord), dtype='intc')
        ntrans = alm.get_atom_mapping_by_pure_translations(self._id, map_p2s)
        return map_p2s.reshape((ntrans, -1))
    

    def get_displacement_patterns(self, fc_order):
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
    

    def get_fc(self, fc_order, mode='origin'): # harmonic: fc_order=1
        if self._id is None:
            self._show_error_message()

        if mode == 'origin':
        
            fc_length = self._get_number_of_fc_elements(fc_order)
            fc_values = np.zeros(fc_length, dtype='double')
            elem_indices = np.zeros((fc_length, fc_order + 1),
                                    dtype='intc', order='C')

            alm.get_fc_origin(self._id, fc_values, elem_indices)

            return fc_values, elem_indices

        elif mode == 'irreducible' or mode == 'irred':

            fc_length = self._get_number_of_irred_fc_elements(fc_order)
            fc_values = np.zeros(fc_length, dtype='double')
            elem_indices = np.zeros((fc_length, fc_order + 1),
                                    dtype='intc', order='C')

            alm.get_fc_irreducible(self._id, fc_values, elem_indices)

            return fc_values, elem_indices

        elif mode == 'all':
            
            map_p2s = np.zeros(len(self._xcoord), dtype='intc')
            ntrans = alm.get_atom_mapping_by_pure_translations(self._id, map_p2s)
            fc_length = self._get_number_of_fc_elements(fc_order) * ntrans
            fc_values = np.zeros(fc_length, dtype='double')
            elem_indices = np.zeros((fc_length, fc_order + 1),
                                    dtype='intc', order='C')

            alm.get_fc_all(self._id, fc_values, elem_indices)

            return fc_values, elem_indices

        else:
            print("Invalid mode in get_fc.")
            exit(1)
        

    def set_fc(self, fc_in):
        if self._id is None:
            self._show_error_message()

        norder = self._norder
        fc_length_irred = 0
        for i in range(norder):
            fc_length_irred += self._get_number_of_irred_fc_elements(i + 1)
        
        if fc_length_irred != len(fc_in):
            print("The size of the given force constant array is incorrect.")
            exit(1)

        alm.set_fc(self._id, np.array(fc_in, dtype='double', order='C'))


    def get_matrix_elements(self):
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
        
        return np.reshape(amat, (3 * nat* ndata_used, fc_length), order='F'), bvec

    
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
    
    def get_id(self):
        return self._id

    def _get_number_of_displacement_patterns(self, fc_order):
        return alm.get_number_of_displacement_patterns(self._id, fc_order)
    
    def _get_numbers_of_displacements(self, fc_order):
        num_disp_patterns = self._get_number_of_displacement_patterns(
            fc_order)
        numbers = np.zeros(num_disp_patterns, dtype='intc')
        alm.get_numbers_of_displacements(self._id, numbers, fc_order)
        return numbers
    
    def _get_number_of_fc_elements(self, fc_order): # harmonic: fc_order=1
        return alm.get_number_of_fc_elements(self._id, fc_order)

    def _get_number_of_irred_fc_elements(self, fc_order): # harmonic: fc_order=1
        return alm.get_number_of_irred_fc_elements(self._id, fc_order)
    
    def _show_error_message(self):
        print("This ALM object has to be initialized by ALM::alm_new()")
        raise
