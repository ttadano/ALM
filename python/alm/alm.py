import numpy as np
from . import _alm as alm

class ALM:
    def __init__(self):
        self._id = None

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
    
    def run_fitting(self):
        if self._id is None:
            self._show_error_message()

        alm.run_fitting(self._id)
    
    def set_cell(self, lavec, xcoord, kd):
        if self._id is None:
            self._show_error_message()

        alm.set_cell(self._id,
                     np.array(lavec, dtype='double', order='C'),
                     np.array(xcoord, dtype='double', order='C'),
                     np.array(kd, dtype='intc', order='C'))
    
    def set_displacement_and_force(self, u, f):
        if self._id is None:
            self._show_error_message()

        alm.set_displacement_and_force(
            self._id,
            np.array(u, dtype='double', order='C'),
            np.array(f, dtype='double', order='C'))
    
    def set_norder(self, norder):
        if self._id is None:
            self._show_error_message()

        alm.set_norder(self._id, norder)
    
    def set_cutoff_radii(self, rcs):
        if self._id is None:
            self._show_error_message()

        alm.set_cutoff_radii(self._id,
                             np.array(rcs, dtype='double', order='C'))
    
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
    
    def get_fc(self, fc_order): # harmonic: fc_order=1
        if self._id is None:
            self._show_error_message()

        fc_length = self._get_number_of_fc_elements(fc_order)
        fc_values = np.zeros(fc_length, dtype='double')
        elem_indices = np.zeros((fc_length, fc_order + 1),
                                dtype='intc', order='C')
        alm.get_fc(self._id, fc_values, elem_indices)
        return fc_values, elem_indices
    
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
    
    def _show_error_message(self):
        print("This ALM object has to be initialized by ALM::alm_new()")
        raise
