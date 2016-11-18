from . import _alm as alm
import numpy as np

def alm_new():
    alm.alm_new()

def alm_delete():
    alm.alm_delete()

def run_suggest():
    alm.run_suggest()

def run_fitting():
    alm.run_fitting()

def set_cell(lavec, xcoord, kd):
    alm.set_cell(np.array(lavec, dtype='double', order='C'),
                 np.array(xcoord, dtype='double', order='C'),
                 np.array(kd, dtype='intc', order='C'))

def set_displacement_and_force(u, f):
    alm.set_displacement_and_force(np.array(u, dtype='double', order='C'),
                                   np.array(f, dtype='double', order='C'))

def set_norder(norder):
    alm.set_norder(norder)

def set_cutoff_radii(rcs):
    alm.set_cutoff_radii(np.array(rcs, dtype='double', order='C'))

def get_fc(fc_order): # harmonic: fc_order=1
    fc_length = _get_fc_length(fc_order)
    fc_values = np.zeros(fc_length, dtype='double')
    elem_indices = np.zeros((fc_length, fc_order + 1), dtype='intc', order='C')
    alm.get_fc(fc_values, elem_indices)
    return fc_values, elem_indices

def _get_fc_length(fc_order): # harmonic: fc_order=1
    return alm.get_fc_length(fc_order)

