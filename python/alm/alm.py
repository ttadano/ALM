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
    alm.set_displacement_and_force(u, f)

def set_norder(norder):
    alm.set_norder(norder)

def set_cutoff_radii(rcs):
    alm.set_cutoff_radii(np.array(rcs, dtype='double', order='C'))
