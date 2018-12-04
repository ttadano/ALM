#!/usr/bin/env python
# coding: utf-8

import numpy as np
from alm import ALM


lavec = np.loadtxt('lavec.dat')
xcoord = np.loadtxt('xcoord.dat')
kd = np.loadtxt('kd.dat')
force = np.loadtxt("sic_force.dat").reshape((-1, 64, 3))
disp = np.loadtxt("sic_disp.dat").reshape((-1, 64, 3))


# alm.alm_new() and alm.alm_delete() are done by 'with' statement
with ALM(lavec, xcoord, kd) as alm:
    alm.define(2)
    alm.set_displacement_and_force(disp, force+5e-4*np.random.random(force.shape))
    alm.set_verbosity(1)
    info = alm.optimize(solver='dense')

    fc_values, elem_indices = alm.get_fc(1) # harmonic fc
    print('RMS FC2: %.4f' % np.sqrt((fc_values**2).mean()))
    fc_values, elem_indices = alm.get_fc(2) # cubic fc
    print('RMS FC3: %.4f' % np.sqrt((fc_values**2).mean()))






