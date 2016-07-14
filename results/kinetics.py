#!/usr/bin/env python

# packages
from __future__ import print_function
import pkes
import numpy as np
import h5py

TO_SECONDS=86400.0

# set material information
mat = pkes.Material()
mat.beta = np.array([0.000218, 0.001023, 0.000605, 0.00131, 0.00220, 
            0.000060, 0.000540, 0.000152])
halflife = np.array([55.6, 24.5, 16.3, 5.21, 2.37, 1.04, 0.424, 0.195])
mat.decay = np.log(2.0)/halflife
mat.pnl = 0.0001866
print(mat)

# get reactivity
fh = h5py.File("reactivity.h5", "r")
time = fh["time"].value
rho = fh["rho"].value
n_points = len(time) 
reactivity = pkes.Solution(n_points + 1)
for i in range(n_points):
    reactivity.add_data_point(i, time[i], rho[i])
reactivity.add_data_point(n_points, 72.0*TO_SECONDS + 5.0, -0.0050*np.sum(mat.beta))

# set up point kinetics solver
pke_solver = pkes.PKEODESolver()
pke_solver.material = mat
pke_solver.end_time = 72.0*TO_SECONDS
pke_solver.reactivity = reactivity
pke_solver.initial_power = 1.e-6
pke_solver.max_step = 1.0*TO_SECONDS

# solve point kinetics
pke_solver.solve()
