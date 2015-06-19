#!/usr/bin/env python

# packages
from __future__ import print_function
import pkes
import numpy as np

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
with open("reactivity.dat") as handle:
    lines = handle.read().splitlines()
n_lines = len(lines)
reactivity = pkes.Solution(n_lines)
for i in range(len(lines)):
    time = float(lines[i].split()[0])
    reactivity_pt = float(lines[i].split()[1])
    reactivity.add_data_point(i, time, reactivity_pt)

# set up point kinetics solver
pke_solver = pkes.PKESolver()
pke_solver.material = mat
pke_solver.end_times = np.asarray([1.0, 2.0, 23.0, 25.0, 72.0, 100.0, 1000.0])
pke_solver.num_time_steps = [10, 10000, 100000, 1000, 1000, 1000, 1000]
pke_solver.reactivity = reactivity
pke_solver.initial_power = 30.0

# solve point kinetics
pke_solver.solve()
