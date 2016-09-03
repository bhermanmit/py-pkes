#!/usr/bin/env python

# packages
from __future__ import print_function
import pkes
import numpy as np
import pkes.decay_solver_new

TO_SECONDS=86400.0

# get power
with open("power.dat") as handle:
    lines = handle.read().splitlines()
n_lines = len(lines)
power = pkes.Solution(n_lines)
for i in range(len(lines)):
    time = float(lines[i].split()[1])
    power_pt = float(lines[i].split()[3])
    power.add_data_point(i, time, power_pt)

# set up solver
solver = pkes.decay_solver_new.DecaySolverNew()
solver.end_times = np.asarray([1.0, 2.0, 23.0, 25.0, 72.0, 100.0, 720.0, 1000.0]) * TO_SECONDS
solver.num_time_steps = [1000, 1000, 3000, 3000, 10000, 10000, 10000, 10000]
solver.power = power

# solve problem
solver.solve()
