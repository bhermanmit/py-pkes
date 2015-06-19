#!/usr/bin/env python

# packages
from __future__ import print_function
import pkes
import numpy as np

# get power
with open("power.dat") as handle:
    lines = handle.read().splitlines()
n_lines = len(lines)
power = pkes.Solution(n_lines)
for i in range(len(lines)):
    time = float(lines[i].split()[0])
    power_pt = float(lines[i].split()[2])
    power.add_data_point(i, time, power_pt)

# set up solver
solver = pkes.DecaySolver()
solver.end_times = np.asarray([1.0, 2.0, 23.0, 25.0, 72.0, 100.0, 1000.0])
solver.num_time_steps = [10, 10000, 100000, 1000, 1000, 1000, 1000]
solver.power = power

# solve problem
solver.solve()
