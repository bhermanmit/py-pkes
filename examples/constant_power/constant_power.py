#!/usr/bin/env python

# packages
from __future__ import print_function
import pkes
import numpy as np

# set material information
mat = pkes.Material()
mat.beta = np.array([0.000218, 0.001023, 0.000605, 0.00131, 0.00220, 
            0.000060, 0.000540, 0.000152])
halflife = np.array([55.6, 24.5, 16.3, 5.21, 2.37, 1.04, 0.424, 0.195])
mat.decay = np.log(2.0)/halflife
mat.pnl = 0.0001866
print(mat)

# set up power trace
power = pkes.Solution(2)
power.add_data_point(0, 0.0, 5.0)
power.add_data_point(1, 1.0, 5.0)
print(power)

# set up solver
solver = pkes.IPKESolver()
solver.material = mat
solver.end_times = [1.0]
solver.num_time_steps = [100]
solver.power = power

# solve problem
solver.solve()
