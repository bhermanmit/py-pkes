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

# set up power trace
power = pkes.Solution(4)
power.add_data_point(0, 0.0,  30.0)
power.add_data_point(1, 1.0*TO_SECONDS,  30.0)
power.add_data_point(2, 24.0*TO_SECONDS,  3000.0)
power.add_data_point(3, 72.0*TO_SECONDS,  3000.0)
print(power)

# set up inverse kinetics solver
ipke_solver = pkes.IPKESolver()
ipke_solver.material = mat
ipke_solver.end_times = np.asarray([1.0, 2.0, 23.0, 25.0, 72.0])*TO_SECONDS
ipke_solver.num_time_steps = [100000, 100000, 3000000, 3000000, 1000]
ipke_solver.power_input = power
ipke_solver.solve()
