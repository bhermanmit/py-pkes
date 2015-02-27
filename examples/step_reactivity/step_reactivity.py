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

# set up reactivity trace
reactivity = pkes.Solution(4)
reactivity.add_data_point(0, 0.0, 0.0)
reactivity.add_data_point(1, 0.5, 0.0)
reactivity.add_data_point(2, 0.51, 0.002)
reactivity.add_data_point(3, 100.0, 0.002)
print(reactivity)

# set up solver
solver = pkes.PKESolver()
solver.material = mat
solver.end_times = [10.0, 90.0]
solver.num_time_steps = [1000, 900]
solver.reactivity = reactivity

# solve problem
solver.solve()
