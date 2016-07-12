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
power = pkes.Solution(6)
power.add_data_point(0, 0.0,  1.0)
power.add_data_point(1, 2.0,  1.0)
power.add_data_point(2, 2.01, 1000.0)
power.add_data_point(3, 10.0, 1000.0)
power.add_data_point(4, 12.0, 100.0)
power.add_data_point(5, 20.0, 100.0)
print(power)

# set up inverse kinetics solver
ipke_solver = pkes.IPKESolver()
ipke_solver.material = mat
ipke_solver.end_times = [20.0]
ipke_solver.num_time_steps = [20000]
ipke_solver.power_input = power

# solve inverse kinetics
ipke_solver.solve()

# set up point kinetics solver
pke_solver = pkes.PKEODESolver()
pke_solver.material = mat
pke_solver.end_time = 20.0
pke_solver.max_step = 1.0
pke_solver.reactivity = ipke_solver.reactivity

# solve point kinetics
pke_solver.solve()