#!/usr/bin/env python

# packages
from __future__ import print_function
import pkes
import pkes.decay_solver_new

# set up power trace
power = pkes.Solution(2)
power.add_data_point(0, 0.0, 3400.)
power.add_data_point(1, 10000.0, 3400.)
power.add_data_point(2, 10000.1, 0.0)
print(power)

# set up solver
solver = pkes.decay_solver_new.DecaySolverNew()
# solver = pkes.DecaySolver()
solver.end_times = [20000.0]
solver.num_time_steps = [1000]
solver.power = power

# solve problem
solver.solve()
