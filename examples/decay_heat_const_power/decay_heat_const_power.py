#!/usr/bin/env python

# packages
from __future__ import print_function
import pkes

# set up power trace
power = pkes.Solution(2)
power.add_data_point(0, 0.0, 3400.e6)
power.add_data_point(1, 10000.0, 3400.e6)
print(power)

# set up solver
solver = pkes.DecaySolver()
solver.end_times = [10000.0]
solver.num_time_steps = [1000]
solver.power = power

# solve problem
solver.solve()