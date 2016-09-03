import pkes

import numpy as np

with open('decay_power_new.dat') as handle:
    lines = handle.read().splitlines()

decay_time = np.empty((len(lines),), dtype=float)
power_interp = np.empty((len(lines),), dtype=float)

for i, aline in enumerate(lines):
    sline = aline.split()
    decay_time[i] = float(sline[0])

pkes.actinide_solve(decay_time)
