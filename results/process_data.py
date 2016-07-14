import numpy as np

with open('power.dat') as handle:
    lines = handle.read().splitlines()

power_time = np.empty((len(lines),), dtype=float)
power_power = np.empty((len(lines),), dtype=float)

for i, aline in enumerate(lines):
    sline = aline.split()
    power_time[i] = float(sline[1])
    power_power[i] = float(sline[3])

with open('decay_power_new.dat') as handle:
    lines = handle.read().splitlines()

decay_time = np.empty((len(lines),), dtype=float)
power_interp = np.empty((len(lines),), dtype=float)

for i, aline in enumerate(lines):
    sline = aline.split()
    decay_time[i] = float(sline[0])
    power_interp[i] = np.interp(decay_time[i], power_time, power_power)

with open('power_interp.dat', 'w') as handle:
   for time, power in zip(decay_time, power_interp):
       handle.write("{} {}\n".format(time, power))
