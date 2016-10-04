import numpy as np

POWER = 3000.

def analytic(time):

    if time <= 400.0:
        res = -6.14575e-3*np.log(time) + 0.060157
    elif 400 < time <= 4.e5:
        res = 1.40680e-1*time**-0.286
    elif 4.e5 < time <= 4.e6:
        res = 8.70300e-1*time**-0.4255
    elif 4.e6 < time <= 4.e7:
        res = 1.28420e1*time**-0.6014
    elif 4.e7 < time <= 4.e8:
        res = 4.03830e4*time**-1.0675
    else:
        res = 3.91130e-5*np.exp(-7.3541e-10*time)

    return res

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
       handle.write("{} {} {}\n".format(time, power, POWER*analytic(time)))
