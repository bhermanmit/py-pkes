"""This module calculations the actinide power after
an irradiation.
"""

import numpy as np

# Data taken from Table 3.6

NAMES = ["Am-241", "Pu-241", "Pu-240", "Pu-239", "Pu-238", "Cm-244", "Cm-242"]

DECAY = np.array([5.078e-11, 1.531e-9, 3.3472e-12, 9.111e-13, 2.504e-10,
        1.213e-9, 4.923e-8])

# Take at 20 MWd/kgU at 3 wt%
BETA = np.array([2.868e-3, 2.147e-3, 8.586e-3, 9.721e-3, 2.778e-2,
                 4.563e-3, 5.457e-1])

BETA_p = BETA
BETA_p[0] = BETA[0] - BETA[1]*(5.629/5.361e-3)*(DECAY[0]/(DECAY[0] - DECAY[1]))
BETA_p[1] = BETA[1]*(1.0 + (5.629/5.361e-3)*(DECAY[0]/(DECAY[0] - DECAY[1])))

# Initial mass of U taken from App. K in T&K
KG_U = 89000.

def actinide_solve(end_time, num_steps):

    time = np.zeros(num_steps + 1)

    dt = end_time / float(num_steps)

    power = np.zeros((len(NAMES), num_steps + 1))

    for i in range(num_steps + 1):

        time[i] = i*dt

        power[:, i] += BETA[:]*np.exp(-DECAY[:]*time[i]) 

    # scale specific power to actual power
    power *= KG_U

    with open("actinide_decay.txt", "w") as fh:
        for i in range(num_steps + 1):
            fh.write("{0} {1}\n".format(time[i], " ".join(map(str, power[:, i]))))
