"""This module calculations the actinide power after
an irradiation.
"""

import numpy as np

# Data taken from Table 3.6

NAMES = ["Am-241", "Pu-241", "Pu-240", "Pu-239", "Pu-238", "Cm-244", "Cm-242"]

DECAY = np.array([5.078e-11, 1.531e-9, 3.3472e-12, 9.111e-13, 2.504e-10,
        1.213e-9, 4.923e-8])

# Take at 50 MWd/kgU at 4 wt%
BETA = np.array([1.344e-2, 5.626e-3, 2.083e-2, 1.159e-2, 2.607e-1,
                 3.108e-1, 4.618e0])

BETA_p = BETA
BETA_p[0] = BETA[0] - BETA[1]*(5.629/5.361e-3)*(DECAY[0]/(DECAY[0] - DECAY[1]))
BETA_p[1] = BETA[1]*(1.0 + (5.629/5.361e-3)*(DECAY[0]/(DECAY[0] - DECAY[1])))

# Initial mass of U taken from App. K in T&K
KG_U = 24000.

def actinide_solve(time):

    print("BETA IS:", BETA)
    print("BETA_prime IS:", BETA_p)

    num_steps = len(time)

    power = np.zeros((num_steps,))

    for i in range(num_steps):

        power[i] = np.sum(BETA_p*np.exp(-DECAY*time[i])) 
        print(time[i], power[i])

    # scale specific power to actual power then convert to MW
    power *= KG_U *1.e-6

    with open("actinide_decay.dat", "w") as fh:
        for i in range(num_steps):
            fh.write("{0} {1}\n".format(time[i], power[i])) 
