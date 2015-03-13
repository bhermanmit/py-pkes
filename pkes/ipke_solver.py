#!/usr/bin/env python

# packages
from __future__ import print_function
import numpy as np
import pkes

###############################################################################

class IPKESolver(object):

    ###########################################################################
    ## Constructor
    ##
    def __init__(self):

        # initialize private class members
        self._material = None
        self._reactivity = None
        self._end_times = None
        self._power = None
        self._num_time_steps = None

    ###########################################################################
    ## Properties
    ##

    def get_material(self):
        return self._material

    def set_material(self, material):

        # check material type
        if not isinstance(material, pkes.Material):
            raise TypeError("Solver material is not a PKE material.")

        # set mats
        self._material = material

    material = property(get_material, set_material)

        #######################################################################

    def get_power(self):
        return self._power

    def set_power(self, power):

        # check power types
        if not isinstance(power, pkes.Solution):
            raise TypeError("Power is not a PKE solution.")

        # set power
        self._power = power

    power = property(get_power, set_power)

        #######################################################################

    def get_end_times(self):
        return self._end_times

    def set_end_times(self, end_times):

        # check end times type
        if not isinstance(end_times, (list, np.ndarray)):
            raise TypeError("End times is not an array.")

        # set end times
        self._end_times = end_times

    end_times = property(get_end_times, set_end_times)

        #######################################################################

    def get_reactivity(self):
        return self._reactivity

    reactivity = property(get_reactivity)

        #######################################################################

    def get_num_time_steps(self):
        return self._num_time_steps

    def set_num_time_steps(self, num_time_steps):

        # check num_time_steps
        if not isinstance(num_time_steps, (list, np.ndarray)):
            raise TypeError("Number of time steps not an array.")

        # set number of time steps
        self._num_time_steps = num_time_steps

    num_time_steps = property(get_num_time_steps, set_num_time_steps)

    ##

    ###########################################################################
    ## Methods
    ##
    ##
    def solve(self):

        # check object
        self._validate()

        # allocate object
        self._allocate()

        # calculate initial precursors 
        C = self._initial_precursors(self.power.data[0])

        # open output file
        fh = open("reactivity.dat", "w")

        # set up initial time step
        time = 0.0
        t_idx = 0
        t_cmp = self.num_time_steps[t_idx]
        dt = self.end_times[t_idx] / float(self.num_time_steps[t_idx])

        # set up last power
        power_last = self.power.data[0]

        # set initial reactivity as zero
        self.reactivity.add_data_point(0, 0.0, 0.0)

        # perform integration
        for i in range(np.sum(self.num_time_steps)):

            # check coarse time index
            if i == t_cmp:
                t_idx += 1
                t_cmp = self.num_time_steps[t_idx]
                dt = self.end_times[t_idx] / float(self.num_time_steps[t_idx])

            # calculate time
            time += dt

            # calculate average power
            power = self.power.interpolate(time)
            power_avg = (power + power_last) / 2.0

            # calculate new precursors
            C = C*np.exp(-self.material.decay*dt) + self.material.beta / \
                (self.material.decay*self.material.pnl) * power_avg * \
                (1.0 - np.exp(-self.material.decay*dt))

            # calculate new power
            rho = self.material.pnl/power*(power - power_last)/dt + \
                  self.material.beta_total - self.material.pnl/power * \
                  np.sum(self.material.decay*C)

            # set new reacitivty 
            self.reactivity.add_data_point(i+1, time, rho)

            # print to screen and write to file
            fh.write("{0} {1} {2}\n".format(time, rho, power))
            print("{0} {1} {2}".format(time, rho, power))

            # move current power to last power
            power_last = power

        fh.close()

        #######################################################################

    def _validate(self):

        # mats
        if self.material is None:
            raise ValueError("Material not set in PKESolver.")

        # power
        if self.power is None:
            raise ValueError("Power not set in PKESolver.")

        # end time
        if self.end_times is None:
            raise ValueError("End times not set in PKESolver.")

        # number of time steps
        if self.num_time_steps is None:
            raise ValueError("Number of time steps not set in PKESolver.")

        #######################################################################

    def _allocate(self):

        # allocate reactivity vector
        self._reactivity = pkes.Solution(np.sum(self.num_time_steps) + 1)

        #######################################################################

    def _initial_precursors(self, power):

        # calculate steady state precs
        C = self.material.beta * power / \
            (self.material.decay * self.material.pnl)

        return C
    ##
