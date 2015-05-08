#!/usr/bin/env python

# packages
from __future__ import print_function
import numpy as np
import pkes

###############################################################################

class DecaySolver(object):

    ###########################################################################
    ## Constructor
    ##
    def __init__(self):

        # initialize private class members
        self._end_times = None
        self._power = None
        self._num_time_steps = None
        self._decay_power = {} 

    ###########################################################################
    ## Properties
    ##

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

    def get_num_time_steps(self):
        return self._num_time_steps

    def set_num_time_steps(self, num_time_steps):

        # check num_time_steps
        if not isinstance(num_time_steps, (list, np.ndarray)):
            raise TypeError("Number of time steps not an array.")

        # set number of time steps
        self._num_time_steps = num_time_steps

    num_time_steps = property(get_num_time_steps, set_num_time_steps)

        #######################################################################

    def get_decay_power(self):
        return self._decay_power

    decay_power = property(get_decay_power)

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

        # open output file
        fh = open("decay_power.dat", "w")

        # set up initial time step
        time = 0.0
        t_idx = 0
        t_cmp = self.num_time_steps[t_idx]
        dt = self.end_times[t_idx] / float(self.num_time_steps[t_idx])

        # perform integration
        for i in range(np.sum(self.num_time_steps)):

            # check coarse time index
            if i == t_cmp:
                t_idx += 1
                t_cmp += self.num_time_steps[t_idx]
                dt = self.end_times[t_idx] / float(self.num_time_steps[t_idx])

            # calculate time
            time += dt

            # set up time information
            time_decay = 0.0
            t_idx_decay = 0
            t_cmp_decay = self.num_time_steps[t_idx_decay]
            dt_decay = self.end_times[t_idx_decay] / float(
                       self.num_time_steps[t_idx_decay])

            # set up last power
            power_last = self.power.data[0]

            # loop around all previous timesteps
            for j in range(i + 1):

                # check decay time index
                if j == t_cmp_decay:
                    t_idx_decay += 1
                    t_cmp_decay += self.num_time_steps[t_idx_decay]
                    dt_decay = self.end_times[t_idx_decay] / float(
                       self.num_time_steps[t_idx_decay])

                # calculate time for decay
                time_decay += dt_decay

                # calculate length of time between times
                t = time - time_decay

                # calculate average power
                power = self.power.interpolate(time_decay)
                power_avg = (power + power_last) / 2.0
            
                # loop around nuclides
                for key, val in pkes.decay_data.iteritems():

                    # extract nuclide data
                    f = val["power_fraction"]
                    Q = val["energy_per_fission"]
                    alpha = val["alpha"]
                    lamb = val["lambda"]

                    # calculate decay power of this nuclide
                    self._decay_power[key][i] += power_avg*f/Q * \
                         np.sum(alpha/lamb*(1.0 - np.exp(-lamb*dt_decay)) * \
                         np.exp(-lamb*t))
                         
                # move current power to last power
                power_last = power

            # print to screen and write to file
            fh.write("{0} {1} {2} {3}\n".format(time, power_avg,
              self._decay_power["U-235"][i], self._decay_power["U-238"][i]))
            print("{0} {1} {2} {3}".format(time, power_avg,
              self._decay_power["U-235"][i], self._decay_power["U-238"][i]))

        fh.close()

        #######################################################################

    def _validate(self):

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

        # set up dictionary of decay powers
        for key, val in pkes.decay_data.iteritems():
            self._decay_power[key] = np.zeros((np.sum(self.num_time_steps)
                                               + 1))
    ##
