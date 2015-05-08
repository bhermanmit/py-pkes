#!/usr/bin/env python

# packages
from __future__ import print_function
import numpy as np
import pkes

###############################################################################

class DecaySolverNew(object):

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

        # create temporary array for decay power
        decay_temp = {}
        for key in self.decay_power.iterkeys():
            decay_temp[key] = np.zeros((23, np.sum(self.num_time_steps) + 1))

        # open output file
        fh = open("decay_power_new.dat", "w")

        # set up initial time step
        time = 0.0
        t_idx = 0
        t_cmp = self.num_time_steps[t_idx]
        dt = self.end_times[t_idx] / float(self.num_time_steps[t_idx])

        # set up last power
        power_last = self.power.data[0]

        # perform integration
        for i in range(np.sum(self.num_time_steps)):

            # check coarse time index
            if i == t_cmp:
                t_idx += 1
                t_cmp += self.num_time_steps[t_idx]
                dt = (self.end_times[t_idx] - self.end_times[t_idx-1]) / \
                    float(self.num_time_steps[t_idx])

            # calculate time
            time += dt

            # calculate average power
            power = self.power.interpolate(time)
            power_avg = (power + power_last) / 2.0
            
            # decay to this time
            for key, val in pkes.decay_data.iteritems():

                # extract decay constants
                lamb = val["lambda"]

                # calculate decay multiplier
                exp = np.exp(-lamb*dt)

                # multiply through all timesteps
                decay_temp[key] = (decay_temp[key].T*exp).T

            # add in new component
            for key, val in pkes.decay_data.iteritems():

                # extract nuclide data
                f = val["power_fraction"]
                Q = val["energy_per_fission"]
                alpha = val["alpha"]
                lamb = val["lambda"]

                # calculate decay power of this nuclide for the current time
                decay_temp[key][:, i] += alpha/lamb*(1.0 - np.exp(-lamb*dt))

                # calculate actual decay power
                self._decay_power[key][i] += power_avg*f/Q * \
                     np.sum(np.sum(decay_temp[key]))

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
