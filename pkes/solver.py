#!/usr/bin/env python

# packages
from __future__ import print_function
from scipy.integrate import ode
import numpy as np
import pkes

###############################################################################

class PKESolver(object):

    ###########################################################################
    ## Constructor
    ##
    def __init__(self):

        # initialize private class members
        self._material = None
        self._reactivity = None
        self._end_time = None
        self._time_step = None
        self._ode = None
        self._power = None
        self._num_time_steps = None
        self._initial_power = 1.0
        self._matrix = None
        self._b = None

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

    def get_reactivity(self):
        return self._reactivity

    def set_reactivity(self, reactivity):

        # check reactivity types
        if not isinstance(reactivity, pkes.Solution):
            raise TypeError("Reactivity is not a PKE solution.")

        # set reactivity
        self._reactivity = reactivity

    reactivity = property(get_reactivity, set_reactivity)

        #######################################################################

    def get_end_time(self):
        return self._end_time

    def set_end_time(self, end_time):

        # check end time type
        if not isinstance(end_time, float):
            raise TypeError("End time is not a float.")

        # set end time
        self._end_time = end_time

    end_time = property(get_end_time, set_end_time)

        #######################################################################

    def get_time_step(self):
        return self._time_step

    time_step = property(get_time_step)

        #######################################################################

    def get_power(self):
        return self._power

    power = property(get_power)

        #######################################################################

    def get_num_time_steps(self):
        return self._num_time_steps

    def set_num_time_steps(self, num_time_steps):

        # check num_time_steps is an integer
        if not isinstance(num_time_steps, int):
            raise TypeError("Number of time steps not an integer.")

        # set number of time steps
        self._num_time_steps = num_time_steps

    num_time_steps = property(get_num_time_steps, set_num_time_steps)

        #######################################################################

    def get_initial_power(self):
        return self._initial_power

    def set_initial_power(self, initial_power):

        # check initial power type
        if not isinstance(initial_power, int):
            raise TypeError("Initial power is not a float.")

        # set initial power
        self._initial_power = initial_power

    initial_power = property(get_initial_power, set_initial_power)

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

        # calculate steady state
        x = self._steady_state()

        # record initial power in vector
        self.power.add_data_point(0, 0.0, x[0])

        # open output file
        fh = open("output.dat", "w")

        # perform integration
        time = 0.0
        for i in range(self.num_time_steps):

            # calculate time step
            dt = self.end_time / float(self.num_time_steps)

            # calculate time
            time += dt

            # get reactivity at this time
            rho = self.reactivity.interpolate(time)

            # create coefficient matrix
            self._create_matrix(dt, rho)

            # set up rhs
            self._rhs(dt, x)

            # solve matrix
            x = np.linalg.solve(self._matrix, self._b)

            # extract power
            self.power.add_data_point(i+1, time, x[0])

            # print to screen and write to file
            fh.write("{0} {1}\n".format(time, x[0]))
            print("{0} {1}".format(time, x[0]))

        fh.close()

        #######################################################################

    def _validate(self):

        # mats
        if self.material is None:
            raise ValueError("Material not set in PKESolver.")

        # reactivity
        if self.reactivity is None:
            raise ValueError("Reactivity not set in PKESolver.")

        # end time
        if self.end_time is None:
            raise ValueError("End time not set in PKESolver.")

        # number of time steps
        if self.num_time_steps is None:
            raise ValueError("Number of time steps not set in PKESolver.")

        #######################################################################

    def _allocate(self):

        # calculate time step
        self._time_step = self.end_time / float(self.num_time_steps)

        # allocate power vector
        self._power = pkes.Solution(self.num_time_steps + 1)

        # coefficient matrix
        self._matrix = np.zeros((self.material.num_precs+1,
                                 self.material.num_precs+1))

        # right hand side vector
        self._b = np.zeros((self.material.num_precs+1))

        #######################################################################

    def _steady_state(self):

        # allocate numpy vector
        x = np.zeros((self.material.num_precs + 1))

        # fill in the initial power
        x[0] = self.initial_power

        # loop around precursors
        for i in range(self.material.num_precs):

            # calculate steady state precs
            x[i+1] = self.material.beta[i] * self.initial_power / \
                     (self.material.decay[i] * self.material.pnl)

        return x

        #######################################################################

    def _create_matrix(self, dt, rho):

        # power
        self._matrix[0,0] = 1.0/dt - (rho - self.material.beta_total) / \
                            self.material.pnl

        # precursors
        for i in range(self.material.num_precs):

            # power row
            self._matrix[0,i+1] = -self.material.decay[i]

            # power column
            self._matrix[i+1,0] = -self.material.beta[i] / self.material.pnl

            # diagonal
            self._matrix[i+1,i+1] = 1.0/dt + self.material.decay[i]

        #######################################################################

    def _rhs(self, dt, x):

        # power
        self._b[0] = x[0] / dt

        # precursors
        for i in range(self.material.num_precs):
            self._b[i+1] = x[i+1] / dt

    ##
