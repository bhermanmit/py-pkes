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
        y0 = self._steady_state()

        # set up ode integrator
        self._ode = ode(_f, _jac).\
                    set_integrator('vode', method='bdf', with_jacobian=True)
        self._ode.set_initial_value(y0)
        self._ode.set_f_params(self.material, self.reactivity)
        self._ode.set_jac_params(self.material, self.reactivity)
        self._ode.set_integrator("vode", order=1, max_step=self.time_step)

        # record initial power in vector
        self.power.add_data_point(0, 0.0, y0[0])

        # open output file
        fh = open("output.dat", "w")

        # perform integration
        for i in range(self.num_time_steps):
            self._ode.integrate(self._ode.t+self.time_step)
            self.power.add_data_point(i+1, self._ode.t, self._ode.y[0])
            print("{0} {1}".format(self._ode.t, self._ode.y[0]))
            fh.write("{0} {1}\n".format(self._ode.t, self._ode.y[0]))
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

        #######################################################################

    def _steady_state(self):

        # allocate numpy vector
        y0 = np.zeros((self.material.num_precs + 1))

        # fill in the initial power
        y0[0] = self.initial_power

        # loop around precursors
        for i in range(self.material.num_precs):

            # calculate steady state precs
            y0[i+1] = self.material.beta[i] * self.initial_power / \
                      (self.material.decay[i] * self.material.pnl)

        return y0
    ##

    ###########################################################################

# Function routine
def _f(time, y, material, reactivity):

    # get reactivity
    rho = reactivity.interpolate(time)

    # allocate vector
    f = np.zeros((material.num_precs + 1))

    # calculate power row
    f[0] = (rho - material.beta_total) / material.pnl * y[0] + \
           np.sum(material.decay * y[1:material.num_precs+1])

    # calculate precursors
    for i in range(material.num_precs):
        f[i+1] = -material.decay[i]*y[i+1] + material.beta[i]*y[0] /\
                  material.pnl

    return f

    ###########################################################################

# Jacobian routine
def _jac(time, y, material, reactivity):

    # get reactivity
    rho = reactivity.interpolate(time)

    # allocate matrix
    jac = np.zeros((material.num_precs + 1, material.num_precs + 1))

    # set power element
    jac[0,0] = (rho - material.beta_total) / material.pnl

    # loop around precursors
    for i in range(material.num_precs):

        # calculate power row
        jac[0,i+1] = material.decay[i]

        # calculate diagonal
        jac[i+1,i+1] = -material.decay[i]

        # calculate first column
        jac[i+1,0] = material.beta[i] / material.pnl

    return jac
