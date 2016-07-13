#!/usr/bin/env python

# packages
from __future__ import print_function

from scipy.integrate import ode
from scipy.sparse import coo_matrix

import numpy as np
import pkes
import h5py
import matplotlib.pyplot as plt


class PKEODESolver(object):

    """Point Kinetics solver using ODE from SciPy.
    """

    def __init__(self):

        # initialize private class members
        self._material = None
        self._reactivity = None
        self._power = None
        self._initial_power = 1.0
        self._end_time = None
        self._time = None
        self._max_step = None

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):

        # check material type
        if not isinstance(material, pkes.Material):
            raise TypeError("Solver material is not a PKE material.")

        # set mats
        self._material = material

    @property
    def reactivity(self):
        return self._reactivity

    @reactivity.setter
    def reactivity(self, reactivity):

        # check reactivity types
        if not isinstance(reactivity, pkes.Solution):
            raise TypeError("Reactivity is not a PKE solution.")

        # set reactivity
        self._reactivity = reactivity

    @property
    def time(self):
        return self._time

    @property
    def power(self):
        return self._power

    @property
    def initial_power(self):
        return self._initial_power

    @initial_power.setter
    def initial_power(self, initial_power):

        # check initial power type
        if not isinstance(initial_power, float):
            raise TypeError("Initial power is not a float.")

        # set initial power
        self._initial_power = initial_power

    @property
    def end_time(self):
        return self._end_time

    @end_time.setter
    def end_time(self, end_time):

        # check end time type
        if not isinstance(end_time, float):
            raise TypeError("End time is not a float.")

        # set end time
        self._end_time = end_time

    @property
    def max_step(self):
        return self._max_step

    @max_step.setter
    def max_step(self, max_step):

        # check max step type
        if not isinstance(max_step, float):
            raise TypeError("Max step is not a float.")

        # set end time
        self._max_step = max_step

    def solve(self):

        self._power = pkes.Solution(1)

        # check object
        self._validate()

        # open output file
        fh = open("power.dat", "w")
        
        # calculate steady state
        y0 = self._steady_state()
        fh.write("{} {} {} {}\n".format(
            0, 0.0, self.reactivity.interpolate(0.0), y0[0]))
        self.power.add_data_point(0, 0.0, y0[0])

        # resolve function handles
        create_f = lambda t, y: self._create_f(self, t, y)
        create_jac = lambda t, y: self._create_jac(self, t, y)
        set_solution = lambda t, y: self._set_solution(self, t, y)

        # create ode solver
        ode_solver = ode(create_f, create_jac).set_integrator(
            'vode', method='adams', rtol=1.e-5, atol=1.e-6, first_step=1.e-4,
            min_step=1.e-5, max_step=self._max_step)

        # set ode parameters and initial value
        ode_solver.set_initial_value(y0, 0.0)

        # integrate to end time
        i = 1
        while ode_solver.successful() and ode_solver.t < self.end_time:
            ode_solver.integrate(self.end_time, step=True)
            self.power.add_data_point(i, ode_solver.t, ode_solver.y[0])

            # print to screen and write to file
            fh.write("{} {} {} {}\n".format(
                i, ode_solver.t, self.reactivity.interpolate(ode_solver.t),
                ode_solver.y[0]))

            i += 1

        fh.close()

        # write hdf5 data file
        fh = h5py.File("power.h5", "w")
        fh["time"] = self.power.time
        fh["power"] = self.power.data
        fh.close()

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

        # max step
        if self.max_step is None:
            raise ValueError("Max step not set in PKESolver.")

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

    @staticmethod
    def _create_jac(self, t, y):

        # create sparse matrix vectors
        row = np.zeros(3*self.material.num_precs + 1)
        col = np.zeros(3*self.material.num_precs + 1)
        val = np.zeros(3*self.material.num_precs + 1)

        # reactivity
        rho = self.reactivity.interpolate(t)

        offset = 0

        # power by power
        row[offset] = 0
        col[offset] = 0
        val[offset] = (rho - self.material.beta_total) / self.material.pnl
        offset += 1

        # power by prec
        for i in xrange(self.material.num_precs):
            row[offset] = 0
            col[offset] = i + 1
            val[offset] = self.material.decay[i]
            offset += 1	

        # prec by power
        for i in xrange(self.material.num_precs):
            row[offset] = i + 1
            col[offset] = 0
            val[offset] = self.material.beta[i] / self.material.pnl
            offset += 1

        # prec by prec
        for i in xrange(self.material.num_precs):
            row[offset] = i + 1
            col[offset] = i + 1
            val[offset] = -self.material.decay[i]
            offset += 1

        # return sparse matrix
        jac = coo_matrix((val, (row, col)),
            shape=(self.material.num_precs + 1, self.material.num_precs + 1)).toarray()

        return jac

    @staticmethod
    def _create_f(self, t, y):

        # create numpy array
        f = np.zeros(self.material.num_precs+1)

        # reactivity
        rho = self.reactivity.interpolate(t)

        # power
        f[0] = (rho - self.material.beta_total) / self.material.pnl * y[0] + \
            np.dot(self.material.decay, y[1:])

        # precursors
        for i in xrange(self.material.num_precs):
            f[i+1] = self.material.beta[i] / self.material.pnl * y[0] - \
                self.material.decay[i] * y[i+1]

        return f

    @staticmethod
    def _set_solution(self, t, y):
        self._solution.append([t, y])
