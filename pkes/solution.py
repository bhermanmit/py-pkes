#!/usr/bin/env python

# packages
from __future__ import print_function
import numpy as np

###############################################################################

class Solution(object):

    ###########################################################################
    ## Constructor
    ##
    def __init__(self, size):

        # initialize private class members
        self._size = None
        self._data = None
        self._time = None

        # check size
        if not isinstance(size, int):
            raise TypeError("Size is not an integer.")

        # set size
        self._size = size

        # allocate vectors
        self._data = np.empty((size))
        self._time = np.empty((size))
    ##

    ###########################################################################
    ## Properties
    ##
    def get_size(self):
        return self._size

    size = property(get_size)

        #######################################################################

    def get_data(self):
        return self._data

    data = property(get_data)

        #######################################################################

    def get_time(self):
        return self._time

    time = property(get_time)

    ##

    ###########################################################################
    ## Methods
    ##
    def add_data_point(self, index, time, data):

        # check if index is an integer
        if not isinstance(index, int): 
            raise TypeError("Index is not an integer.")

        # check if time is a float
        if not isinstance(time, float):
            raise TypeError("Time is not a float.")

        # check if data is a float
        if not isinstance(data, float):
            raise TypeError("Data is not a float.")

        # check if index is valid
        if index < 0:
            raise IndexError("Index not in range of solution.")
            
        # extend if necessary
        if index >= self.size:
            self._time.resize(index+1)
            self._data.resize(index+1)
            self._size = index+1

        # set parameters
        self._time[index] = time
        self._data[index] = data

        #######################################################################

    def interpolate(self, time):

        # check if time is a float
        if not isinstance(time, float):
            raise TypeError("Time is not a float for interpolation.")

        # perform interpolation
        return np.interp(time, self.time, self.data)

        #######################################################################

    def exp_interpolate(self, time):

        # check if time is a float
        if not isinstance(time, float):
            raise TypeError("Time is not a float for interpolation.")

        # perform interpolation
        data = np.interp(time, self.time, np.log(self.data))
        return np.exp(data)

        #######################################################################

    def __repr__(self):

        the_str = "Solution: size={0}\n\n".format(self._size)

        # set headers
        the_str += "    Time        Data    \n"

        # loop around data 
        for time, data in zip(self.time, self.data):
            the_str += "{0:12.5e} {1:12.5e}\n".format(time, data)

        return the_str
    ##
