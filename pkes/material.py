#!/usr/bin/env python

# packages
from __future__ import print_function
import numpy as np

###############################################################################

class Material(object):

    ###########################################################################
    ## Constructor
    ##
    def __init__(self):

        # initialize private class members
        self._beta = None
        self._beta_total = None
        self._decay = None
        self._num_precs = None
        self._pnl = None
    ##

    ###########################################################################
    ## Properties
    ##
    def get_beta(self):
        return self._beta

    def set_beta(self, beta):

        # check that beta is an array
        if not isinstance(beta, (list, tuple, np.ndarray)):
            raise TypeError("Beta must be a list or numpy array.")

        # check that all elements are floats
        for b in beta:
            if not isinstance(b, float):
                raise TypeError("Element in beta not a float.")

        # check number of precursors
        if not self._num_precs is None:
            if len(beta) != self._num_precs:
                raise IndexError("Beta vector is not of num_precs length.")
        else:
            self._num_precs = len(beta)

        # store beta as numpy array
        self._beta = np.asarray(beta)

        # store total beta
        self._beta_total = np.sum(self.beta)

    beta = property(get_beta, set_beta)

        #######################################################################

    def get_beta_total(self):
        return self._beta_total

    beta_total = property(get_beta_total)

        #######################################################################

    def get_decay(self):
        return self._decay

    def set_decay(self, decay):

        # check that decay is an array
        if not isinstance(decay, (list, tuple, np.ndarray)):
            raise TypeError("Decay must be a list or numpy array.")

        # check that all elements are floats
        for d in decay:
            if not isinstance(d, float):
                raise TypeError("Element in decay not a float.")

        # check number of precursors
        if not self._num_precs is None:
            if len(decay) != self._num_precs:
                raise IndexError("Decay vector is not of num_precs length.")
        else:
            self._num_precs = len(decay)

        # store decay as numpy array
        self._decay = decay

    decay = property(get_decay, set_decay)

        #######################################################################

    def get_num_precs(self):
        return self._num_precs

    num_precs = property(get_num_precs)

        #######################################################################

    def get_pnl(self):
        return self._pnl

    def set_pnl(self, pnl):

        # check that pnl is a float
        if not isinstance(pnl, float):
            raise TypeError("Pnl is not a float.")

        # set pnl
        self._pnl = pnl
    ##

    ###########################################################################
    ## Methods
    ##
    def __repr__(self):

        the_str = "PKE Material\n"

        # number of precursors
        the_str += "    Number of precursors: {0}\n".format(self._num_precs)

        # prompty neutron lifetime
        the_str += "    Prompt neutron lifetime: {0}\n".format(self.pnl)

        # loop around precursor groups
        the_str += "    Total beta: {0}\n".format(self.beta_total)
        for i in range(self._num_precs):
            the_str += "        Group: {0:2d}  Beta: {1:7.5f}  " \
                       "Decay: {2:7.5f}\n".format(i+1, self.beta[i],
                       self.decay[i])

        return the_str
    ##
