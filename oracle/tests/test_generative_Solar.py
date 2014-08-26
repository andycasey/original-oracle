# coding: utf-8

""" Infer the properties of the Sun using the GenerativeModel class. """

import logging
import os
import unittest

from glob import glob
from textwrap import dedent

import specutils
import models

logger = logging.getLogger("unnamed")

class Infer_Solar_GenerativeModel(unittest.TestCase):

    def setUp(self):

        # Load the data.
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/")
        self.data = map(specutils.Spectrum1D.load, 
            glob(os.path.join(path, "Solar_*_blue.fits")))

        assert len(self.data) > 0

        configuration = dedent("""
            model:
              redshift: yes
              outliers: yes
              continuum: no
              doppler_broadening: yes 
            """)
        self.model = models.GenerativeModel(configuration)
        print("Model parameters: {0}".format(", ".join(self.model.parameters)))


    def runTest(self):
        
        # Make some initial guess of the model parameters based on the prior
        # distributions and some suitable approximations.
        # NB: This is performed automatically by star.optimise if we supply nothing,
        #     but we would like to keep track of what the initial guess was.
        initial_guess = self.model.initial_guess(self.data)
        print("Initial guess for model parameters:")
        for parameter, value in initial_guess.iteritems():
            print("{0}: {1:.2f}".format(parameter, value))

        # Create figure showing initial guess?

        # Perform the optimisation from the initial guess point.
        optimised_theta = star.optimise(spectra, initial_guess)

        raise a




if __name__ == "__main__":
    foobar = Infer_Solar_GenerativeModel()
    foobar.setUp()
    foobar.runTest()