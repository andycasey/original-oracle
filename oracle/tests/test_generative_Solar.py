# coding: utf-8

""" Infer the properties of the Sun using the GenerativeModel class. """

import logging
import os
import unittest
import matplotlib.pyplot as plt

from glob import glob
from textwrap import dedent

import specutils
import models
import plot

logger = logging.getLogger("oracle")

class Infer_Solar_GenerativeModel(unittest.TestCase):

    output_prefix = "solar-generative"

    def setUp(self):

        # Load the data.
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/")
        self.data = map(specutils.Spectrum1D.load, 
            glob(os.path.join(path, "Solar_*_blue.fits")))

        assert len(self.data) > 0

        configuration = dedent("""
            model:
              redshift: no
              outliers: no
              continuum: no
              doppler_broadening: yes 

            settings:
              max_synth_threads: 10

            mask:
              - [4860.95, 4861.70]

            priors:
              teff: normal(5800, 150)
            """)
        self.model = models.GenerativeModel(configuration)
        

    def runTest(self):
        
        # Make some initial guess of the model parameters based on the prior
        # distributions and some suitable approximations.
        # NB: This is performed automatically by star.optimise if we supply nothing,
        #     but we would like to keep track of what the initial guess was.
        initial_theta = self.model.initial_guess(self.data)
        print("Initial guess for model parameters:")
        for parameter, value in initial_theta.iteritems():
            print("{0}: {1:.2f}".format(parameter, value))

        raise a
        # Create a figure showing the initial point.
        fig = plot.comparison(self.data, self.model, theta=initial_theta)
        fig.savefig("{0}-initial.pdf".format(self.output_prefix))
        
        # Perform the optimisation from the initial guess point.
        optimised_theta = self.model.optimise(self.data, initial_theta)

        # Create a figure showing the optimised point
        optimised_model_spectra = self.model(
            dispersions=[spectrum.disp for spectrum in self.data],
            **optimised_theta)

        fig, axes = plt.subplots(len(self.data), figsize=(25, 2*len(self.data)))
        axes = [axes] if len(self.data) == 1 else axes
        for ax, observed_spectrum, model_spectrum in zip(axes, self.data,
            optimised_model_spectra):

            ax.plot(observed_spectrum.disp, observed_spectrum.flux, 'k')
            ax.plot(model_spectrum[:, 0], model_spectrum[:, 1], 'b')
            ax.set_xlabel("Wavelength, $\lambda$ ($\AA$)")
            ax.set_ylabel("Flux, $F_\lambda$")

        fig.savefig("{0}-optimised.pdf".format(self.output_prefix))

        raise a




if __name__ == "__main__":
    foobar = Infer_Solar_GenerativeModel()
    foobar.setUp()
    foobar.runTest()