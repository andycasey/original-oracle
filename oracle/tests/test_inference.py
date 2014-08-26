# coding: utf-8

""" Self-consistent inference test """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import unittest

import si
import models
import specutils

class InferenceTest(unittest.TestCase):

	def setUp(self):
		""" Create a fake noisy spectrum to use for the inference test. """

		self.star = models.GenerativeModel(os.path.join(os.path.dirname(__file__), "config.yaml"))

		truth = {
			"Teff": 4503,
			"logg": 2.104,
			"log_Fe": -1.403,
			"xi": 2.10,
			"z": 0.
			# Fake some continuum coefficients.

		}

		# Let's replicate the HERMES band passes.
		hermes_wavelengths = [
			(4718, 4903),
			(5649, 5873),
			(6481, 6739),
			(7590, 7890)
		]
		self.data = [specutils.Spectrum1D(disp=spectrum[:, 0], flux=spectrum[:, 1], \
			variance=np.ones(len(spectrum))*0.001) for spectrum in \
			self.star(wavelength_ranges=hermes_wavelengths, **truth)]

		"""
		for channel in self.data:
			N = len(channel)
			# [TODO] Check these calculations.
			flux_err = (0.1 + 0.1 * np.random.randn(N)) * channel[:, 1]
			channel[:, 1] += flux_err * np.random.randn(N)
		

		fig, axes = plt.subplots(len(self.data))
		if len(self.data) == 1: axes = [axes]
		for ax, channel in zip(axes, self.data):
			ax.plot(channel[:, 0], channel[:, 1], 'k')
			ax.set_ylabel("Flux, $F_\lambda$")
			ax.set_yticklabels([])
			ax.set_xlim(channel[0, 0], channel[-1, 0])
			ax.set_xlabel("Wavelength, $\lambda$ [$\AA$]")
		fig.savefig("spectrum.pdf")
		"""


	def runTest(self):
		""" Infer the model parameters given the fake data. """

		result = self.star.optimise(self.data)


	def tearDown(self):
		pass