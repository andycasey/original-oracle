# coding: utf-8

from __future__ import absolute_import, print_function

""" Theremin Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import logging

import numpy as np
from scipy import optimize as op

from oracle import si, specutils, utils
from .model import Model, _log_prior_eval_environment

logger = logging.getLogger("oracle")


class ThereminModel(Model):

    def __init__(self, configuration):
        """
        A class to probabilistically model stellar spectra. This class performs
        on-the-fly synthesis and fits individual absorption lines in order to
        perform an excitation and ionisation balance, while accounting for
        blends.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, or a string containing the configuration. Only YAML-style
            formats are accepted.

        :type configuration:
            str
        """
        super(ThereminModel, self).__init__(configuration)
        return None


    @property
    def parameters(self):
        """
        Return the parameters of the model.
        """

        return ("teff", "logg", "[M/H]", "xi")


    def __call__(self, **theta):

        # For this set of stellar parameters, what are the expected line
        # abundances?
        # This just needs the atomic lines used for SP determination + some
        # solar abundance reference.

        return None

        
    def initial_guess(self):
        """
        Generate an initial guess of the model parameters.
        """

        environment = dict(zip(["locals", "globals", "__name__", "__file__",
            "__builtins__"], [None] * 5))
        environment["np"] = np

        default_rules = collections.OrderedDict([
            ("teff", "np.random.uniform(3000, 7000)"),
            ("logg", "np.random.uniform(0, 5)"),
            ("[M/H]", "np.random.uniform(-2, 0)"),
            ("xi", "(1.28 + 3.3e-4 * (teff - 6000) - 0.64 * (logg - 4.5)) "\
                "if logg > 3.5 else 2.70 - 0.509 * logg")
        ])

        parameter_guess = {}
        for parameter, rule in default_rules.iteritems():

            local_environment = environment.copy()
            local_environment.update(parameter_guess)
            parameter_guess[parameter] = rule if isinstance(rule, (int, float)) \
                else eval(rule, local_environment)

        return [parameter_guess[p] for p in self.parameters]


    def optimise(self, data):
        """
        Optimise the model parameters given some data.

        :param data:
            The observed stellar spectra.

        :type data:
            list of :class:`specutils.Spectrum1D` objects
        """

        # At each stellar parameter test we need to re-fit all the lines
        # So we need another model class.

        model_spectra = ThereminSpectrumModel(self.configuration, data)

        initial_theta = model_spectra.initial_guess()

        def min(teff, logg, mh, xi):

            model_data = model_spectra.optimise(teff, logg, mh, xi)
            expected = self(teff, logg, mh, xi)

            # Use jacobian approximation for 




        return None



class ThereminSpectrumModel(Model):

    def __init__(self, configuration, data):
        super(ThereminSpectrumModel, self).__init__(configuration)

        if not isinstance(data, (tuple, list)):
            raise TypeError("data must be a list-type of Spectrum1D objects")

        self.data = data
        self._optimised_theta = None
        return None


    @property
    def parameters(self):
        """
        Return the model parameters.
        """

        try:
            return self._parameters
        except AttributeError:
            None

        parameters = []
        num_channels = len(data)
        for transition in self.config["ThereminModel"]["atomic_lines"]:
            wavelength, species = transition[:2]

            # No point adding the transition if it is outside of our observed
            # spectral range.
            for i, spectrum in enumerate(data):
                if (spectrum.disp[-1] >= wavelength >= spectrum.disp[0]):
                    break
            else:
                # This transition is not in any of the observed spectral ranges
                # so let's ignore it.
                continue

            parameters.append("log({0}_{1:.0f})".format(utils.element(species),
                wavelength))

        # Dumb check.
        assert len(parameters) == len(set(parameters)), "At least one atomic "\
            "transition was listed more than once."

        # Any redshift parameters?
        if self.config["model"]["redshift"]:
            # [TODO] Allow individual channel redshifts
            parameters.extend(["z_{0}".format(i) for i in range(len(data))])

        # Any continuum parameters?
        if self.config["model"]["continuum"]:
            # Don't create unnecessary parameters. Only create parameters for
            # channels that are defined *and* have data.
            num_channels_with_continuum = min([
                len(self.config["continuum"]["order"]),
                num_channels
            ])

            for i, each in enumerate(self.config["continuum"]["order"]):
                parameters.extend(["c_{0}_{1}".format(i, j) \
                    for j in range(num_channels_with_continuum)])

        # Any doppler broadening?
        if config["model"]["doppler_broadening"]:
            # [TODO] You idiot.
            parameters.extend(["doppler_sigma_{0}".format(i) \
                for i in range(num_channels)])

        self._parameters = parameters
        return parameters


    def __call__(self, effective_temperature, surface_gravity, metallicity, xi,
        si_kwargs=None, **theta):
        """
        Produce some synthetic spectra to compare against the data.
        """

        expected_fluxes = [np.ones((len(spectrum.disp))) for spectrum in data]

        # Synthesise the required spectra around each line.
        for i, (observed, expected_flux) in enumerate(zip(self.data, expected_fluxes)):

            # For each line:
            # Synthesise some spectrum given the abundance.


            # Apply continuum, redshift, and doppler broadening,
            # The continuum enters multiplicatively with the expected fluxes.
            if self.config["model"]["continuum"]:
                
                # Apply continuum transformation to the model spectra.
                j, coefficients = 0, []
                while "c_{0}_{1}".format(i, j) in theta:
                    coefficients.append(theta["c_{0}_{1}".format(i, j)])
                    j += 1

                if len(coefficients) > 0:
                    expected_flux *= np.polyval(coefficients, observed.disp)

            # Redshift needs to be done before convolution because we observe
            # certain lines at redshifted wavelengths, so we must place those
            # lines at the right wavelengths so they are convolved correctly.
            if self.config["model"]["redshift"]:
                z = theta["z"] if "z" in theta else theta["z_{0}".format(i)]
                redshifted_disp = (1. + z) * observed.disp

                # Interpolate the fluxes onto the observed dispersion.
                expected_flux = np.interp(observed.disp, redshifted_disp,
                    expected_flux, left=np.nan, right=np.nan)

            if self.config["model"]["doppler_broadening"]:
                sigma = theta["doppler_sigma_{0}".format(i)]

                # Convolve the expected fluxes
                # [TODO] This should be done by resolution
                raise NotImplementedError

        return expected_fluxes


    def optimise(self, effective_temperature, surface_gravity, metallicity, xi,
        initial_guess=None, use_previously_optimised=True):
        """
        Given some fixed stellar parameters, optimise the model parameters.
        """

        if initial_guess is None and use_previously_optimised \
        and self._optimised_theta is not None:
            logger.info("Using previously optimised theta as initial guess:"\
                "{0}".format(self._optimised_theta))
            initial_guess = self._optimised_theta


        # Save the optimised
        return None



