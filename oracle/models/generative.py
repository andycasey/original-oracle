# coding: utf-8

from __future__ import absolute_import, print_function

""" Generative Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import cPickle as pickle
import collections
import logging

import numpy as np
from scipy import optimize as op

from oracle import specutils, utils
from oracle.models.model import Model

logger = logging.getLogger("oracle")


# Process is:
#model = models.GenerativeModel("my_filename.yaml")
# model filename contains:
# masks to use
# normalisation to use (in terms of a list of rules per observed channel)
# redshift to apply: single/different redshifts per channel
# instrumental broadening: yes/no since this is assumed to be different per arm
# radiative transfer code to use
# line list location to use
# model atmospheres to use.
# elemental abundances to consider

# model parameters are determined when data is fed to it, since it depends
# on *how* many channels are provided,
#data = []
#initial_theta, initial_r_chi_sq, info = model.initial_theta(data)
#optimised_theta, model.fit(data, initial_theta=What)

def load_grid(filename):
    with open(filename, "rb") as fp:
        grid_description, grid_points, grid_dispersion, grid_fluxes = \
            pickle.load(fp)
        
    # Reshape the grid fluxes accordingly
    grid_fluxes = grid_fluxes.reshape(-1, grid_dispersion.size)
    return (grid_description, grid_points, grid_dispersion, grid_fluxes)


class GenerativeModel(Model):

    """
    A class to forward model stellar spectra. This class does multi-dimensional
    interpolation of model spectra and on-the-fly synthesis to generate stellar
    spectra for some given set of model parameters theta.
    """

    def __init__(self, configuration):
        """
        Initialise a GenerativeModel class.
        
        :param configuration:
            The path to the configuration file that contains settings for the
            class, a dictionary describing the configuration, or a string 
            containing the configuration in a YAML-dump format. Only YAML-style
            formats are accepted.

        :type configuration:
            str or dict
        """

        super(GenerativeModel, self).__init__(configuration)
        return None


    def parameters(self, num_data_channels):
        """
        Return the model parameters for some data. The model configuration is
        applicable for a large number of observed channels, so the number of 
        observed channels is required to determine the exact number of model
        parameters.

        :param num_data_channels:
            The number of observed data channels.

        :type num_data_channels:
            int

        :returns:
            A list of model parameters.

        :rtype:
            tuple
        """

        try:
            num_data_channels = int(num_data_channels)
        except (ValueError, TypeError):
            raise TypeError("number of data channels must be an integer")

        if 1 > num_data_channels:
            raise ValueError("number of data channels must be a positive integer")

        parameters = ["effective_temperature", "surface_gravity", "metallicity",
            "microturbulence"]

        # Single radial velocity for all channels
        if self.config["model"]["redshift"] == True:
            parameters.append("v_rad")

        # Different radial velocity for each channel?
        elif isinstance(self.config["model"]["redshift"], (tuple, list, )):
            parameters.extend(["v_rad.{}".format(i) for i, channel_v_rad in \
                zip(range(num_data_channels), self.config["model"]["redshift"])\
                if channel_v_rad])

        # Instrumental broadening
        if self.config["model"]["instrumental_resolution"]:
            parameters.extend(["instrumental_resolution.{}".format(i) \
                for i in range(num_data_channels)])

        # Continuum treatment
        if isinstance(self.config["model"]["continuum"], (tuple, list)):
            # List contains order for each channel
            for i, order in \
                zip(range(num_data_channels), self.config["model"]["continuum"]):
                parameters.extend(["continuum.{0}.{1}".format(i, j) \
                    for j in range(order + 1)])

        return tuple(parameters)


    def initial_theta(self, data, grid_filename):
        """
        Return an initial guess of the model parameters theta using no prior
        information.

        :param data:
            The observed data.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects

        :param grid_filename:
            A filename pointing to a pickled grid of model spectra.

        :type grid_filename:
            str
        """

        # Need some cacher object to do fast-on-the-fly-comparisons
        parameters = self.parameters(len(data))

        # Load the pickled cacher object thing, and work out some estimates
        # of the model parameters
        # Cacher must contain:
        # dict with info about the grid, etc., grid points record array, 
        # dispersion map, gigantic flux array
        grid_description, grid_points, grid_dispersion, grid_fluxes = \
            load_grid(grid_filename)

        theta = {}
        continuum_coefficients = []
        chi_sqs = np.zeros(grid_points.size)

        # Parallelise channels
        logger.warn("andy parallelise this part")
        for i, channel in enumerate(data):

            # Splice the wavelengths
            indices = grid_dispersion.searchsorted([channel.disp[0],
                channel.disp[-1]])

            # Temporarily transform the data to the model dispersion points
            rebinned_channel_disp = grid_dispersion[indices[0]:indices[1]]
            rebinned_channel_flux = np.interp(rebinned_channel_disp,
                channel.disp, channel.flux, left=np.nan, right=np.nan)
            rebinned_channel_ivar = np.interp(rebinned_channel_disp,
                channel.disp, channel.ivariance, left=np.nan, right=np.nan)

            # Reference only finite pixels
            finite = np.isfinite(rebinned_channel_flux)

            # Cross-correlate the observed data against the grid
            if ("v_rad" in parameters) or ("v_rad.{}".format(i) in parameters):
                v_rads, v_errs, ccf_maxes = specutils.cross_correlate_grid(
                    rebinned_channel_flux, grid_fluxes[:, indices[0]:indices[1]])

                # Identify one with largest CCF max
                v_rad = v_rads[ccf_maxes.argmax()]
                if "v_rad" in parameters:
                    # Global, so add it to a list which we will use to take an
                    # average from later
                    if "v_rad" in theta:
                        theta["v_rad"].append(v_rad)
                    else:
                        theta["v_rad"] = [v_rad]
                else:
                    # Measured radial velocity applies to this channel only
                    theta["v_rad.{}".format(i)] = v_rad

            # Calculate continuum coefficients
            try:
                order = self.config["model"]["normalise"][i]
            
            except (KeyError, TypeError):
                continuum = np.ones(finite.size)
            
            else:

                # Calculate continuum coefficients for each model grid point
                dispersion_matrix = np.ones((order + 1, finite.size))
                for j in range(order + 1):
                    dispersion_matrix[j] *= rebinned_channel_disp[finite]**j

                A = (rebinned_channel_flux[finite] \
                    / grid_fluxes[:, indices[0]:indices[1]][:, finite]).T
                coefficients = np.linalg.lstsq(dispersion_matrix.T, A)[0]

                # Save the continuum coefficients (we will use them later)
                continuum_coefficients.append(coefficients)
                continuum = np.dot(coefficients.T, dispersion_matrix)
            
            # Calculate the expected fluxes and the chi-sq value at each point
            expected_fluxes = grid_fluxes[:, indices[0]:indices[1]][:, finite]\
                * continuum

            # Add to the chi-sq values
            chi_sqs += np.nansum(
                (rebinned_channel_flux[finite] - expected_fluxes)**2\
                    * rebinned_channel_ivar[finite], axis=1)

            # Estimate instrumental broadening
            logger.warn("haven't done instrumental broadening")

        # OK, now determine the nearest point as that which has the lowest chi-sq
        grid_index = chi_sqs.argmin()

        # Smash The Cannon at it!



        raise a
        """

            for j, value in enumerate(coefficients[::-1, index]):
        theta["normalise.{0}.b{1:.0f}".format(channel, j)] = value

    # Use values to estimate microturbulence


        ("xi", "(1.28 + 3.3e-4 * (teff - 6000) - 0.64 * (logg - 4.5)) "\
        "if logg > 3.5 else 2.70 - 0.509 * logg"),
            
        """


        # ON A PER CHANNEL BASIS:
        # - Splice the pickled cacher thing:
        # If redshift is modelled for this channel (global or otherwise):
        # - Cross-correlate the observed data against the grid (in parallel)
        # - Use the best-estimate of the CCF as the redshift estimate
        # Then:
        # - Calculate continuum coefficients
        # - Estimate microturbulence from standard relation
        # - Estimate smoothing value from the observed pixel sizes
        # Then:
        # - Estimate any abundance parameters scale with the overall metallicity


    def fit(self, data, initial_theta=None, grid=None):

        # Use a grid to do interpolation if provided, otherwise we will have to
        # do radiative transfer on the fly


        return None



