# coding: utf-8

""" Probabilistic Inference for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
from functools import partial
from scipy import stats

# [TODO] remove these lines
import matplotlib.pyplot as plt
from glob import glob

import utils

# Initialise logger
logger = logging.getLogger("unnamed")
logger.setLevel(logging.INFO)

# Specify a prior environment so that we can evaluate distributions on-the-fly.
prior_environment = dict(zip(["locals", "globals", "__name__", "__file__",
    "__builtins__"], [None]*5))
prior_environment.update({
    "uniform": lambda a, b: partial(stats.uniform.logpdf, **{"loc": a, "scale": b - a}),
    "normal": lambda a, b: partial(stats.norm.logpdf, **{"loc": a, "scale": b})
})

def log_prior(theta, model):
    """
    Return the logarithm of prior probability for the model parameters 
    :math:`\\theta`.

    :param theta:
        The list of values for the ``model.parameters`` :math:`\\theta`.

    :type theta:
        list

    :param model:
        The stellar spectrum model.

    :type model:
        :class:`models.Spectrum`
    """

    prior = 0
    configuration = model._configuration.get("priors", False)
    for parameter, value in zip(model.parameters, theta):
        if configuration and parameter in configuration:
            prior += eval(configuration[parameter], prior_environment)
            continue

        # Check for outlier parameters:
        if parameter == "Po" and not (1. > value > 0.) \
        or parameter == "Vo" and 0 > value:
            prior = -np.inf
            break

    logger.debug("Returning prior of {0} for {1}".format(
        prior, ", ".join(["{0} = {1}".format(k, v) \
            for k, v in zip(model.parameters, theta)])))
    return prior


def log_likelihood(theta, model, data):
    """
    Return the logarithm of likelihood for the model parameters :math:`\\theta`
    given the data.

    :param theta:
        The list of values ofr the ``model.parameters`` :math:`\\theta`.

    :type theta:
        list

    :param model:
        The stellar spectrum model.

    :type model:
        :class:`models.Spectrum`

    :param data:
        A list of the observed spectra.

    :type data:
        list of :class:`specutils.Spectrum1D` objects

    :returns:
        The logarithm of likelihood of the model parameters :math:`\\theta`
        given the data.

    :rtype float:
    """

    theta_dict = dict(zip(model.parameters, theta))

    try:
        model_spectra, model_continua = model(return_continuum=True, **theta_dict)
    except:
        return -np.inf

    likelihood = 0
    for i, (model_spectrum, continuum, observed_spectrum) \
    in enumerate(zip(model_spectra, model_continua, data)):

        expected = model_spectrum[:, 1]
        z = theta_dict.get("z", theta_dict.get("z_{0}".format(i), 0))
        mask = model.mask(observed_spectrum.disp, z)
        chi_sq = (observed_spectrum.flux - expected)**2 * observed_spectrum.ivariance * mask

        if "Po" not in theta_dict:
            likelihood += -0.5 * np.nansum(chi_sq)

        else:
            Po, Ys = theta_dict["Po"], theta_dict["Ys"]
            
            model_likelihood = -0.5 * chi_sq
            outlier_ivariance = 1.0/(theta_dict["Vo"] + observed_spectrum.variance) # No jitter (yet).
            outlier_likelihood = -0.5 * ((observed_spectrum.flux - Ys * continuum)**2 \
                * outlier_ivariance - np.log(outlier_ivariance)) * mask

            likelihood += np.nansum(np.logaddexp(
                np.log(1. - Po) + model_likelihood,
                np.log(Po)      + outlier_likelihood))

    logger.debug("Returning likelihood of {0} for {1}".format(
        likelihood, ", ".join(["{0} = {1:.2f}".format(k, v) \
            for k, v in theta_dict.iteritems()])))

    if likelihood == 0:
        raise ValueError("No pixels sampled in likelihood. Posterior = Prior.")
    return likelihood


def log_probability(theta, model, data):
    """
    Return the logarithm of probability for the model parameters :math:`\\theta`
    given the data.

    :param theta:
        The list of values for the ``model.parameters`` :math:`\\theta`.

    :type theta:
        list

    :param model:
        The stellar spectrum model.

    :type model:
        :class:`models.Spectrum`

    :param data:
        A list of the observed spectra.

    :type data:
        list of :class:`specutils.Spectrum1D' objects

    :returns:
        The logarithm of probability of the model parameters :math:`\\theta`
        given the data.

    :rtype float:
    """

    prior = log_prior(theta, model)
    if not np.isfinite(prior):
        return prior
    probability = prior + log_likelihood(theta, model, data)
    logger.debug("Returning probability of {0} for {1}".format(
        probability, ", ".join(["{0} = {1}".format(k, v) \
            for k, v in zip(model.parameters, theta)])))
    return probability