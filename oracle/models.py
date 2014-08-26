# coding: utf-8

""" Probabilistic Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import logging
import multiprocessing as mp
import numpy as np
import os
import scipy.optimize as op
import scipy.stats as stats
import yaml

from functools import partial
from scipy.ndimage import gaussian_filter1d
from time import time


import emcee

import inference
import si
import specutils
import utils

__all__ = ["Star", "StellarSpectrum", "GenerativeModel"]

logger = logging.getLogger("oracle")
logger.setLevel(logging.DEBUG)

_eval_environment_ = { 
    "locals": None, "globals": None, "__name__": None, "__file__": None, "__builtins__": None,
    "uniform": lambda a, b: partial(stats.uniform.logpdf, **{"loc": a, "scale": b - a}),
    "normal": lambda a, b: partial(stats.norm.logpdf, **{"loc": a, "scale": b})
}


class AbsorptionProfile(object):

    def __init__(self, wavelength, method="gaussian", mask=None, outliers=False,
        wavelength_tolerance=0.0, wavelength_contribution=0.2):
        """
        Model an absorption profile in a spectrum.

        :param wavelength:
            Approximate wavelength of the absorption feature.

        :type wavelength:
            float

        :param method: [optional]
            The type of profile to use. Available profiles are Gaussian or Voigt.

        :type method:
            str

        :param mask: [optional]
            A mask to use for the observed channel. This should be of the same
            length as the observed data points.

        :type mask:
            :class:`numpy.array`

        :param outliers: [optional]
            Model outlier pixels with a Gaussian mixture model.

        :type outliers:
            bool

        :wavelength_tolerance: [optional]
            Allow some tolerance in the exact rest wavelength of the absorption
            profile. If set to zero, no tolerance is allowed.

        :type wavelength_tolerance:
            float

        :wavelength_contribution: [optional]
            The spectral region to consider around the absorption feature. This
            must be a positive quantity.

        :type wavelength_contribution:
            float

        :param upper_smoothing_sigma: [optional]
            The upper limit allowed for the Gaussian full-width half-maximum of
            the absorption profile. Set as 0 for no limit.

        :type upper_smoothing_sigma:
            float

        :raises ValueError:
            If the ``method`` specified is not Voigt or Gaussian, or if the 
            ``wavelength_tolerance``, ``wavelength_contribution`` values are negative.
        """

        self.method = method.lower()
        if self.method not in ("gaussian", "voigt"):
            raise ValueError("method must be either Gaussian or Voigt")

        if 0 > wavelength_tolerance:
            raise ValueError("wavelength tolerance must be a positive quantity")

        if 0 > wavelength_contribution:
            raise ValueError("wavelength contribution must be a positive quantity")

        self.mask = mask if mask is not None else 1.
        self.outliers = outliers
        self.approx_wavelength = wavelength
        self.wavelength_tolerance = wavelength_tolerance
        self.wavelength_contribution = wavelength_contribution
        
        return None


    @property
    def parameters(self):
        """ Return the model parameters. """

        parameters = ["ld"]
        if self.method == "voigt":
            parameters.extend(["fwhm", "shape"])
        elif self.method == "gaussian":
            parameters.append("sigma")

        if self.wavelength_tolerance > 0:
            parameters.append("wl")
        if self.outliers:
            parameters.extend(["Po", "Vo", "Yo"])
        return parameters


    def __call__(self, dispersion, continuum=1., **theta):
        """
        Return the continuum with the absorption profile for the given theta.

        :param dispersion:
            The dispersion array to calculate the spectrum for.

        :type dispersion:
            :class:`numpy.array`

        :param continuum: [optional]
            The continuum for the spectrum. Defaults to unity at all points.

        :type continuum:
            :class:`numpy.array` or float

        :param theta:
            The model parameters (keys) and their values.

        :type theta:
            dict

        :returns:
            An array containing the dispersion and continuum with the requested
            absorption profile.

        :rtype:
            :class:`numpy.array`
        """

        wavelength = theta.get("wl", self.approx_wavelength)
        if self.method == "voigt":
            depth, shape, fwhm = [theta[p] for p in ("ld", "shape", "fwhm")]
            flux = continuum * (1. - utils.voigt(wavelength, fwhm, depth, shape,
                dispersion))
        elif self.method == "gaussian":
            depth, sigma = theta["ld"], theta["sigma"]
            flux = continuum * (1. - depth * utils.gaussian(wavelength, sigma,
                dispersion))
        return np.vstack([dispersion, flux]).T


    def _log_prior(self, theta):
        """
        Return the logarithmic prior probability of the parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :returns:
            The logarithmic prior probability given the parameters ``theta``.

        :rtype:
            float
        """

        theta_dict = dict(zip(self.parameters, theta))
        if (self.wavelength_tolerance > 0 \
        and abs(self.approx_wavelength - theta_dict["wl"]) > self.wavelength_tolerance) \
        or not (1 > theta_dict.get("Po", 0.5) > 0) or 0 > theta_dict.get("Vo", 1) \
        or not (1 >= theta_dict.get("ld", 0.5) >= 0):
            return -np.inf

        ln_prior = 0
        for parameter, distribution in self._configuration.get("priors", {}).iteritems():
            f = eval(distribution, _eval_environment_)
            ln_prior += f(theta_dict[parameter])
        return ln_prior
        

    def _log_likelihood(self, theta, data, continuum):
        """
        Return the logarithmic likelihood of the parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :param continuum:
            The continuum shape to apply.

        :type continuum:
            float or :class:`numpy.array` with the same length as ``data``

        :returns:
            The logarithmic likelihood of the parameters ``theta``.

        :rtype:
            float
        """

        theta_dict = dict(zip(self.parameters, theta))
        expected = self(dispersion=data.disp,
            continuum=continuum, **theta_dict)[:, 1]

        # Build an extra mask if only nearby points contribute to the model.
        if self.wavelength_contribution > 0:
            wavelength = theta_dict.get("wavelength", self.approx_wavelength)
            indices = data.disp.searchsorted([
                wavelength - self.wavelength_contribution,
                wavelength + self.wavelength_contribution
            ])
            additional_mask = np.array([np.nan]*len(data.disp))
            additional_mask[indices[0]:indices[1]] = 1.
        
        else:
            additional_mask = 1.

        total_mask = self.mask * additional_mask
        chi_sq = (data.flux - expected)**2 * data.ivariance * total_mask

        if not self.outliers:
            likelihood = -0.5 * np.nansum(chi_sq)

        else:
            Po, Vo, Yo = [theta_dict.get(each) for each in ("Po", "Vo", "Yo")]
            model_likelihood = -0.5 * (chi_sq - np.log(data.ivariance))
            outlier_ivariance = 1.0/(Vo + data.variance)
            outlier_likelihood = -0.5 * ((data.flux - Yo)**2 * outlier_ivariance \
                * total_mask - np.log(outlier_ivariance))
            likelihood = np.nansum(np.logaddexp(
                np.log(1. - Po) + model_likelihood,
                np.log(Po)      + outlier_likelihood))

        return likelihood


    def _log_probability(self, theta, data, continuum):
        """
        Return the logarithmic probability of the model parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :param continuum:
            The continuum shape to apply.

        :type continuum:
            float or :class:`numpy.array` with the same length as ``data``

        :returns:
            The logarithmic probability of the parameters ``theta``.

        :rtype:
            float
        """

        log_prior = self._log_prior(theta)
        if not np.isfinite(log_prior):
            return log_prior
        return log_prior + self._log_likelihood(theta, data, continuum)
        

    def initial_guess(self, data, continuum):
        """
        Return an initial guess for the model parameters.

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :returns:
            An initial guess for the model parameters.

        :rtype:
            dict
        """

        index = data.disp.searchsorted(self.approx_wavelength)
        continuum_at_line = continuum if isinstance(continuum, (int, float)) \
            else continuum[index]

        theta = {
            "wl": self.approx_wavelength,
            "ld": 1.0 - data.flux[index]/continuum_at_line,
            "fwhm": 0.10,
            "sigma": 0.10,
            "shape": 0.,
            "Po": 0.01,
            "Vo": 0.01,
            "Yo": 1.0
        }
        return dict(zip(self.parameters, [theta[p] for p in self.parameters]))


    def optimise(self, data, continuum, maxfun=10e3, maxiter=10e3, xtol=1e-8,
        ftol=1e-8, full_output=False):
        """
        Optimise the model parameters given the data.

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :param continuum:
            The continuum shape to apply.

        :type continuum:
            float or :class:`numpy.array` with the same length as ``data``

        :param maxfun: [optional]
            The maximum number of function evaluations to perform.

        :type maxfun:
            int

        :param maxiter: [optional]
            The maximum number of function iterations to perform.

        :type maxiter:
            int

        :param xtol: [optional]
            The tolerance required in the model parameters.

        :type xtol:
            float

        :param ftol: [optional]
            The tolerance required in logarithmic probability.

        :type xtol:
            float

        :param full_output: [optional]
            Return a tuple containing the optimised points theta, the logarithmic
            probability of the optimised point, the number of function iterations,
            the number of function evaluations, and an integer warning flag.

        :type full_output:
            bool

        :returns:
            The optimised model parameters.
        """

        initial_theta = self.initial_guess(data, continuum)
        op_kwargs = {
            "maxfun": maxfun,
            "maxiter": maxiter,
            "xtol": xtol,
            "ftol": ftol,
            "disp": False,
            "full_output": True # Necessary for introspection and provenance.
        }
        t_init = time()
        op_theta, op_fopt, op_niter, op_nfunc, op_warnflag = op.fmin(
            lambda t, d, c: -self._log_probability(t, d, c),
            [initial_theta[p] for p in self.parameters],
            args=(data, continuum), **op_kwargs)
        Model._opt_warn_message(op_warnflag, op_niter, op_nfunc)

        t_elapsed = time() - t_init
        logger.info("Optimisation took {0:.2f} seconds".format(t_elapsed))

        if full_output:
            return (op_theta, op_fopt, op_niter, op_funcalls, op_warnflag)
        return op_theta


class Model(object):

    def __init__(self, configuration):
        """
        A general class to probabilistically model stellar spectra.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, a dictionary describing the configuration, or a string 
            containing the configuration in a YAML-dump format. Only YAML-style
            formats are accepted.

        :type configuration:
            str or dict
        """

        if isinstance(configuration, dict):
            # Make a copy of the configuration.
            self._configuration = {}
            self._configuration.update(configuration)

        else:
            if not os.path.exists(configuration):
                # Probably a string configuration.
                try:
                    self._configuration = yaml.load(configuration)
                except:
                    raise IOError("configuration file does not exist or the "\
                        "string configuration is not valid YAML")
            else:
                with open(configuration, "r") as fp:
                    self._configuration = yaml.load(fp)

        # Check the configuration is valid.
        self._validate()

        # Set the masks.
        self.masks = np.array(self._configuration["mask"]) \
            if "mask" in self._configuration else None

        return None


    def _validate(self):
        """ Check that the configuration is valid. """

        if not self._configuration.has_key("model"):
            raise KeyError("no model information specified")

        validation_functions = {
            "continuum": self._validate_continuum,
            "elements": self._validate_elements
        }
        for item, state in self._configuration["model"].iteritems():
            if not state or item not in validation_functions: continue
            validation_functions[item]()

        return True


    def _validate_continuum(self):
        """ Check that the continuum configuration is valid. """
        
        # We actually need *something* specified to model the continuum.
        if not self._configuration.has_key("continuum"):
            raise KeyError("no information specified for continuum modelling")

        order = self._configuration["continuum"]["order"]
        assert self._configuration["continuum"]["method"] == "polynomial"
        assert order

        # Order should be integer or list-like of integers.
        try: order = [int(order)]
        except:
            try: order = map(int, order)
            except:
                raise TypeError("continuum order must be an integer or a list-\
                    like of integers")
        self._configuration["continuum"]["order"] = order
        
        return True


    def _validate_elements(self):
        """ Check that the elements actually exist. """

        map(utils.atomic_number, self._configuration["model"]["elements"])
        return True


    @classmethod
    def _opt_warn_message(cls, warnflag, niter, nfunc):
        """ Log a warning message after optimisation. """

        if warnflag > 0:
            message = [
                "Che problem?",
                "Maximum number of function evaluations ({0}) made".format(nfunc),
                "Maximum number of iterations ({0}) made".format(niter)
            ]
            logging.warn("{0}. Optimised values may be inaccurate.".format(
                message[warnflag]))
        return None


    def mask(self, dispersion, z, fill_value=np.nan):
        """
        Return an array mask for a given dispersion array and redshift, based on
        the mask information provided in the model configuration file.

        :param dispersion:
            An array of dispersion values to mask.

        :type dispersion:
            :class:`numpy.array`

        :param z:
            The redshift to apply.

        :type z:
            float

        :param mask_value: [optional]
            The value to use for the masked dispersion points.

        :type mask_value:
            float-like

        :returns:
            An array of the same length of ``dispersion``, filled with ones,
            except for masked values which contain ``mask_value``.

        :rtype:
            :class:`numpy.array`
        """

        mask = np.ones(len(dispersion))
        if self.masks is not None:
            for start, end in self.masks * (1. + z):
                indices = mask.searchsorted([start, end])
                mask[indices[0]:indices[1] + 1] = fill_value
        return mask


def _optimiser(_class, *args):
    return _class.optimise(*args)


def _inferencer(theta, _class, *args):
    return _class._log_probability(theta, *args)
        

class SpectralChannel(Model):

    wavelength_tolerance = 0.10

    def __init__(self, configuration, channel_index=None,
        absorption_profile_kwargs=None):
        """
        A class to model a stellar spectrum with a combination of continuum and
        absorption profiles. No spectra are synthesised in this model.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, a dictionary describing the configuration, or a string 
            containing the configuration in a YAML-dump format. Only YAML-style
            string formats are accepted.

        :type configuration:
            str or dict

        :param channel_index: [optional]
            An integer index to use for all the channel model parameters.

        :type channel_index:
            integer or str

        :param absorption_profile_kwargs: [optional]
            Keywords to pass directly to :class:`AbsorptionProfile`.

        :type absorption_profile_kwargs:
            dict
        """

        super(SpectralChannel, self).__init__(configuration)
        self._data = None
        self.channel_index = (str(channel_index), "")[channel_index is None]
        self.absorption_profile_kwargs = absorption_profile_kwargs \
            if absorption_profile_kwargs is not None else {}
        return None


    @property
    def parameters(self, data=None):
        """ Return the model parameters. """
        
        if hasattr(self, "_parameters"):
            return self._parameters

        parameters = []
        config = self._configuration
        
        # Any redshift parameters?
        if config["model"]["redshift"]:
            parameters.append("z")

        # Any continuum parameters?
        if config["model"]["continuum"]:
            for i, order in enumerate(config["continuum"]["order"]):
                if (i + 1) == self.channel_index \
                or (len(self.channel_index) == 0 and i == 0):
                    parameters.extend(["c_{0}_{0}".format(self.channel_index, j) \
                        for j in range(order)])

        # Any outlier parameters?
        if config["model"]["outliers"]:
            parameters.extend(["Po", "Vo", "Yo"])

        # Get absorption profile parameters for each transition in this channel.
        absorption_profile_parameters = set(AbsorptionProfile(0,
            **self.absorption_profile_kwargs).parameters).difference(["Po", "Vo", "Yo"])
        
        # We can't get the exact model parameters unless we know the observed
        # dispersion range.
        for i, (wavelength, species) in enumerate(config["balance"]["atomic_lines"]):
            if self._data is None \
            or (self._data.disp[-1] >= wavelength >= self._data.disp[0]):
                parameters.extend(["{0}_{1}".format(each, i) \
                    for each in absorption_profile_parameters])

        if self._data is not None:
            self._parameters = parameters
        return parameters


    def _get_num_coefficients(self):
        try:
            num_coefficients = self._num_coefficients

        except AttributeError:
            self._num_coefficients = 0
            while "c_{0}_{1}".format(self.channel_index, self._num_coefficients) \
            in self.parameters:
                self._num_coefficients += 1
            num_coefficients = self._num_coefficients

        return num_coefficients


    def __call__(self, dispersion, **theta):
        """
        Return a stellar spectrum for the given theta.

        :param dispersion:
            The dispersion array to calculate the spectrum for.

        :type dispersion:
            :class:`numpy.array`

        :param theta:
            The model parameters (keys) and their values.

        :type theta:
            dict

        :returns:
            An array containing the dispersion and flux points.

        :rtype:
            :class:`numpy.array`
        """
        
        num_coefficients = self._get_num_coefficients()

        # Create continuum (but we'll call it flux)
        if num_coefficients > 0:
            coefficients = [theta["c_{0}_{1}".format(self.channel_index, i)] \
                for i in range(num_coefficients)]
            flux = np.polyval(coefficients, dispersion)

        else:
            flux = np.ones(len(dispersion))

        # Model the absorption profiles!
        method = self.absorption_profile_kwargs.get("method", "gaussian")
        for i, (wavelength, species) in enumerate(self._configuration["balance"]["atomic_lines"]):
            # Should we be modelling this line?
            if "ld_{0}".format(i) in self.parameters:
                wavelength = theta.get("wl_{0}".format(i), wavelength)

                indices = dispersion.searchsorted([
                    wavelength - 0.5,
                    wavelength + 0.5
                ])
                x = dispersion.__getslice__(*indices)
                y = flux.__getslice__(*indices)
                
                if method == "voigt":
                    depth, shape, fwhm = [theta["_".join([p, str(i)])] \
                        for p in ("ld", "shape", "fwhm")]
                    flux.__setslice__(indices[0], indices[1],
                        y * (1. - utils.voigt(wavelength, fwhm, depth, shape, x)))

                elif method == "gaussian":
                    depth, sigma = [theta["_".join([p, str(i)])] \
                        for p in ("ld", "sigma")]
                    flux.__setslice__(indices[0], indices[1],
                        y * (1. - depth * utils.gaussian(wavelength, sigma, x)))

        # [TODO] No redshift modelling.
        return np.vstack([dispersion, flux]).T


    def _log_prior(self, theta):
        """
        Return the logarithmic prior probability of the parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :returns:
            The logarithmic prior probability given the parameters ``theta``.

        :rtype:
            float
        """

        theta_dict = dict(zip(self.parameters, theta))

        if not (1 > theta_dict.get("Po", 0.5) > 0) or 0 > theta_dict.get("Vo", 1):
            return -np.inf

        for i, (wavelength, species) in enumerate(self._configuration["balance"]["atomic_lines"]):
            if ("wl_{0}".format(i) in self.parameters \
            and abs(wavelength - theta_dict["wl_{0}".format(i)]) > self.wavelength_tolerance) \
            or not (1 >= theta_dict.get("ld_{0}".format(i), 0) >= 0) \
            or 0 > theta_dict.get("sigma_{0}".format(i), 1) \
            or 0 > theta_dict.get("fwhm_{0}".format(i), 1):
                return -np.inf
        
        ln_prior = 0
        for parameter, distribution in self._configuration.get("priors", {}).iteritems():
            f = eval(distribution, _eval_environment_)
            ln_prior += f(theta_dict[parameter])
        return ln_prior

    def _log_likelihood(self, theta, data):
        """
        Return the logarithmic likelihood of the parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :returns:
            The logarithmic likelihood of the parameters ``theta``.

        :rtype:
            float
        """

        theta_dict = dict(zip(self.parameters, theta))
        expected = self(dispersion=data.disp, **theta_dict)[:, 1]

        z = theta_dict.get("z", 0.)
        mask = self.mask(data.disp, z)

        chi_sq = (data.flux - expected)**2 * data.ivariance * mask
        if "Po" in theta_dict:
            Po, Vo, Yo = [theta_dict.get(each) for each in ("Po", "Vo", "Yo")]
            model_likelihood = -0.5 * (chi_sq - np.log(data.ivariance))
            outlier_ivariance = 1.0/(Vo + data.variance)
            outlier_likelihood = -0.5 * ((data.flux - Yo)**2 * outlier_ivariance \
                * mask - np.log(outlier_ivariance))
            likelihood = np.nansum(np.logaddexp(
                np.log(1. - Po) + model_likelihood,
                np.log(Po)      + outlier_likelihood))

        else:
            likelihood = -0.5 * np.nansum(chi_sq)

        return likelihood


    def _log_probability(self, theta, data):
        """
        Return the logarithmic probability of the model parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :returns:
            The logarithmic probability of the parameters ``theta``.

        :rtype:
            float
        """

        log_prior = self._log_prior(theta)
        if not np.isfinite(log_prior):
            return log_prior
        return log_prior + self._log_likelihood(theta, data)


    def initial_guess(self, data, initial_clip_iterations=5,
        initial_clip_limits=(0.2, 3.0), absorption_profile_parameters=True):
        """
        Return an initial guess for the model parameters.

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :param initial_clip_iterations: [optional]
            The number of sigma-clipping iterations to initially perform to get
            a good starting estimate of the continuum.

        :type initial_clip_iterations:
            int

        :param initial_clip_limits: [optional]
            The lower and upper sigma-clipping limits fo apply to the data.

        :type initial_clip_limits:
            tuple

        :param absorption_profile_parameters: [optional]
            Estimate the absorption profile parameters as well.

        :type absorption_profile_parameters:
            bool

        :returns:
            An initial guess for the model parameters.

        :rtype:
            dict
        """

        # Make an internal reference to the data so that we can get the right
        # model parameters.
        self._data = data

        # Estimate the continuum parameters.
        num_coefficients = self._get_num_coefficients()

        initial_theta = {}
        if num_coefficients > 0:
            mask = np.isfinite(self.mask(data.disp, z))
            continuum_coefficients = np.polyfit(data.disp[mask],
                data.flux[mask], num_coefficients)

            lower_clip, upper_clip = initial_clip_limits
            for j in range(initial_clip_iterations):

                difference = (np.polyval(continuum_coefficients, data.disp) \
                    - data.flux)[mask]
                sigma = np.std(difference)
                exclude = (difference > lower_clip * sigma) \
                    + (difference < -upper_clip * sigma)

                # Remove the excluded points.
                continuum_coefficients = np.polyfit(data.disp[mask][~exclude],
                    data.flux[mask][~exclude], num_coefficients - 1)

            # Save the continuum coefficients.
            initial_theta.update(dict(zip(
                ["c_{0}_{1}".format(self.channel_index, j) for j in range(num_coefficients)],
                continuum_coefficients)))
            continuum = np.polyval(continuum_coefficients, self.data.disp)

        else:
            continuum = 1.

        # Estimate line absorption parameters.
        if absorption_profile_parameters:
            for i, (wavelength, species) in enumerate(config["balance"]["atomic_lines"]):
                if "ld_{0}".format(i) in self.parameters:
                    profile = AbsorptionProfile(
                        self.parameters.get("wl_{0}".format(i), wavelength),
                        **self.absorption_profile_kwargs)
                    for parameter, value in profile.initial_guess(data).iteritems():
                        if parameter not in ("Po", "Vo", "Yo"):
                            initial_theta["{0}_{1}".format(parameter, i)] = value

        # Some final defaults.
        initial_theta.update({
            "Po": 0.1,
            "Vo": 0.01,
            "Yo": np.median(continuum),
            "z": 0. # [TODO]
        })
        return dict(zip(self.parameters, [initial_theta.get(p, np.nan) \
            for p in self.parameters]))


    def optimise(self, data, maxfun=10e4, maxiter=10e4, xtol=1e-8, ftol=1e-8,
        threads=24, full_output=False):
        """
        Optimise the model parameters given the data.

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :param continuum:
            The continuum shape to apply.

        :type continuum:
            float or :class:`numpy.array` with the same length as ``data``

        :param maxfun: [optional]
            The maximum number of function evaluations to perform.

        :type maxfun:
            int

        :param maxiter: [optional]
            The maximum number of function iterations to perform.

        :type maxiter:
            int

        :param xtol: [optional]
            The tolerance required in the model parameters.

        :type xtol:
            float

        :param ftol: [optional]
            The tolerance required in logarithmic probability.

        :type xtol:
            float

        :param threads: [optional]
            The maximum number of parallel threads to use.

        :type threads:
            int

        :param full_output: [optional]
            Return a tuple containing the optimised points theta, the logarithmic
            probability of the optimised point, the number of function iterations,
            the number of function evaluations, and an integer warning flag.

        :type full_output:
            bool

        :returns:
            The optimised model parameters.
        """

        initial_theta = self.initial_guess(data, absorption_profile_parameters=False)

        # Create the continuum.
        num_coefficients = self._get_num_coefficients()
        if num_coefficients > 0:
            coefficients = [initial_theta["c_{0}_{1}".format(self.channel_index, k)] \
                for k in range(num_coefficients)]
            continuum = np.polyval(coefficients, self.data.disp)

        else:
            continuum = 1.

        xopts = []
        processes = []
        pool = mp.Pool(threads)
        for i, (wavelength, species) in enumerate(self._configuration["balance"]["atomic_lines"]):
            if "ld_{0}".format(i) in self.parameters:
                
                # Optimise the parameters of the absorption profile.
                profile = AbsorptionProfile(wavelength, **self.absorption_profile_kwargs)
                process = pool.apply_async(_optimiser, args=(profile, data, continuum))
                processes.append([i, profile, process])
            
        for i, profile, process in processes:
            xopt = dict(zip(profile.parameters, process.get()))
            for parameter in set(xopt).difference(["Po", "Vo", "Yo"]):
                initial_theta["{0}_{1}".format(parameter, i)] = xopt[parameter]
            
        pool.close()
        pool.join()


        posteriors, sampler, info = self.infer(data, np.array([initial_theta[p] for p in self.parameters]))


        import matplotlib.pyplot as plt

        raise a

        # OK, now optimise (minimise) the negative log probability.
        op_kwargs = {
            "maxfun": maxfun,
            "maxiter": maxiter,
            "xtol": xtol,
            "ftol": ftol,
            "disp": False,
            "full_output": True # Necessary for introspection and provenance.
        }
        t_init = time()
        op_theta, op_fopt, op_niter, op_nfunc, op_warnflag = op.fmin(
            lambda theta, data: -self._log_probability(theta, data),
            [initial_theta.get(parameter) for parameter in self.parameters],
            args=(data, ), **op_kwargs)
        self._opt_warn_message(op_warnflag, op_niter, op_nfunc)


        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        ax.plot(data.disp, data.flux, 'k')
        model_spectra = self(dispersion=data.disp, **initial_theta)
        ax.plot(model_spectra[:,0], model_spectra[:,1], 'b')

        model_spectra = self(dispersion=data.disp, **dict(zip(self.parameters, op_theta)))
        ax.plot(model_spectra[:,0], model_spectra[:,1], 'r')

        raise a


        t_elapsed = time() - t_init
        logger.info("Optimisation took {0:.2f} seconds".format(t_elapsed))

        if full_output:
            return (op_theta, op_fopt, op_niter, op_funcalls, op_warnflag)
        return op_theta        


    def infer(self, data, opt_theta, walkers=200, burn_in=200, sample=100,
        **kwargs):
        """
        Infer the model parameters given the data.

        :param opt_theta:
            The initial starting point of theta.

        :type opt_theta:
            :class:`numpy.ndarray`

        :param walkers:
            The number of Goodman & Weare (2010) ensemble walkers.

        :type walkers:
            int

        :param burn_in:
            The number of MCMC steps to discard as burn-in.

        :type burn_in:
            int

        :param sample:
            The number of MCMC steps to sample the posterior with.
        
        :type sample:
            int

        :returns:
            The posterior quantiles in all parameters, the model sampler, 
            and a ``dict`` containing the mean acceptance fractions, concatenated
            chains and log-probability values for the burn-in and posterior.
        """

        #stds = [std.get(parameter, 1e-8) for parameter in self.parameters]
        p0 = [opt_theta + 1e-4*np.random.randn(len(self.parameters)) for i in range(walkers)]
        sampler_kwargs = {
            "threads": self._configuration["settings"].get("max_sampler_threads", 1) \
                if "settings" in self._configuration else 1
        }
        sampler_kwargs.update(kwargs)
        sampler = emcee.EnsembleSampler(walkers, len(self.parameters),
            _inferencer, args=(self, data, ), **sampler_kwargs)

        mean_acceptance_fractions = np.zeros((burn_in + sample))

        for i, (pos, lnprob, rstate) in enumerate(sampler.sample(p0, \
            iterations=burn_in)):
            mean_acceptance_fractions[i] = np.mean(sampler.acceptance_fraction)
            print("Sampler has finished step {0:.0f} of burn-in with "\
                "<a_f> = {1:.3f}, maximum log probability in last step was "\
                "{2:.3e}".format(i + 1, mean_acceptance_fractions[i],
                    np.max(sampler.lnprobability[:, i])))

            if mean_acceptance_fractions[i] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i]))
            
        # Save the chain and calculated log probabilities for later
        chain, lnprobability = sampler.chain, sampler.lnprobability
        sampler.reset()

        logger.info("Sampling posterior...")

        for j, state in enumerate(sampler.sample(pos, iterations=sample)):

            mean_acceptance_fractions[i + j + 1] = np.mean(sampler.acceptance_fraction)
            print("Sampler has finished step {0:.0f} of sampling with "\
                "<a_f> = {1:.3f}, maximum log probability in last step was "\
                "{2:.3e}".format(j + 1, mean_acceptance_fractions[i + j + 1],
                    np.max(sampler.lnprobability[:, j])))

            if mean_acceptance_fractions[i + j + 1] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i + j + 1]))

        # Concatenate the existing chain and lnprobability with the posteriors
        chain = np.concatenate([chain, sampler.chain], axis=1)
        lnprobability = np.concatenate([lnprobability, sampler.lnprobability], axis=1)

        # Get the quantiles.
        posteriors = dict(zip(self.parameters, map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
            zip(*np.percentile(sampler.chain.reshape(-1, len(self.parameters)), [16, 50, 84], axis=0)))))

        info = {
            "chain": chain,
            "lnprobability": lnprobability,
            "mean_acceptance_fractions": mean_acceptance_fractions
        }

        return (posteriors, sampler, info)



class StellarSpectrum(Model):

    def __init__(self, configuration, data):
        """
        A class to model stellar spectra with continuum and absorption profiles.
        The measured absorption profiles can then be used to infer astrophysical
        quantities like stellar parameters and abundances.

        This model uses :class:`SpectralChannel` objects to model each observed 
        spectral channel.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, or a string containing the configuration. Only YAML-style
            formats are accepted.

        :type configuration: str
        """
        super(StellarSpectrum, self).__init__(configuration)

        self.data = data
        # Create channels.
        self.channels = [SpectralChannel(configuration, i) for i in range(len(data))]
        return None


    @property
    def parameters(self):
        """ Return the parameters of the model. """
        parameters = []
        for channel in self.channels:
            parameters.extend(channel.parameters)
        return list(set(parameters))


    def __call__(self, return_continuum=False, **theta):
        """
        Generate data for the given :math:`\\theta`.

        :param return_continuum: [optional]
            Specify whether to return the continuum in addition to the expected
            spectra.

        :type return_continuum:
            bool

        :keyword theta:
            A dictionary mapping of the model parameters (keys) and values.

        :type theta:
            dict
        """

        return [channel(**theta) for channel in self.channels]


    def initial_guess(self, **kwargs):
        return [channel.initial_guess(**kwargs) for channel in self.channels]



    def optimise(self, data, initial_theta=None, max_global_iterations=1,
        threads=-1, op_line_kwargs=None):
        """
        Optimise the model parameters given the data.

        """

        for channel, spectrum in zip(self.channels, data):
            result = channel.optimise(spectrum)

        if initial_theta is None:
            initial_theta = self.initial_guess(data)

        line_kwargs = { "xtol": 1e-8, "ftol": 1e-8, "maxfun": 10e4, "maxiter": 10e4,
            "full_output": True, "disp": False }
        channel_kwargs = { "ftol": 1e-3, "maxfun": 10e4, "maxiter": 10e4,
            "full_output": True, "disp": False } 

        if op_line_kwargs is not None:
            line_kwargs.update(op_line_kwargs)
       
        opt_theta = {}
        opt_theta.update(initial_theta)
        z_parameter = "z" if "z" in initial_theta else "z_{0}"
        model_redshift = z_parameter in self.parameters
        threads = mp.cpu_count() if threads < 0 else threads
        
        # Create SpectralChannel models for each observed spectrum.
        # Then we will optimise each model.
        
        opt_theta = {}
        for i in range(max_global_iterations):
            for j, channel in enumerate(self.channels):

                result = channel.optimise()
                
                # Change the channel continuum coefficient parameter names

                # Change the redshift parameter name if necessary



        config = self._configuration["balance"]
        
        for j in range(max_global_iterations):

            # Fit all the lines in each spectrum.
            num_total_coefficients = 0
            for i, spectrum in enumerate(data):
                num_coefficients, z = 0, opt_theta.get(z_parameter.format(i), 0)
                mask = self.mask(spectrum.disp, z)
                while "c_{0}_{1}".format(i, num_coefficients) in self.parameters:
                    num_coefficients += 1
                num_total_coefficients += num_coefficients

                # Create a continuum.
                continuum = mask.copy()
                if num_coefficients > 0:
                    coefficients = [opt_theta["c_{0}_{1}".format(i, k)] \
                        for k in range(num_coefficients)]
                    continuum *= np.polyval(coefficients, spectrum.disp)

                xopts = []
                processes = []
                pool = mp.Pool(threads)
                for k, (wavelength, species) in enumerate(config["atomic_lines"]):

                    # Should we even be fitting this line?
                    if not (spectrum.disp[-1] >= wavelength >= spectrum.disp[0]):
                        continue

                    # Fit the absorption line.
                    p0 = [initial_theta["{0}_{1}".format(each, k)] \
                        for each in ("ld", "sigma", "wl")]

                    processes.append(pool.apply_async(_absorption_line_fitter,
                        args=(p0, spectrum, continuum, wavelength, line_kwargs)))
                    
                for k, process in enumerate(processes):
                    xopt = process.get()
                    opt_theta.update({
                        "ld_{0}".format(k): np.clip(xopt[0], 0, 1),
                        "sigma_{0}".format(k): xopt[1],
                    })
                    if "wl_{0}".format(k) in opt_theta:
                        opt_theta["wl_{0}".format(k)] = xopt[2]

                pool.close()
                pool.join()

          

                # Do we need to optimise global parameters now?
                if num_total_coefficients > 0 or (model_redshift and z_parameter == "z"):

                    # We *just* want these channel parameters.
                    single_channel_data = [None] * len(data)
                    single_channel_data[i] = spectrum



                    p0 = np.array([opt_theta[p] for p in self.parameters])
                    xopt, fopt, niter, nfunc, warnflag = op.fmin(
                        lambda t, *args: -inference.log_probability(t, *args), p0,
                        args=(self, single_channel_data))
                    self._opt_warn_message(warnflag, niter, nfunc)

                    # Update the dictionary of optimised values.
                    opt_theta.update(dict(zip(self.parameters, xopt)))

                elif (z_parameter == "z" and not model_redshift):
                    # No point re-doing the fitting.
                    break


        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(4)
        best_models = self(data=data, **opt_theta)

        for ax, spectrum, model in zip(axes, data, best_models):
            ax.plot(spectrum.disp, spectrum.flux, 'k')
            ax.plot(model_spectra[:,0], model_spectra[:,1], 'b')


        raise a


        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(spectrum.disp, spectrum.flux, 'k')

        p1 = [opt_theta.get(p) for p in self.parameters]
        foobar = self(p1, data=data)

        raise a
        # Should we be doing a global optimisation (all channels simultaneously)
        # if we have a single redshift?

        return opt_theta
        


class GenerativeModel(Model):

    def __init__(self, configuration):
        """
        A class to probabilistically model stellar spectra. This class performs
        on-the-fly synthesis to generate data for some given theta.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, or a string containing the configuration. Only YAML-style
            formats are accepted.

        :type configuration: str
        """
        super(GenerativeModel, self).__init__(configuration)
        return None


    @property
    def parameters(self):
        """ Return the model parameters. """

        # [TODO] Check for 1D/<3D> model atmospheres to determine free stellar
        #        parameters.
        config = self._configuration
        parameters = ["teff", "logg", "log_Fe", "xi"]

        # Any redshift parameters?
        if config["model"]["redshift"]:
            # [TODO] Allow individual channel redshifts
            parameters.append("z")

        # Any continuum parameters?
        if config["model"]["continuum"]:
            for i, each in enumerate(config["continuum"]["order"]):
                parameters.extend(["c_{0}_{1}".format(i, j) for j in range(each)])

        # Any doppler broadening?
        if config["model"]["doppler_broadening"]:
            # [TODO] You idiot.
            parameters.extend(["doppler_sigma_{0}".format(i) for i in range(1)])

        # Any outlier parameters?
        if config["model"]["outliers"]:
            parameters.extend(["Po", "Vo", "Ys"])

        # Any element parameters (other than Fe)?
        #if "elements" in config["model"]:
        #    parameters.extend(["log_{0}".format(each) \
        #        for each in set(config["model"]["elements"]).difference("Fe")])

        # [TODO] Any telluric parameters?

        return parameters


    def _log_prior(self, theta):
        """
        Return the logarithmic prior probability of the parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :returns:
            The logarithmic prior probability given the parameters ``theta``.

        :rtype:
            float
        """

        theta_dict = dict(zip(self.parameters, theta))

        if not (1 > theta_dict.get("Po", 0.5) > 0) \
        or 0 > theta_dict.get("Vo", 1) \
        or 0 > theta_dict.get("xi", 1) \
        or 0 > theta_dict.get("teff"):
            return -np.inf
        
        ln_prior = 0
        for parameter, distribution in self._configuration.get("priors", {}).iteritems():
            f = eval(distribution, _eval_environment_)
            ln_prior += f(theta_dict[parameter])
        return ln_prior


    def _log_likelihood(self, theta, data):
        """
        Return the logarithmic likelihood of the parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :returns:
            The logarithmic likelihood of the parameters ``theta``.

        :rtype:
            float
        """

        theta_dict = dict(zip(self.parameters, theta))
        expected = self(dispersions=[spectrum.disp for spectrum in data],
            **theta_dict)

        z = theta_dict.get("z", 0.)
        likelihood = 0
        for observed, model in zip(data, expected):

            mask = self.mask(observed.disp, z)
            chi_sq = (observed.flux - model[:, 1])**2 * observed.ivariance * mask

            if "Po" not in theta_dict:
                finite = np.isfinite(chi_sq)
                likelihood += -0.5 * np.sum(chi_sq[finite])

            else:
                Po, Vo, Yo = [theta_dict.get(each) for each in ("Po", "Vo", "Ys")]
                model_likelihood = -0.5 * (chi_sq - np.log(observed.ivariance))
                outlier_ivariance = 1.0/(Vo + observed.variance)
                outlier_likelihood = -0.5 * ((observed.flux - Yo)**2 \
                    * outlier_ivariance * mask - np.log(outlier_ivariance))
                mixture_likelihood = np.logaddexp(
                    np.log(1. - Po) + model_likelihood,
                    np.log(Po)      + outlier_likelihood)
                finite = np.isfinite(mixture_likelihood)
                likelihood += np.sum(mixture_likelihood[finite])

        return likelihood


    def _log_probability(self, theta, data):
        """
        Return the logarithmic probability of the model parameters ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :returns:
            The logarithmic probability of the parameters ``theta``.

        :rtype:
            float
        """
        
        log_prior = self._log_prior(theta)
        if not np.isfinite(log_prior):
            return log_prior
        log_likelihood = self._log_likelihood(theta, data)
        log_probability = log_prior + log_likelihood
        logger.debug("Calculated logarithmic prior, likelihood, and probability"\
            " for {0} to be {1:.3e}, {2:.3e}, and {3:.3e}".format(
            ", ".join(["{0} = {1:.3f}".format(p, v) \
                for p, v in zip(self.parameters, theta)]),
            log_prior, log_likelihood, log_probability))
        return log_probability


    def initial_guess(self, data, synth_kwargs=None):
        """
        Make an initial guess of the model parameters for the given data. The
        initial guess is based on information available in this order: (1) if
        there is an "initial_guess" section in the model configuration file,
        (2) any prior distributions specified in the model configuration file,
        (3) some (sensible) default rules and (4) logic.

        :param data:
            A list of the observed spectra.

        :type data:
            list of :class:`specutils.Spectrum1D' objects

        :param synth_kwargs: [optional]
            Keyword arguments to pass directly to :func:`si.synthesise`

        :type synth_kwargs:
            dict
        """

        if synth_kwargs is None:
            synth_kwargs = {}

        synth_kwargs.setdefault("threads", self._configuration["settings"].get(
            "max_synth_threads", 1) if "settings" in self._configuration else 1)

        parameter_guess = {}
        environment = dict(zip(["locals", "globals", "__name__", "__file__",
            "__builtins__"], [None] * 5))

        # Add our distributions and other explicit functions we want to use.
        environment.update({
            "abs": abs,
            "normal": np.random.normal,
            "uniform": np.random.uniform
        })

        # Any explicit initial guess information?
        if "initial_guess" in self._configuration:

            sampled_filenames = {}
            for parameter, rule in self._configuration["initial_guess"].iteritems():
                if parameter not in model.parameters: continue

                # If we are drawing from a filename, then we need each parameter
                # to be drawn from the same line in the filename, not as separate
                # draws.
                if "sample_from" in rule:
                    # See if we have already sampled from this filename.
                    filename, column = eval(rule, {"sample_from": lambda f, c: (f, c)})
                    if filename not in sampled_filenames:
                        content = np.loadtxt(filename, dtype=str)
                        sampled_filenames[filename] = content[np.random.randint(
                            0, content.shape[0], 1)]

                    parameter_guess[parameter] = float(sampled_filenames[filename][column])

                else:
                    # The rule could be an expression (e.g., normal, uniform),
                    # a value, or some logic expression. So we will evaluate it
                    # in a "safe" environment.
                    parameter_guess[parameter] = eval(rule, environment)

        # Any explicit prior distributions set?
        if "priors" in self._configuration:
            for parameter, rule in self._configuration["priors"].iteritems():
                if parameter in self.parameters \
                and parameter not in parameter_guess:
                    parameter_guess[parameter] = eval(rule, environment)

        # Apply some default rules in case we haven't guessed these parameters
        # yet.
        default_rules = collections.OrderedDict([
            ("teff", np.random.uniform(3000, 7000)),
            ("logg", np.random.uniform(0, 5)),
            ("log_Fe", np.random.uniform(-2, 0)),
            ("xi", "(1.28 + 3.3e-4 * (teff - 6000) - 0.64 * (logg - 4.5)) "\
                "if logg > 3.5 else 2.70 - 0.509 * logg"),
            ("Po", np.random.uniform(0, 1)),
            ("Vo", abs(np.random.normal(0, 1))),
            ("Ys", np.random.normal(1, 0.01))
        ])
        # [TODO] Vo and Ys above. What are you doing?!

        for parameter, rule in default_rules.iteritems():
            if parameter not in self.parameters \
            or parameter in parameter_guess: continue

            # Here we use a local environment so we can pass parameter guesses
            # so far into it.
            local_environment = environment.copy()
            local_environment.update(parameter_guess)
            parameter_guess[parameter] = rule if isinstance(rule, (int, float)) \
                else eval(rule, local_environment)      

        # OK, now just use logic for the remaining parameters.
        # [TODO] Set the abundances of all other elements to be zero:
        remaining_parameters = set(self.parameters).difference(parameter_guess.keys())

        # OK, what about parameter stubs
        default_rule_stubs = collections.OrderedDict([
            ("doppler_sigma_", "uniform(0, 0.2)")
        ])

        for parameter in remaining_parameters:
            for parameter_stub, rule in default_rule_stubs.iteritems():
                if parameter.startswith(parameter_stub):
                    local_environment = environment.copy()
                    local_environment.update(parameter_guess)

                    parameter_guess[parameter] = rule if isinstance(rule, (int, float)) \
                        else eval(rule, local_environment)

        z_parameter = "z" if "z" in remaining_parameters else \
            ("z_{0}" if "z_0" in remaining_parameters else False)
        any_continuum = len(self._configuration["continuum"].get("order", [])) > 0 \
            if hasattr(self._configuration, "continuum") else False

        # Remaining parameters could be related to:
        # tellurics, continuum, redshift.
        if z_parameter or any_continuum:
            
            # Synthesis may be required.
            for i, observed_spectrum in enumerate(data):

                # Synthesis may be required if:
                # (1) this channel requires continuum coefficients, or
                # (2) there is an individual redshift for this channel, or
                # (3) there is a global redshift and this is the first channel.
                synthesis_required = ("c_{0}_0".format(i) in self.parameters) \
                    or (z_parameter == "z_{0}") or (z_parameter == "z" and i == 0)
                if not synthesis_required: continue

                model_spectrum = si.synthesise(
                    parameter_guess["teff"],
                    parameter_guess["logg"],
                    parameter_guess["log_Fe"],
                    parameter_guess["xi"],
                    (observed_spectrum.disp[0], observed_spectrum.disp[-1]),
                    **synth_kwargs)

                # Put the model spectrum onto the same pixels as the observed
                # spectrum.
                resampled_model_flux = np.interp(observed_spectrum.disp,
                    model_spectrum[:, 0], model_spectrum[:, 1], left=np.nan,
                    right=np.nan)

                # Do we need a redshift estimate?
                if (z_parameter == "z" and i == 0) or z_parameter == "z_{0}":
                    parameter_guess[z_parameter.format(i)] = specutils._cross_correlate(
                        observed_spectrum.disp, observed_spectrum.flux,
                        resampled_model_flux)

                # Do we need continuum coefficient estimates?
                continuum_order = -1
                while "c_{0}_{1}".format(i, 1 + continuum_order) in self.parameters:
                    continuum_order += 1

                if continuum_order >= 0:
                    continuum = observed_spectrum.flux/resampled_model_flux
                    finite = np.isfinite(continuum)
                    continuum_parameters = np.polyfit(
                        observed_spectrum.disp[finite], continuum[finite],
                        continuum_order)

                    # Use .setdefault, not .update, in case the user has
                    # specified initial_guess or some prior distribution on this
                    # parameter
                    for j, continuum_parameter in enumerate(continuum_parameters):
                        parameter_guess.setdefault("c_{0}_{1}".format(i, j), continuum_parameter)

        # Make sure we haven't missed any parameters
        assert len(set(self.parameters).difference(parameter_guess)) == 0
        return parameter_guess


    def optimise(self, data, initial_theta=None, synth_kwargs=None, maxfun=10e2,
        maxiter=10e2, xtol=0.001, ftol=0.001, full_output=False):
        """
        Optimise the logarithmic probability of the model parameters theta given
        the data.

        :param data:
            A list of the observed spectra.

        :type data:
            list of :class:`specutils.Spectrum1D` objects

        :param initial_theta: [optional]
            The initial guess for :math:`\\theta` to optimise from.

        :type initial_theta:
            dict

        :param synth_kwargs: [optional]
            Keyword arguments to pass directly to :func:`si.synthesise`.

        :type synth_kwargs:
            dict

        :param maxfun: [optional]
            The maximum number of function evaluations to perform.

        :type maxfun:
            int

        :param maxiter: [optional]
            The maximum number of function iterations to perform.

        :type maxiter:
            int

        :param xtol: [optional]
            The tolerance required in the model parameters.

        :type xtol:
            float

        :param ftol: [optional]
            The tolerance required in logarithmic probability.

        :type xtol:
            float

        :param full_output: [optional]
            Return a tuple containing the optimised points theta, the logarithmic
            probability of the optimised point, the number of function iterations,
            the number of function evaluations, and an integer warning flag.

        :type full_output:
            bool
        """

        if initial_theta is None:
            initial_theta = self.initial_guess(data, synth_kwargs=synth_kwargs)

        # Check to make sure that the initial_theta contains all of the model
        # parameters.
        missing_parameters = set(self.parameters).difference(initial_theta)
        if len(missing_parameters) > 0:
            raise KeyError("initial guess is missing parameter(s) {0}".format(
                ", ".join(missing_parameters)))

        # OK, now optimise (minimise) the negative log probability.
        op_kwargs = {
            "maxfun": maxfun,
            "maxiter": maxiter,
            "xtol": xtol,
            "ftol": ftol,
            "disp": False,
            "full_output": True # Necessary for introspection and provenance.
        }
        t_init = time()
        op_theta, op_fopt, op_niter, op_nfunc, op_warnflag = op.fmin(
            lambda theta, data: -self._log_probability(theta, data),
            [initial_theta.get(parameter) for parameter in self.parameters],
            args=(data, ), **op_kwargs)
        self._opt_warn_message(op_warnflag, op_niter, op_nfunc)

        t_elapsed = time() - t_init
        logger.info("Optimisation took {0:.2f} seconds".format(t_elapsed))

        op_theta_dict = dict(zip(self.parameters, op_theta))

        raise a

        if full_output:
            return (op_theta_dict, op_fopt, op_niter, op_funcalls, op_warnflag)
        return op_theta_dict


    def __call__(self, dispersions, synth_kwargs=None, **theta):
        """
        Generate data for the given :math:`\\theta`.

        :param dispersions:
            A two-length tuple containing the range of wavelengths to synthesise,
            or a list of two-length tuples containing start and end wavelengths 
            to synthesise.

        :type dispersions:
            list 

        :param synth_kwargs: [optional]
            Keyword arguments to pass directly to :func:`si.synthesise`.

        :type synth_kwargs:
            dict

        :keyword theta:
            A dictionary mapping of the model parameters (keys) and values.

        :type theta:
            dict

        :raises KeyError:
            If a model parameter is missing from ``**theta.``

        :returns:
            A list of two-column arrays that represent dispersion and flux for
            each observed channel in ``data``, or each specified region in
            ``wavelength_ranges``. 

        :rtype:
            list
        """
        
        # [TODO] Redshift the *requested* synthesis range based on z in theta?

        missing_parameters = set(self.parameters).difference(theta)
        if len(missing_parameters) > 0:
            logging.warn("Generate call is missing some theta parameters: {0}".format(
                ", ".join(missing_parameters)))

        # Get the number of threads we should use.
        threads = self._configuration["settings"].get("max_synth_threads", 1) \
            if "settings" in self._configuration else 1

        if synth_kwargs is None:
            synth_kwargs = {}

        # Request the wavelength step to be ~twice the observed pixel sampling.
        #if data is not None:
        #    synth_kwargs.setdefault("wavelength_steps",
        #        [(0., np.min(np.diff(each.disp))/2., 0.) for each in data])
            
        synthesis_ranges = [(each[0], each[-1]) for each in dispersions]
    
        # Synthesise the model spectra first (in parallel where applicable) and
        # then apply the cheap transformations in serial.
        model_spectra = si.synthesise(theta["teff"], theta["logg"], theta["log_Fe"],
            theta["xi"], synthesis_ranges, threads=threads, **synth_kwargs)

        continuum_spectra = []
        for i in range(len(model_spectra)):

            model_spectrum = model_spectra[i]
            
            # Apply macroturbulent, vsini and instrumental broadening to the
            # spectrum.
            # [TODO] For the moment we only have instrumental broadening.
            kernel = theta.get("doppler_sigma_{0}".format(i), 0.)
            if kernel > 0:
                pixel_kernel = (kernel/(2 * (2*np.log(2))**0.5))/np.mean(np.diff(model_spectrum[:, 0]))
                model_spectrum[:, 1] = gaussian_filter1d(model_spectrum[:, 1],
                    pixel_kernel)

            # Apply redshift corrections to the model spectra.
            z = theta["z"] if "z" in theta else theta.get("z_{0}".format(i), 0)
            model_spectrum[:, 0] *= (1. + z)

            # Interpolate the new model spectra onto the observed pixels.
            # We do this before transforming the continuum such that the
            # continuum spectra (if it is returned with return_continuum)
            # will always be on the same dispersion map as the returned
            # spectra.
            if len(dispersions[i]) > 2:
                expected = np.interp(dispersions[i], model_spectrum[:, 0],
                    model_spectrum[:, 1], left=np.nan, right=np.nan)
                model_spectrum = np.vstack([dispersions[i], expected]).T

            # Apply continuum transformation to the model spectra.
            j, coefficients = 0, []
            while "c_{0}_{1}".format(i, j) in theta:
                coefficients.append(theta["c_{0}_{1}".format(i, j)])
                j += 1

            if len(coefficients) > 0:
                continuum_spectra.append(np.polyval(coefficients, model_spectrum[:, 0]))
                model_spectrum[:, 1] *= continuum_spectra[-1]

            else:
                # No continuum for this order. Fill it with a None to maintain
                # order.
                continuum_spectra.append(None)

            model_spectra[i] = model_spectrum

        #if return_continuum:
        #    return (model_spectra, continuum_spectra)
        return model_spectra


    def infer(self, data, p0=None, walkers=200, burn_in=100, sample=100,
        **kwargs):
        """
        Infer the model parameters given the data.

        :param data:
            A list of the observed spectra.

        :type data:
            list of :class:`specutils.Spectrum1D' objects

        :param p0:
            The initial starting point of shape [TODO].

        :type p0:
            :class:`numpy.ndarray`

        :param walkers:
            The number of Goodman & Weare (2010) ensemble walkers.

        :type walkers:
            int

        :param burn_in:
            The number of MCMC steps to discard as burn-in.

        :type burn_in:
            int

        :param sample:
            The number of MCMC steps to sample the posterior with.
        
        :type sample:
            int

        :returns:
            The posterior quantiles in all parameters, the model sampler, 
            and a ``dict`` containing the mean acceptance fractions, concatenated
            chains and log-probability values for the burn-in and posterior.
        """

        sampler_kwargs = {
            "threads": self._configuration["settings"].get("max_sampler_threads", 1) \
                if "settings" in self._configuration else 1
        }
        sampler_kwargs.update(kwargs)
        sampler = emcee.EnsembleSampler(walkers, len(self.parameters),
            inference.log_probability, args=(self, data), **sampler_kwargs)

        mean_acceptance_fractions = np.zeros((burn_in + sample))

        for i, (pos, lnprob, rstate) in enumerate(sampler.sample(p0, \
            iterations=burn_in)):
            mean_acceptance_fractions[i] = np.mean(sampler.acceptance_fraction)
            logger.info("Sampler has finished step {0:.0f} of burn-in with "\
                "<a_f> = {1:.3f}, maximum log probability in last step was "\
                "{2:.3e}".format(i + 1, mean_acceptance_fractions[i],
                    np.max(sampler.lnprobability[:, i])))

            if mean_acceptance_fractions[i] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i]))

        # Save the chain and calculated log probabilities for later
        chain, lnprobability = sampler.chain, sampler.lnprobability
        sampler.reset()

        logger.info("Sampling posterior...")

        for j, state in enumerate(sampler.sample(pos, iterations=sample)):

            mean_acceptance_fractions[i + j + 1] = np.mean(sampler.acceptance_fraction)
            logger.info("Sampler has finished step {0:.0f} of sampling with "\
                "<a_f> = {1:.3f}, maximum log probability in last step was "\
                "{2:.3e}".format(j + 1, mean_acceptance_fractions[i + j + 1],
                    np.max(sampler.lnprobability[:, i + j + 1])))

            if mean_acceptance_fractions[i + j + 1] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i + j + 1]))

        # Concatenate the existing chain and lnprobability with the posteriors
        chain = np.concatenate([chain, sampler.chain], axis=1)
        lnprobability = np.concatenate([lnprobability, sampler.lnprobability], axis=1)

        # Get the quantiles.
        posteriors = dict(zip(self.parameters, map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
            zip(*np.percentile(sampler.chain.reshape(-1, len(self.parameters)), [16, 50, 84], axis=0)))))

        info = {
            "chain": chain,
            "lnprobability": lnprobability,
            "mean_acceptance_fractions": mean_acceptance_fractions
        }

        return (posteriors, sampler, info)