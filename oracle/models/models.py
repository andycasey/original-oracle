# coding: utf-8

""" Probabilistic Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import itertools
import logging
import multiprocessing as mp
import os
import warnings
import yaml
from functools import partial
from time import time

import numpy as np
from scipy import optimize as op, integrate, stats

from scipy.ndimage import gaussian_filter1d
import emcee

import plot
import moog
import si
import specutils

import line
import profiles
import utils
from profiles import AbsorptionProfile

__all__ = ["Star", "StellarSpectrum", "GenerativeModel"]

logger = logging.getLogger("oracle")

# This dictionary will be used to evaluate prior distributions provided in the
# model configuration file. This needs to be pickleable and accessible to all
# functions, which is why it's out here.
_log_prior_eval_environment = { 
    "locals": None,
    "globals": None,
    "__name__": None,
    "__file__": None,
    "__builtins__": None,
    "uniform": lambda a, b: partial(stats.uniform.logpdf, **{"loc": a, "scale": b - a}),
    "normal": lambda a, b: partial(stats.norm.logpdf, **{"loc": a, "scale": b})
}


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
            self.config = {}
            self.config.update(configuration)

        else:
            if not os.path.exists(configuration):
                # Probably a string configuration.
                try:
                    self.config = yaml.load(configuration)

                except:
                    raise IOError("configuration file does not exist or the"\
                        " YAML string provided does not describe a valid "\
                        "dictionary")
                else:
                    # We expect a dictionary.
                    if not isinstance(self.config, dict):
                        raise IOError("configuration file does not exist or the"\
                            " YAML string provided does not describe a valid "\
                            "dictionary")
            else:
                with open(configuration, "r") as fp:
                    self.config = yaml.load(fp)

        # Check the configuration is valid.
        self._validate()
        return None


    def _validate(self):
        """ Check that the configuration is valid. """

        if not self.config.has_key("model"):
            raise KeyError("no model information specified")

        validation_functions = {
            "continuum": self._validate_continuum,
            "elements": self._validate_elements
        }
        for item, state in self.config["model"].iteritems():
            if not state or item not in validation_functions: continue
            validation_functions[item]()

        return True


    def _validate_continuum(self):
        """ Check that the continuum configuration is valid. """
        
        # We actually need *something* specified to model the continuum.
        if not self.config.has_key("continuum"):
            raise KeyError("no information specified for continuum modelling")

        order = self.config["continuum"]["order"]
        assert self.config["continuum"]["method"] == "polynomial"
        assert order

        # Order should be integer or list-like of integers.
        try: order = [int(order)]
        except:
            try: order = map(int, order)
            except:
                raise TypeError("continuum order must be an integer or a list-\
                    like of integers")
        self.config["continuum"]["order"] = order
        
        return True


    def _validate_elements(self):
        """ Check that the elements actually exist. """

        map(utils.atomic_number, self.config["model"]["elements"])
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


    def mask(self, dispersion, z=0., fill_value=np.nan):
        """
        Return an array mask for a given dispersion array and redshift, based on
        the mask information provided in the model configuration file.

        :param dispersion:
            An array of dispersion values to mask.

        :type dispersion:
            :class:`numpy.array`

        :param z: [optional]
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

        if z == 0:
            try:
                return self._rest_mask

            except AttributeError:
                None

        mask_regions = self.config.get("mask", None)
        if mask_regions is not None:
            mask = np.ones(len(dispersion))
            mask_regions = np.array(mask_regions)
            for start, end in mask_regions * (1. + z):
                indices = dispersion.searchsorted([start, end])
                mask[indices[0]:indices[1] + 1] = fill_value

            if z == 0:
                self._rest_mask = mask

            return mask

        else:
            return 1.


def _optimiser(_class, *args):
    return _class.optimise(*args)


def _inferencer(theta, _class, *args):
    return _class._log_probability(theta, *args)
        

class SpectralChannel(Model):

    wl_contribution = 2.
    
    def __init__(self, configuration, channel_index):
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

        """

        super(SpectralChannel, self).__init__(configuration)
        self._data = None
        self.channel_index = channel_index
        return None


    @property
    def parameters(self, data=None):
        """ Return the model parameters. """
        
        if hasattr(self, "_parameters"):
            return self._parameters

        parameters = []
        config = self.config
        
        # Any redshift parameters?
        if config["model"]["redshift"]:
            parameters.append("z")

        # Any continuum parameters?
        if config["model"]["continuum"]:
            order = config["continuum"]["order"]
            order = order if isinstance(order, int) \
                else order[self.channel_index]
            parameters.extend(["c_{0}_{0}".format(self.channel_index, j) \
                for j in range(order)])

        # Any outlier parameters?
        if config["model"]["outliers"]:
            parameters.extend(["Po", "Vo", "Yo"])

        # Underestimated variance?
        if config["model"]["noise_model"]:
            parameters.append("f")

        # We can't get the exact model parameters unless we know the observed
        # dispersion range.
        transition_parameters = ["shape", "ld"]
        if self.config["classical"].get("use_instrumental_profile", False):
            parameters.append("fwhm_instrumental")
        else:
            transition_parameters.append("fwhm")

        for i, transition_data in enumerate(config["classical"]["atomic_lines"]):
            wavelength = transition_data[0]
            if self._data is None \
            or (self._data.disp[-1] >= wavelength >= self._data.disp[0]):
                parameters.extend(["{0}_{1}".format(each, i) \
                    for each in transition_parameters])

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

        if isinstance(dispersion, (list, tuple)):
            dispersion = dispersion[0]
        
        num_coefficients = self._get_num_coefficients()
        
        # Create continuum (but we'll call it flux)
        if num_coefficients > 0:
            coefficients = [theta["c_{0}_{1}".format(self.channel_index, i)] \
                for i in range(num_coefficients)]
            flux = np.polyval(coefficients, dispersion)

        else:
            flux = np.ones(len(dispersion))
        
        # Model the absorption profiles
        fwhm_key = ["fwhm_{0}", "fwhm_instrumental"][self.config["classical"].get("use_instrumental_profile", False)]
        indices = [int(p.split("_")[1]) for p in theta if p.startswith("ld_")]
        for index in indices:
            depth = theta["ld_{0}".format(index)]
            if depth > 0:
                fwhm = theta[fwhm_key.format(index)]
                shape = theta["shape_{0}".format(index)]
                wavelength = self.config["classical"]["atomic_lines"][index][0]

                wl_indices = dispersion.searchsorted([
                    wavelength - self.wl_contribution,
                    wavelength + self.wl_contribution
                ])
                x, y = dispersion.__getslice__(*wl_indices), flux.__getslice__(*wl_indices)
                flux.__setslice__(wl_indices[0], wl_indices[1],
                    y * (1 - depth * profiles.voigt(wavelength, fwhm, shape, x)))

        # [TODO] No redshift modelling.

        # Can probably remove this line..
        result = np.vstack([dispersion, flux]).T
        return result


    def _log_prior(self, theta):
        """
        Return the logarithmic prior probability of the parameters ``theta``
        for the StellarSpectrum model.

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
        or not (10 > theta_dict.get("Vo", 1) > 0):
            return -np.inf

        for i, transition_data in enumerate(self.config["classical"]["atomic_lines"]):
            if not (1 >= theta_dict.get("ld_{0}".format(i), 0) >= -1) \
            or not (1 >= theta_dict.get("fwhm_{0}".format(i), 0.5) > 0) \
            or not (1 >= theta_dict.get("shape_{0}".format(i), 0.5) >= 0):
                return -np.inf
        
        # Priors on fwhm values.
        ln_prior = 0
        if "fwhm_instrumental" not in theta_dict:
            # Scale them to resolution
            line_indices = [int(p.split("_")[1]) for p in self.parameters]
            spectral_resolutions = np.array(
                [self.config["classical"]["atomic_lines"][i][0]/theta_dict["fwhm_{0}".format(i)] \
                    for i in line_indices])
            ln_prior += stats.norm.logpdf(spectral_resolutions,
                loc=np.mean(spectral_resolutions),
                scale=np.std(spectral_resolutions)).sum()

        for parameter, distribution in self.config.get("priors", {}).iteritems():
            f = eval(distribution, _log_prior_eval_environment)
            ln_prior += f(theta_dict[parameter])

        return ln_prior

    
    def _log_likelihood(self, theta, data):
        """
        Return the logarithmic likelihood of the parameters ``theta`` for the
        SpectralChannel model.

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
        t_init = time()

        theta_dict = dict(zip(self.parameters, theta))
        expected = self(dispersion=data.disp, **theta_dict)[:, 1]

        mask = self.mask(data.disp, theta_dict.get("z", 0))

        if "f" in self.parameters:
            ivariance = 1.0/(data.variance + expected**2 * np.exp(2. * theta_dict["f"]))
        else:
            ivariance = data.ivariance

        chi_sq = (data.flux - expected)**2 * ivariance * mask
        stellar_likelihood = -0.5 * (chi_sq - np.log(ivariance))

        if "Po" in theta_dict:
            Po, Vo, Yo = [theta_dict.get(each) for each in ("Po", "Vo", "Yo")]

            outlier_ivariance = 1.0/(Vo + data.variance)
            outlier_likelihood = -0.5 * ((data.flux - Yo)**2 * outlier_ivariance \
                * mask - np.log(outlier_ivariance))
            likelihood = np.logaddexp(
                np.log(1. - Po) + stellar_likelihood,
                np.log(Po)      + outlier_likelihood)

        else:
            likelihood = stellar_likelihood
        
        finite = np.isfinite(likelihood)
        if finite.sum() == 0:
            raise RuntimeError("no observed pixels used to calculate likelihood"\
                " -- are all of your pixels masked?")

        return likelihood[finite].sum()


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
        probability = log_prior + self._log_likelihood(theta, data)
        return probability


    def initial_guess(self, data, initial_clip_iterations=5, resolution=28000,
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
            z = 0. # [TODO]
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
            continuum = np.polyval(continuum_coefficients, data.disp)

        else:
            continuum = 1.

        # Estimate line absorption parameters.
        use_instrumental_profile = self.config["classical"].get("use_instrumental_profile", False)
        if absorption_profile_parameters:
            for i, transition_data in enumerate(self.config["classical"]["atomic_lines"]):
                if "ld_{0}".format(i) in self.parameters:
                    
                    wavelength = transition_data[0]    
                    index = data.disp.searchsorted(wavelength)
                    line_continuum = continuum if isinstance(continuum, (int, float)) \
                        else continuum[index]

                    line_depth = line_continuum - data.flux[index]
                    if not np.isfinite(line_depth):
                        line_depth = 0.
                    initial_theta["ld_{0}".format(i)] = line_depth
                    initial_theta["shape_{0}".format(i)] = 0.
                    if not use_instrumental_profile:
                        initial_theta["fwhm_{0}".format(i)] = wavelength / float(resolution)

        # Some final defaults.
        initial_theta.update({
            "fwhm_instrumental": np.mean(data.disp) / float(resolution),
            "Po": 0.50,
            "Vo": 0.01,
            "Yo": np.median(continuum),
            "f": -2., # log(f)
            "z": 0. # [TODO]
        })

        return dict(zip(self.parameters, [initial_theta.get(p, 0) \
            for p in self.parameters]))


    def optimise(self, data, initial_theta=None, maxfun=10e4, maxiter=10e4,
        xtol=1e-4, ftol=0.1, threads=-1, full_output=False):
        """
        Optimise the SpectralChannel model parameters given the data.

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

        t_init = time()
        if initial_theta is None:
            initial_theta = self.initial_guess(data,
                absorption_profile_parameters=True)

        print("initial theta", initial_theta)

        # Create the continuum.
        num_coefficients = self._get_num_coefficients()
        if num_coefficients > 0:
            coefficients = [initial_theta["c_{0}_{1}".format(self.channel_index, k)] \
                for k in range(num_coefficients)]
            continuum = np.polyval(coefficients, data.disp)

        else:
            continuum = 1.
        
        xopts = []
        processes = []
        threads = mp.cpu_count() if threads < 0 else threads
        pool = mp.Pool(threads)
        mask = self.mask(data.disp)

        pre_opt_theta = {}
        pre_opt_theta.update(initial_theta)

        """
        
        # Find what transitions we have and bunch them by wl_contribution so that
        # any lines within +/- wl_contribution of each other are fit simultaneously
        indices = [int(p.split("_")[1]) for p in self.parameters \
            if p.startswith("ld_")]
        wavelengths = [self.config["classical"]["atomic_lines"][i][0] for i in indices]

        sort_indices = np.argsort(wavelengths)
        indices = np.array(indices)[sort_indices]
        wavelengths = np.array(wavelengths)[sort_indices]
        """

        for i, transition_data in enumerate(self.config["classical"]["atomic_lines"]):
            if "ld_{0}".format(i) in self.parameters:
                
                wavelength = transition_data[0]
                # Optimise the parameters of the absorption profile.

                indices = data.disp.searchsorted([
                    wavelength - self.wl_contribution,
                    wavelength + self.wl_contribution
                ])
                _data = specutils.Spectrum1D(disp=data.disp.__getslice__(*indices),
                    flux=data.flux.__getslice__(*indices).copy(),
                    variance=data.variance.__getslice__(*indices))
                _mask = mask.__getslice__(*indices)
                _data.flux[~np.isfinite(_mask)] = np.nan

                if isinstance(continuum, (int, float)):
                    _continuum = continuum

                else:
                    _continuum = continuum.__getslice__(*indices)

                profile = profiles.AbsorptionProfile(wavelength, fwhm=pre_opt_theta.get("fwhm_instrumental", None))
                process = pool.apply_async(_optimiser,
                    args=(profile, _data, _continuum))
                processes.append([i, profile, process])
        

        ignore_parameters = ["Po", "Vo", "Yo"]
        if self.config["classical"].get("use_instrumental_profile", False):
            ignore_parameters.append("fwhm")
        for i, profile, process in processes:
            #op_theta, op_fopt, op_niter, op_funcalls, op_warnflag = process.get()
            xopt = dict(zip(profile.parameters, process.get()))

            for parameter in set(xopt).difference(ignore_parameters):
                pre_opt_theta["{0}_{1}".format(parameter, i)] = xopt[parameter]
        
        pool.close()
        pool.join()
    
        if True or num_coefficients > 0 \
        or "z" in self.parameters \
        or "Po" in self.parameters:
            op_kwargs = {
                "maxfun": maxfun,
                "maxiter": maxiter,
                "xtol": xtol,
                "ftol": ftol,
                "disp": False,
                "full_output": True # Necessary for introspection and provenance.
            }
            t_init = time()
            op_x, op_fopt, op_niter, op_nfunc, op_warnflag = op.fmin(
                lambda theta, data: -self._log_probability(theta, data),
                [pre_opt_theta.get(parameter) for parameter in self.parameters],
                args=(data, ), **op_kwargs)
            self._opt_warn_message(op_warnflag, op_niter, op_nfunc)
            opt_theta = dict(zip(self.parameters, op_x))


        else:
            op_x = [pre_opt_theta.get(parameter) for parameter in self.parameters]
            opt_theta = pre_opt_theta

        t_elapsed = time() - t_init
        logger.info("Spectral channel optimisation took {0:.2f} seconds".format(t_elapsed))

        assert not np.any(~np.isfinite(np.array(opt_theta.values())))

        if full_output:
            return (op_theta, op_fopt, op_niter, op_funcalls, op_warnflag)
        return opt_theta


    def _initial_proposal(self, opt_theta, walkers):
        """
        Return an initial proposal array for the sampler.

        :param opt_theta:
            The optimised starting point theta.

        :type opt_theta:
            dict

        :param walkers:
            The number of Goodman & Weare (2010) ensemble walkers.

        :type walkers:
            int
        """

        p0 = emcee.utils.sample_ball(
            [opt_theta.get(p) for p in self.parameters], [1e-4] * len(opt_theta),
            size=walkers)

        for i, parameter in enumerate(self.parameters):
            if parameter.startswith("shape_"):
                p0[:, i] = np.random.uniform(0, 1, size=walkers)

        # Do reflections so that we have a reasonable initial starting point.
        limits = {
            "Po": [0, 1],
            "Vo": [0, None],
            "ld": [-1, 1],
            "fwhm": [0, 1]
        }
        for i, parameter in enumerate(self.parameters):
            limit = limits.get(parameter.split("_")[0], False)
            if limit:
                p0[:, i] = utils.reflect_about(p0[:, i], limit)

        return p0


    def infer(self, data, opt_theta, walkers=200, burn=5000, sample=100, threads=1,
        **kwargs):
        """
        Infer the SpectralChannel parameters given the data.

        :param opt_theta:
            The initial starting point of theta.

        :type opt_theta:
            dict

        :param walkers:
            The number of Goodman & Weare (2010) ensemble walkers. If a negative
            value is given then the number of walkers will be set by ``-walkers 
            * len(parameters)``. 

        :type walkers:
            int

        :param burn:
            The number of MCMC steps to discard as burn-in.

        :type burn:
            int

        :param sample:
            The number of MCMC steps to sample the posterior with.
        
        :type sample:
            int

        :param threads: [optional]
            The number of threads to use for the :class:`emcee.EnsembleSampler`.
            Specifying -1 will set ``threads`` to the number of available CPUs.

        :type threads:
            int

        :param kwargs: [optional]
            Keyword arguments to supply directly to :class:`emcee.EnsembleSampler`

        :returns:
            The posterior quantiles in all parameters, the model sampler, 
            and a ``dict`` containing the mean acceptance fractions, concatenated
            chains and log-probability values for the burn-in and posterior.
        """

        assert len(self.parameters) > 0

        walkers = walkers if walkers > 0 else -walkers*len(opt_theta)
        logger.info("Initialising {0} walkers for {1} parameters".format(
            walkers, len(self.parameters)))

        p0 = self._initial_proposal(opt_theta, walkers)

        sampler = emcee.EnsembleSampler(walkers, len(self.parameters),
            _inferencer, args=(self, data), threads=threads, **kwargs)

        t_init = time()
        mean_acceptance_fractions = np.zeros((burn + sample))
        for i, (pos, lnprob, rstate) in enumerate(sampler.sample(p0, \
            iterations=burn)):
            mean_acceptance_fractions[i] = np.mean(sampler.acceptance_fraction)
            logger.info("Sampler has finished step {0:.0f} of burn-in with "\
                "<a_f> = {1:.3f}, maximum log probability in last step was "\
                "{2:.3e}. Mean sample time {3:.1f}s".format(i + 1,
                    mean_acceptance_fractions[i],
                    np.max(sampler.lnprobability[:, i]),
                    (time()-t_init)/float(i + 1)))

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
                    np.max(sampler.lnprobability[:, j])))

            if mean_acceptance_fractions[i + j + 1] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i + j + 1]))

        # Concatenate the existing chain and lnprobability with the posteriors
        chain = np.concatenate([chain, sampler.chain], axis=1)
        lnprobability = np.concatenate([lnprobability, sampler.lnprobability], axis=1)

        # Get the quantiles.
        posteriors = dict(zip(self.parameters, map(lambda v: (v[1], v[2]-v[1], v[0]-v[1]),
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


    def __call__(self, dispersion, theta):
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

        return [c(d, t) for c, d, t in zip(self.channels, dispersion, theta)]


    def initial_guess(self, data, **kwargs):
        return [channel.initial_guess(spectrum, **kwargs) \
            for channel, spectrum in zip(self.channels, data)]



    def optimise(self, data, initial_theta=None,
        threads=-1, op_line_kwargs=None):
        """
        Optimise the model parameters given the data.

        """

        # TODO CHANGE INIITAL_THETA TO SINGLE DICT AND NOT LIST OF DICTS
        #if initial_theta is None:
        #    initial_theta = self.initial_guess(data)

        return [c.optimise(d) for c, d in zip(self.channels, data)]


    def infer(self, data, optimised_theta):
        """
        Infer the model parameters (e.g., continuum, line profile parameters) of
        all the :class:`SpectralChannel` models associated with this :class:`StellarSpectrum`
        class.

        :param data:
            The observed data.

        :type data:
            list of :class:`specutils.Spectrum1D` objects

        :param optimised_theta:
            The optimised values of the model parameters :math:`\\Theta` for each
            of the :class:`SpectralChannel` channels.

        :type optimised_theta:
            list of dict objects
        """

        # [TODO] Parallelise this? One arm per core?
        return [c.infer(s, t) for c, s, t in zip(self.channels, data, optimised_theta)]


    def integrate_profiles(self, channel_thetas, estimated_variance=0.10,
        xlimits=1):
        """
        Integrate measured absorption profiles to find equivalent widths for all
        atomic absorption lines in the current model.

        :param channel_thetas:
            The optimised model parameters for each channel. The list contains
            dictionaries with parameters (as keys) and values for each channel.
            The values can either be a single float (e.g., the optimised value)
            or a two- or three-length tuple first containing the value and then
            the associated positive and negative uncertainties.

        :type channel_thetas:
            list of dicts

        :param estimated_variance: [optional]
            Fraction of variance in measured equivalent widths. This value will
            be used to estimate the variance in equivalent widths if no 
            uncertainties in theta parameters is provided.

        :type estimated_variance:
            float

        :returns:
            An array of atomic transitions and integrated equivalent widths.
        """

        num_atomic_lines = len(self.config["classical"]["atomic_lines"])
        tabular_results = np.zeros((num_atomic_lines, 5))

        # Fill 'er up.
        tabular_results[:, :2] = np.array(self.config["classical"]["atomic_lines"])[:, :2]
        tabular_results[:, 2:] = np.nan

        uip = self.config["classical"].get("use_instrumental_profile", False)
        fwhm_key = ["fwhm_{0}", "fwhm_instrumental"][uip]
        for channel_theta in channel_thetas:
            for parameter, parameter_value in channel_theta.iteritems():
                # Identify a line.
                if parameter.startswith("ld_"):
                    transition_index = int(parameter.split("_")[1])
                    wavelength, atomic_number = tabular_results[transition_index, :2]
                    
                    # Has this transition been measured already in a different
                    # channel?
                    if np.isfinite(tabular_results[transition_index, 2]):
                        logger.warn("{0} transition at {1:.3f} already has a "\
                            "measured equivalent width of {2:.2f} mA".format(
                                utils.element(atomic_number), wavelength,
                                tabular_results[transition_index, -1]))

                    # Do we have a single value, or do we have quantiles?
                    if isinstance(parameter_value, (float, int)):
                        # Integrate the profile.
                        depth = parameter_value
                        shape = channel_theta["shape_{0}".format(transition_index)]
                        fwhm = channel_theta[fwhm_key.format(transition_index)]
                        y = lambda x: depth * profiles.voigt(wavelength, fwhm, shape, x)
                        equivalent_width, integration_error_est = integrate.quad(y,
                            wavelength - xlimits, wavelength + xlimits)

                        # Save it (in mA)
                        equivalent_width *= 1000.
                        tabular_results[transition_index, 2] = equivalent_width

                        # Use an estimated variance if it exists.
                        if estimated_variance > 0:
                            tabular_results[transition_index, 3:] = [
                                +equivalent_width * estimated_variance,
                                -equivalent_width * estimated_variance,
                            ]

                    else:
                        depth = parameter_value[0]
                        shape = channel_theta["shape_{0}".format(transition_index)][0]
                        fwhm = channel_theta[fwhm_key.format(transition_index)][0]
                        y = lambda x: depth * profiles.voigt(wavelength, fwhm, shape, x)
                        equivalent_width, integration_error_est = integrate.quad(y,
                            wavelength - xlimits, wavelength + xlimits)
                        equivalent_width *= 1000.

                        tabular_results[transition_index, 2] = equivalent_width

                        # Create permutations of the profile, given the parameters.
                        depths = np.array(parameter_value)
                        fwhms = np.array(channel_theta[fwhm_key.format(transition_index)])
                        shapes = np.array(channel_theta["shape_{0}".format(transition_index)])
                        permutations = itertools.product(
                            depths[0] + depths[1:],
                            fwhms[0] + fwhms[1:],
                            shapes[0] + shapes[1:]
                        )
                        equivalent_widths = []
                        for depth, fwhm, shape in permutations:
                            y = lambda x: depth * profiles.voigt(wavelength, fwhm, shape, x)
                            integral, integration_error_est = integrate.quad(
                                y, wavelength - xlimits, wavelength + xlimits)
                            equivalent_widths.append(1000. * integral)

                        # Take the minimum and maximum values as projected quantiles.
                        # Restrict equivalent widths to be [0, infinity)
                        tabular_results[transition_index, 3:] = np.array([
                            np.max(equivalent_widths),
                            np.min(equivalent_widths)
                        ]) - equivalent_width
        
        return np.core.records.fromarrays(np.hstack([
                np.array(self.config["classical"]["atomic_lines"]),
                tabular_results[:, 2:].reshape(-1, 3)
            ]).T,
            names=("wavelength", "species", "excitation_potential", "loggf",
                "equivalent_width", "u_pos_equivalent_width", 
                "u_neg_equivalent_width"),
            formats=["f8"] * 7)


    def initial_guess_stellar_parameters(self):
        """
        Provide a completely random guess of the stellar parameters.

        [TODO] Introduce prior information into this function.
        """


        default_rules = collections.OrderedDict([
            ("teff", np.random.uniform(4000, 7000)),
            ("logg", np.random.uniform(0, 5)),
            ("[M/H]", np.random.uniform(-2, 0)),
            ("xi", "(1.28 + 3.3e-4 * (teff - 6000) - 0.64 * (logg - 4.5)) "\
                "if logg > 3.5 else 2.70 - 0.509 * logg"),
            ("Po", np.random.uniform(0, 1)),
            ("Vo", abs(np.random.normal(0, 1))),
            ("Ys", np.random.normal(1, 0.01)),
        ])

        initial_guess = {}
        for parameter, rule in default_rules.iteritems():

            if isinstance(rule, float):
                initial_guess[parameter] = rule

            else:
                initial_guess[parameter] = eval(rule, initial_guess)

        parameter_order = ("teff", "xi", "logg", "[M/H]")

        return [initial_guess.get(p) for p in parameter_order]


    def optimise_stellar_parameters(self, atomic_data, species=(26.0, 26.1), 
        initial_stellar_parameters=None, maxiter=1, ftol=4e-3, full_output=False):
        """
        Optimise stellar parameters (Teff, logg, [Fe/H], xi) given some 
        measured equivalent widths.
        """

        if initial_stellar_parameters is None:
            initial_stellar_parameters = self.initial_guess_stellar_parameters()

        # It *really* bothers me immensely that the singular and plural of the 
        # word 'species' are identical.
        indices = np.zeros(len(atomic_data), dtype=bool)
        for each in species:
            indices[atomic_data["species"] == each] += 1

        if 2 > indices.sum():
            raise ValueError("less than two transitions identified for stellar "\
                "parameter determination")

        atomic_data = atomic_data[indices]
        logger.info("Identified {0} transitions of species {1} to use for stellar"\
            " parameter determination".format(
                len(atomic_data), ", ".join(map(str, species))))

        # Create an instance of MOOGSILENT.
        with moog.instance(debug=False) as moogsilent:

            # Write atomic data to disk
            finite = np.isfinite(atomic_data["equivalent_width"]) \
                * (atomic_data["equivalent_width"] > 0)
            line_list_filename = os.path.join(moogsilent.twd, "lines")
            with open(line_list_filename, "w") as fp:
                fp.write(moogsilent._format_ew_input(atomic_data[finite]))
            
            def excitation_ionisation_balance(theta):
                temperature, xi, logg, metallicity = theta

                try:
                    data = moogsilent.abfind(
                        temperature, logg, metallicity, xi, line_list_filename,
                        clobber=True)
                except:
                    logger.exception("Failed to find abundances for stellar "\
                        "parameters Teff = {0}, logg = {1:.2f}, [M/H] = {2:.2f}"\
                        ", xi = {3:.2f}".format(temperature, logg, metallicity,
                            xi))
                    return np.array([np.inf]*4)

                # Create abundance uncertainties
                y_uncertainty = np.nanmax(np.abs(np.vstack([
                    data["u_pos_abundance"], data["u_neg_abundance"]])), axis=0)
                assert len(y_uncertainty) == len(data)

                neutral = (data["species"] % 1) == 0
                
                # Excitation balance (with neutral lines only)
                excitation_balance = line.fit(
                    data["excitation_potential"][neutral],
                    data["abundance"][neutral],
                    y_uncertainty=y_uncertainty[neutral])

                # Line strength balance (with neutral lines only)
                line_strength_balance = line.fit(
                    np.log(data["equivalent_width"]/data["wavelength"])[neutral],
                    data["abundance"][neutral],
                    y_uncertainty=y_uncertainty[neutral])

                # Ionisation balance.
                if not np.all(np.isfinite(y_uncertainty)):
                    # No uncertainties; all points have equal weight
                    mean_neutral_abundance = np.mean(data[neutral]["abundance"])
                    mean_ionised_abundance = np.mean(data[~neutral]["abundance"])

                else:
                    # Weighted mean
                    mean_neutral_abundance = (data[neutral]["abundance"] \
                        * (y_uncertainty[neutral]**-2)).sum()/(y_uncertainty[neutral]**-2).sum()
                    mean_ionised_abundance = (data[~neutral]["abundance"] \
                        * (y_uncertainty[~neutral]**-2)).sum()/(y_uncertainty[~neutral]**-2).sum()

                ionisation_state = mean_neutral_abundance - mean_ionised_abundance

                # Metallicity state.
                # [TODO] Remove hard coding.
                abundance_state = mean_neutral_abundance - metallicity - 7.50
                results = np.array([excitation_balance, line_strength_balance,
                    ionisation_state, abundance_state])
                
                return results

            op_kwargs = {
                "col_deriv": 1,
                "epsfcn": 0,
                "xtol": 1e-16,
                "maxfev": 100,
                "fprime": utils.stellar_jacobian,
                "full_output": True # Necessary for introspection and provenance.
            }

            ftol_achieved = False
            for i in range(maxiter):
                logger.info("Starting on iteration {0}".format(i + 1))

                logger.info("Initial stellar parameters are Teff = {0:.0f} K, logg = "\
                    "{1:.3f}, [M/H] = {2:.3f}, xi = {3:.3f} km/s".format(
                        initial_stellar_parameters[0], initial_stellar_parameters[2],
                        initial_stellar_parameters[3], initial_stellar_parameters[1]))

                t_init = time()
                op_x, infodict, ier, mesg = op.fsolve(excitation_ionisation_balance,
                    initial_stellar_parameters, **op_kwargs)

                t_elapsed = time() - t_init
                f_op_x = excitation_ionisation_balance(op_x)
                f_op_x_sum = np.abs(f_op_x).sum()

                if np.isfinite(f_op_x_sum) and ftol >= f_op_x_sum:
                    ftol_achieved = True
                    break

                else:
                    initial_stellar_parameters = self.initial_guess_stellar_parameters()
                    logger.warn("Convergence tolerance not achieved {0:.2e} > "\
                        "{1:.2e} in iteration {2}".format(f_op_x_sum, ftol,
                            i + 1))

            op_atomic_data = moogsilent.abfind(op_x[0], op_x[2], op_x[3], op_x[1], 
                line_list_filename, clobber=True)

            if not ftol_achieved:
                logger.warn("Convergence tolerance not achieved after maximum "\
                    "number of iterations ({0})".format(maxiter))
                logger.warn("Optimised parameters are: Teff = N/A K, logg = N/A, "\
                    "[M/H] = N/A, xi = N/A km/s")

                if full_output:
                    return op_x, op_atomic_data, ftol_achieved
                return op_x

            logger.info("Tolerance ({0:.2e} > {1:.2e}) achieved after {2} iteration{3}".format(
                ftol, f_op_x_sum, i + 1, ["", "s"][i > 0]))
            logger.info("Stellar parameter optimisation took {0:.2f} seconds".format(t_elapsed))
            logger.info("Optimised parameters are: Teff = {0:.0f} K, logg = {2:.3f} "\
                "[M/H] = {3:.3f}, xi = {1:.3f} km/s".format(*op_x))
            logger.info("Total squared difference: {0:.2e}".format((f_op_x**2).sum()))

        if full_output:
            return op_x, op_atomic_data, ftol_achieved
        return op_x



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
        self._num_observed_channels = -1
        return None


    @property
    def parameters(self):
        """ Return the model parameters. """

        config = self.config
        any_channel_parameters = (config["model"]["continuum"] \
            or config["model"]["doppler_broadening"])
        if any_channel_parameters and 0 >= self._num_observed_channels:
            raise ValueError("Cannot determine total number of model parameters"\
                " until the model has been supplied some data. We have to give "\
                "model parameters to each observed channel, but we don't know " \
                "how many observed channels there are yet.")
            # [TODO] This is a complicated thing. Perhaps link to online doc page
        
        # [TODO] Check for 1D/<3D> model atmospheres to determine free stellar
        #        parameters.
        parameters = ["teff", "logg", "[M/H]", "xi"]

        # Any redshift parameters?
        if config["model"]["redshift"]:
            # [TODO] Allow individual channel redshifts
            parameters.extend(["z_{0}".format(i) for i in range(self._num_observed_channels)])

        # Any continuum parameters?
        if config["model"]["continuum"]:
            # Don't create unnecessary parameters. Only create parameters for
            # channels that are defined *and* have data.
            num_channels_with_continuum = min([
                len(config["continuum"]["order"]),
                self._num_observed_channels
            ])

            for i, each in enumerate(config["continuum"]["order"]):
                parameters.extend(["c_{0}_{1}".format(i, j) \
                    for j in range(num_channels_with_continuum)])

        # Any doppler broadening?
        if config["model"]["doppler_broadening"]:
            # [TODO] You idiot.
            parameters.extend(["doppler_sigma_{0}".format(i) \
                for i in range(self._num_observed_channels)])

        # Any outlier parameters?
        if config["model"]["outliers"]:
            parameters.extend(["Po", "Vo", "Ys"])

        # [TODO] Abundance parameters.
        # Any element parameters (other than Fe)?
        #if "elements" in config["model"]:
        #    parameters.extend(["log_{0}".format(each) \
        #        for each in set(config["model"]["elements"]).difference("Fe")])

        # [TODO] Any telluric parameters?
        
        return parameters


    def __call__(self, dispersion, synth_kwargs=None, **theta):
        """
        Generate data for the given :math:`\\theta`.

        :param dispersion:
            A list of two-length tuples containing the range of wavelengths to
            synthesise.

        :type dispersion:
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

        if synth_kwargs is None:
            synth_kwargs = {}

        synth_kwargs["full_output"] = False

        # Request the wavelength step to be ~twice the observed pixel sampling.
        #if data is not None:f
        #    synth_kwargs.setdefault("wavelength_steps",
        #        [(0., np.min(np.diff(each.disp))/2., 0.) for each in data])
            
        synthesis_ranges = tuple([(each[0], each[-1]) for each in dispersion])
    
        # Synthesise the model spectra first (in parallel where applicable) and
        # then apply the cheap transformations in serial.

        # Note: Pass keyword arguments so that the cacher works.
        model_spectra = si.synthesise(theta["teff"], theta["logg"], theta["[M/H]"], 
            theta["xi"], wavelengths=synthesis_ranges, **synth_kwargs)

        continuum_spectra = []
        transformed_spectra = []
        for i in range(len(model_spectra)):

            model_spectrum = model_spectra[i].copy()
            
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
            if len(dispersion[i]) > 2:
                expected = np.interp(dispersion[i], model_spectrum[:, 0],
                    model_spectrum[:, 1], left=np.nan, right=np.nan)

                if not np.any(np.isfinite(expected)):
                    raise si.SIException("only nans")

                model_spectrum = np.vstack([dispersion[i], expected]).T

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

            transformed_spectra.append(model_spectrum)

        #if return_continuum:
        #    return (model_spectra, continuum_spectra)
        return transformed_spectra


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

        synth_kwargs = self._get_speedy_synth_kwargs(data, synth_kwargs)

        parameter_guess = {}
        environment = dict(zip(["locals", "globals", "__name__", "__file__",
            "__builtins__"], [None] * 5))

        # Add our distributions and other explicit functions we want to use.
        environment.update({
            "abs": abs,
            "normal": np.random.normal,
            "uniform": np.random.uniform
        })

        self._num_observed_channels = len(data)

        # Any explicit initial guess information?
        if "initial_guess" in self.config:

            sampled_filenames = {}
            for parameter, rule in self.config["initial_guess"].iteritems():
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
        if "priors" in self.config:
            for parameter, rule in self.config["priors"].iteritems():
                if parameter in self.parameters \
                and parameter not in parameter_guess:
                    parameter_guess[parameter] = eval(rule, environment)

        # Apply some default rules in case we haven't guessed these parameters
        # yet.
        default_rules = collections.OrderedDict([
            ("teff", np.random.uniform(3000, 7000)),
            ("logg", np.random.uniform(0, 5)),
            ("[M/H]", np.random.uniform(-2, 0)),
            ("xi", "(1.28 + 3.3e-4 * (teff - 6000) - 0.64 * (logg - 4.5)) "\
                "if logg > 3.5 else 2.70 - 0.509 * logg"),
            ("Po", np.random.uniform(0, 1)),
            ("Vo", abs(np.random.normal(0, 1))),
            ("Ys", np.random.normal(1, 0.01)),
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
            ("doppler_sigma_", "0.20")
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
        any_continuum = len(self.config["continuum"].get("order", [])) > 0 \
            if hasattr(self.config, "continuum") else False

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
                    parameter_guess["[M/H]"],
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
        or not (10000 > theta_dict.get("Vo", 1) > 0) \
        or 0 > theta_dict.get("xi", 1) \
        or 0 > theta_dict.get("teff"):
            return -np.inf

        # Put in priors for channel stubs
        for i in range(self._num_observed_channels):
            if 0 > theta_dict.get("doppler_sigma_{0}".format(i), 1):
                return -np.inf
        
        ln_prior = 0
        for parameter, distribution in self.config.get("priors", {}).iteritems():
            f = eval(distribution, _log_prior_eval_environment)
            ln_prior += f(theta_dict[parameter])
        return ln_prior


    def _log_likelihood(self, theta, data, synth_kwargs):
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

        :param synth_kwargs: [optional]
            Keyword arguments to pass directly to :func:`si.synthesise`.

        :type synth_kwargs:
            dict

        :returns:
            The logarithmic likelihood of the parameters ``theta``.

        :rtype:
            float
        """

        theta_dict = dict(zip(self.parameters, theta))
        try:
            expected = self(dispersion=[spectrum.disp for spectrum in data],
                synth_kwargs=synth_kwargs, **theta_dict)

        except si.SIException:
            # Probably unrealistic astrophysical parameters requested.
            # SI fell over, so we will too.
            # This is OK though: if *all* walkers fell over simultaneously then
            # we would still end up raising an exception
            return -np.inf

        else:
            # No exception, but no spectra either!
            if expected is None:
                return -np.inf

        z = theta_dict.get("z", 0.)
        likelihood = 0
        for observed, model in zip(data, expected):

            mask = self.mask(observed.disp, z)
            chi_sq = (observed.flux - model[:, 1])**2 * observed.ivariance

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

        if likelihood == 0:
            raise a
        return likelihood


    def _log_probability(self, theta, data, synth_kwargs=None):
        """
        Return the logarithmic probability of the GenerativeModel parameters 
        ``theta``.

        :param theta:
            The values of the ``model.parameters``.

        :type theta:
            list-type

        :param data:
            The observed data.

        :type data:
            :class:`specutils.Spectrum1D`

        :param synth_kwargs: [optional]
            Keyword arguments to pass directly to :func:`si.synthesise`.

        :type synth_kwargs:
            dict

        :returns:
            The logarithmic probability of the parameters ``theta``.

        :rtype:
            float
        """
        
        print(theta)
        log_prior = self._log_prior(theta)
        if not np.isfinite(log_prior):
            return log_prior
        log_likelihood = self._log_likelihood(theta, data, synth_kwargs)
        log_probability = log_prior + log_likelihood
        logger.debug("Calculated logarithmic prior, likelihood, and probability"\
            " for {0} to be {1:.3e}, {2:.3e}, and {3:.3e}".format(
            ", ".join(["{0} = {1:.3f}".format(p, v) \
                for p, v in zip(self.parameters, theta)]),
            log_prior, log_likelihood, log_probability))
        return log_probability


    def _get_speedy_synth_kwargs(self, data, synth_kwargs=None):
        """ Get default synth_kwargs that are optimised for speediness. """

        if synth_kwargs is None:
            synth_kwargs = {}

        undersample_rate = self.config["settings"].get("undersample_rate", 1)
        synth_kwargs.setdefault("chunk", True)
        synth_kwargs.setdefault("threads", self.config["settings"].get(
            "max_synth_threads", 1) if "settings" in self.config else 1)
        synth_kwargs.setdefault("wavelength_steps",
            tuple([(wls, wls, wls) for wls \
                in [undersample_rate*np.median(np.diff(s.disp)) for s in data]]))

        return synth_kwargs


    def scatter(self, data, num, synth_kwargs=None):
        """
        Randomly draw points ``num`` from the parameter space and calculate the
        logarithmic probability at each point. The most probable point is
        returned.

        :param data:
            A list of the observed spectra.

        :type data:
            list of :class:`specutils.Spectrum1D` objects

        :param num:
            The number of points to draw from the parameter space.

        :type num:
            int

        :param synth_kwargs: [optional]
            Keyword arguments to pass directly to :func:`si.synthesise`

        :type synth_kwargs:
            dict

        :returns:
            The most probable sampled point.

        :rtype:
            dict
        """

        """
        # [TODO] parallelise
        :param threads: [optional]
            The number of parallel threads to use for random scattering. If
            ``threads`` is a negative number then the number of threads will be
            set by :func:`multiprocessing.cpu_count`.

        :type threads:
            int
        """

        synth_kwargs = self._get_speedy_synth_kwargs(data, synth_kwargs)

        best_theta, highest_log_prob = None, None
        for i in xrange(num):
            theta_dict = self.initial_guess(data, synth_kwargs)
            theta = [theta_dict[p] for p in self.parameters]
            log_prob = self._log_probability(theta, data, synth_kwargs)

            # Is this the best so far?
            if np.isfinite(log_prob) \
            and (highest_log_prob is None or log_prob > highest_log_prob):
                best_theta, highest_log_prob = theta_dict, log_prob

        # If we have gotten to this point and found no reasonable starting
        # point then we should just keep trying until we get *a* start point
        while highest_log_prob is None:
            theta_dict = self.initial_guess(data, synth_kwargs)
            theta = [theta_dict[p] for p in self.parameters]
            log_prob = self._log_probability(theta, data, synth_kwargs)

            if (np.isfinite(log_prob) and log_prob > highest_log_prob):
                best_theta, highest_log_prob = theta_dict, log_prob
                break

        return best_theta


    def optimise(self, data, initial_theta=None, synth_kwargs=None, maxfun=10e3,
        maxiter=10e3, xtol=0.05, ftol=0.01, full_output=False):
        """
        Optimise the logarithmic probability of the GenerativeModel parameters
        given some data.

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

        synth_kwargs = self._get_speedy_synth_kwargs(data, synth_kwargs)

        if initial_theta is None:
            initial_theta = self.initial_guess(data, synth_kwargs=synth_kwargs)

        else:
            self._num_observed_channels = len(data)

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
        """op_theta, op_fopt, op_niter, op_nfunc, op_warnflag = op.fmin(
            lambda theta, data, skw: -self._log_probability(theta, data, skw),
            [initial_theta.get(parameter) for parameter in self.parameters],
            args=(data, synth_kwargs), **op_kwargs)
        """
        op_kwargs = {
            "maxiter": maxiter,
            "epsilon": np.array([10, 0.1, 0.1, 0.1]),
            "norm": -np.inf,
            "gtol": 1e-8
        }
        #op.fmin_bfgs(f, x0, fprime=None, args=(), gtol=1e-05, norm=inf, epsilon=1.4901161193847656e-08, maxiter=None, full_output=0, disp=1, retall=0, callback=None)
        

        #['teff', 'logg', '[M/H]', 'xi'
        #teff, vt, logg, feh
        def min_func(theta, data, skw):
            nlp = -self._log_probability(theta, data, skw)
            #return nlp
            return nlp if np.isfinite(nlp) else 1e3

        result = op.fmin_cg(
            min_func,
            [initial_theta.get(parameter) for parameter in self.parameters],
            args=(data, synth_kwargs), **op_kwargs)
        result = op_theta, op_fopt, op_gopt, op_Bopt, op_nfunc, op_ngrad, warnflag
        #self._opt_warn_message(op_warnflag, op_niter, op_nfunc)

        t_elapsed = time() - t_init
        logger.info("Generative model optimisation took {0:.2f} seconds".format(t_elapsed))

        op_theta_dict = dict(zip(self.parameters, op_theta))

        if full_output:
            return (op_theta_dict, op_fopt, op_niter, op_nfunc, op_warnflag)
        return op_theta_dict


    def infer(self, data, optimised_theta=None, p0=None, walkers=-2, burn=100, 
        sample=100, threads=1, synth_kwargs=None, **kwargs):
        """
        Infer the GenerativeModel parameters given the data.

        :param data:
            A list of the observed spectra.

        :type data:
            list of :class:`specutils.Spectrum1D' objects

        :param optimised_theta: [optional]
            The optimised point :math:`\Theta_{opt}`. The walkers will start from
            a multi-dimensional ball centered around this point. You must supply
            either ``optimised_theta`` or ``p0``, but not both.

        :type optimised_theta:
            dict

        :param p0: [optional]
            The initial starting point for the walkers. You must supply either
            ``optimised_theta`` or ``p0``, but not both.

        :type p0:
            :class:`numpy.ndarray`

        :param walkers: [optional]
            The number of Goodman & Weare (2010) ensemble walkers. 

        :type walkers:
            int

        :param burn: [optional]
            The number of MCMC steps to discard as burn-in.

        :type burn:
            int

        :param sample: [optional]
            The number of MCMC steps to sample the posterior with.
        
        :type sample:
            int

        :param threads: [optional]
            The number of threads to specify to :class:`emcee.EnsembleSampler`.
            Specifying -1 will set the threads to the number of available CPUs.

        :type threads:
            int

        :param kwargs: [optional]
            Keyword arguments to pass directly to :class:`emcee.EnsembleSampler`

        :type kwargs:
            dict

        :returns:
            The posterior quantiles in all parameters, the model sampler, 
            and a ``dict`` containing the mean acceptance fractions, concatenated
            chains and log-probability values for the burn-in and posterior.
        """

        if (optimised_theta is None and p0 is None) \
        or (optimised_theta is not None and p0 is not None):
            raise ValueError("either optimised_theta *or* p0 must be supplied")

        if 0 > threads:
            threads = mp.cpu_count()

        walkers = walkers if walkers > 0 else 2 * abs(walkers) * len(optimised_theta)

        synth_kwargs = self._get_speedy_synth_kwargs(data, synth_kwargs)
        mean_acceptance_fractions = np.zeros((burn + sample))
        sampler = emcee.EnsembleSampler(walkers, len(self.parameters),
            _inferencer, args=(self, data, synth_kwargs), threads=threads, **kwargs)
        
        # Do we need to create p0 from optimised_theta?
        if optimised_theta is not None:
            std = {
                "teff":  10.,
                "logg":  0.05,
                "[M/H]": 0.05,
                "xi":    0.05,
            }
            stds = np.array([std.get(p, 0.01) for p in self.parameters])
            theta = np.array([optimised_theta[p] for p in self.parameters])
            p0 = emcee.utils.sample_ball(theta, stds, size=walkers)

            #p0 = np.array([theta + 1e-4*np.random.randn(len(self.parameters)) \
            #    for i in range(walkers)])
        print("ready")
        for i, (pos, lnprob, rstate) in enumerate(sampler.sample(p0, iterations=burn+sample)):
            mean_acceptance_fractions[i] = np.mean(sampler.acceptance_fraction)
            logger.info("Sampler has finished step {0:.0f} ({1}) of {2:.0f} with "\
                "<a_f> = {3:.3f}, maximum log probability in last step was "\
                "{4:.3e}".format(i + 1, ["burn-in", "sample"][i >= burn], burn + sample,
                    mean_acceptance_fractions[i], np.max(sampler.lnprobability[:, i])))

            if mean_acceptance_fractions[i] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i]))

        """
        # Save the chain and calculated log probabilities for later
        chain, lnprobability = sampler.chain, sampler.lnprobability
        sampler.reset()

        logger.info("Sampling posterior...")

        for j, state in enumerate(sampler.sample(pos, iterations=sample)):
            mean_acceptance_fractions[i + j + 1] = np.mean(sampler.acceptance_fraction)
            logger.info("Sampler has finished step {0:.0f} of sampling with "\
                "<a_f> = {1:.3f}, maximum log probability in last step was "\
                "{2:.3e}".format(j + 1, mean_acceptance_fractions[i + j + 1],
                    np.max(sampler.lnprobability[:, j])))

            if mean_acceptance_fractions[i + j + 1] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i + j + 1]))

        # Concatenate the existing chain and lnprobability with the posteriors
        chain = np.concatenate([chain, sampler.chain], axis=1)
        lnprobability = np.concatenate([lnprobability, sampler.lnprobability], axis=1)
        """

        # Get the quantiles.
        posteriors = dict(zip(self.parameters, map(lambda v: (v[1], v[2]-v[1], v[0]-v[1]),
            zip(*np.percentile(sampler.chain.reshape(-1, len(self.parameters))[-sample*walkers:, :],
                [16, 50, 84], axis=0)))))

        info = {
            "walkers": walkers,
            "burn": burn,
            "sample": sample,
            "mean_acceptance_fractions": mean_acceptance_fractions
        }

        return (posteriors, sampler, info)
