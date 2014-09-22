# coding: utf-8

from __future__ import absolute_import, print_function

""" Classical Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import itertools
import logging
import multiprocessing as mp
import os
from functools import partial
from time import time

import emcee
import numpy as np
from scipy import optimize as op, integrate, stats
from scipy.ndimage import gaussian_filter1d

from . import line, profiles
from oracle import si, specutils, utils
from oracle.models import Model

logger = logging.getLogger("oracle")

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
        if self.config["classical"].get("use_instrumental_profile", False):
            fwhm_key = "fwhm_{0}"
        else:
            fwhm_key = "fwhm_instrumental"

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
        for the ClassicalModel model.

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



class ClassicalModel(Model):

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
        super(ClassicalModel, self).__init__(configuration)

        self.data = data
        # Create channels.
        self.channels = [SpectralChannel(configuration, i) for i in range(len(data))]

        # Load the line list
        self.atomic_lines = si.io.read_line_list(
            self.config["ClassicalModel"]["clean_line_list_filename"])
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
        all the :class:`SpectralChannel` models associated with this :class:`ClassicalModel`
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
