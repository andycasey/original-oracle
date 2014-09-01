# coding: utf-8

""" Probabilistic Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
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
import si
import specutils

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
                    raise IOError("configuration file does not exist or the "\
                        "string configuration is not valid YAML")
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

        mask_regions = self.config.get("mask", None)
        if mask_regions is not None:
            mask = np.ones(len(dispersion))
            mask_regions = np.array(mask_regions)
            for start, end in mask_regions * (1. + z):
                indices = dispersion.searchsorted([start, end])
                mask[indices[0]:indices[1] + 1] = fill_value

            return mask

        else:
            return 1.


def _optimiser(_class, *args):
    return _class.optimise(*args)


def _inferencer(theta, _class, *args):
    return _class._log_probability(theta, *args)
        

class SpectralChannel(Model):
    
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

        # We can't get the exact model parameters unless we know the observed
        # dispersion range.
        transition_parameters = ("shape", "ld", "sigma")
        for i, (wavelength, species) in enumerate(config["balance"]["atomic_lines"]):
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

        # Model the absorption profiles!
        for i, (wavelength, species) in enumerate(self.config["balance"]["atomic_lines"]):
            # Negative depth means don't model the line.
            depth = theta.get("ld_{0}".format(i), -1)
            if depth > 0: 
                sigma = theta["sigma_{0}".format(i)]
                shape = theta["shape_{0}".format(i)]
                wavelength = theta.get("wl_{0}".format(i), wavelength)

                # Over what dispersion do we care?
                indices = dispersion.searchsorted([
                    wavelength - 1, 
                    wavelength + 1
                ])
                x, y = dispersion.__getslice__(*indices), flux.__getslice__(*indices)
                flux.__setslice__(indices[0], indices[1],
                    y * (1 - depth * profiles.voigt(wavelength, sigma, shape,x)))


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

        if not (1 > theta_dict.get("Po", 0.5) > 0) \
        or not (10 > theta_dict.get("Vo", 1) > 0):
            return -np.inf

        for i, (wavelength, species) in enumerate(self.config["balance"]["atomic_lines"]):
            #if ("wl_{0}".format(i) in self.parameters \
            #and abs(wavelength - theta_dict["wl_{0}".format(i)]) > wavelength_tolerance) \
            if not (1 >= theta_dict.get("ld_{0}".format(i), 0) >= -1) \
            or not (1 >= theta_dict.get("shape_{0}".format(i), 0.5) >= 0) \
            or not (1 >= theta_dict.get("sigma_{0}".format(i), 0.5) > 0):

                return -np.inf
        
        ln_prior = 0
        for parameter, distribution in self.config.get("priors", {}).iteritems():
            f = eval(distribution, _log_prior_eval_environment)
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
        
        mask = self.mask(data.disp, theta_dict.get("z", 0))
        chi_sq = (data.flux - expected)**2 * data.ivariance * mask
        
        stellar_likelihood = -0.5 * (chi_sq - np.log(data.ivariance))

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


    def _log_probability(self, theta, data, verbose=False):
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
            if verbose:
                print(log_prior)
            return log_prior
        probability = log_prior + self._log_likelihood(theta, data)
        if verbose:
            print(probability)
        return probability


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
        if absorption_profile_parameters:
            for i, (wavelength, species) in enumerate(self.config["balance"]["atomic_lines"]):
                if "ld_{0}".format(i) in self.parameters:
                    
                    index = data.disp.searchsorted(wavelength)
                    line_continuum = continuum if isinstance(continuum, (int, float)) \
                        else continuum[index]

                    initial_theta["ld_{0}".format(i)] = \
                        line_continuum - data.flux[index]
                    initial_theta["sigma_{0}".format(i)] = 0.05
                    initial_theta["shape_{0}".format(i)] = 0.50

        # Some final defaults.
        initial_theta.update({
            "Po": 0.50,
            "Vo": 0.01,
            "Yo": np.median(continuum),
            "z": 0. # [TODO]
        })
        return dict(zip(self.parameters, [initial_theta.get(p, np.nan) \
            for p in self.parameters]))


    def optimise(self, data, initial_theta=None, maxfun=10e3, maxiter=10e3,
        xtol=0.01, ftol=0.01, threads=-1, full_output=False):
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

        t_init = time()
        if initial_theta is None:
            initial_theta = self.initial_guess(data,
                absorption_profile_parameters=True)

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
        for i, (wavelength, species) in enumerate(self.config["balance"]["atomic_lines"]):
            if "ld_{0}".format(i) in self.parameters:
                
                # Optimise the parameters of the absorption profile.

                # Only supply +/- 1 or 2 Angstroms
                # [TODO]
                indices = data.disp.searchsorted([
                    wavelength - 0.3,
                    wavelength + 0.3
                ])
                _data = specutils.Spectrum1D(disp=data.disp.__getslice__(*indices),
                    flux=data.flux.__getslice__(*indices),
                    variance=data.variance.__getslice__(*indices))
                _mask = mask.__getslice__(*indices)
                _data.flux[~np.isfinite(_mask)] = np.nan

                if isinstance(continuum, (int, float)):
                    _continuum = continuum

                else:
                    _continuum = continuum.__getslice__(*indices)

                profile = profiles.AbsorptionProfile(wavelength)
                process = pool.apply_async(_optimiser, args=(profile, _data, _continuum))
                processes.append([i, profile, process])
        
        for i, profile, process in processes:
            xopt = dict(zip(profile.parameters, process.get()))
            for parameter in set(xopt).difference(["Po", "Vo", "Yo"]):
                pre_opt_theta["{0}_{1}".format(parameter, i)] = xopt[parameter]
        
        pool.close()
        pool.join()
        

        # If some optimisation failed and those values are at the edge of what
        # is acceptable, then we will need to alter those values. If we don't,
        # then half of the points will have starting values greater than what
        # is allowed by the priors. For each parameter that has an initial value
        # on the edge of what's allowed, we will lose 50% of the walkers.

        #initial_theta = self._clip_points_for_sampling(initial_theta)


        import matplotlib.pyplot as plt
        
        print("MODEL PARAMETERS")
        for k, v in initial_theta.iteritems():
            print("{0} = {1:.3f}".format(k, v))

        print(self.parameters)

        fig = plot.spectrum_comparison([data], self, initial_theta)
        ax = fig.axes[0]
        

        # OK, now optimise (minimise) the negative log probability, maybe?
        if num_coefficients > 0 \
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
                lambda theta, data, verbose: -self._log_probability(theta, data, verbose),
                [pre_opt_theta.get(parameter) for parameter in self.parameters],
                args=(data, True, ), **op_kwargs)
            self._opt_warn_message(op_warnflag, op_niter, op_nfunc)
            opt_theta = dict(zip(self.parameters, op_x))

            #opt_spectra = self(dispersion=data.disp, **op_theta)
            #ax.plot(data.disp, opt_spectra[:,1], 'g', label='opt')

        else:
            op_x = [pre_opt_theta.get(parameter) for parameter in self.parameters]
            opt_theta = pre_opt_theta

        opt_spectra = self(dispersion=data.disp, **pre_opt_theta)
        ax.plot(data.disp, opt_spectra[:,1], 'r', label='opt')



        t_elapsed = time() - t_init
        logger.info("Optimisation took {0:.2f} seconds".format(t_elapsed))

        if full_output:
            return (op_theta, op_fopt, op_niter, op_funcalls, op_warnflag)
        return opt_theta


    def _clip_points_for_sampling(self, theta, tolerance=0.01):
        """
        Alter the starting (usually optimised) values of theta such that the
        initial walker values will not be mostly junk. If we don't do this then
        the MCMC will have extremely low acceptance fractions and likely not
        converge in a Hubble time.
        """

        clipped_theta = {}
        clipped_theta.update(theta)
        wavelength_tolerance = self.absorption_profile_kwargs["wavelength_tolerance"]
        for i, (parameter, value) in enumerate(theta.iteritems()):
            if parameter.startswith("wl_"):
                index = int(parameter.split("_")[1])
                approximate_wavelength = self.config["balance"]["atomic_lines"][index][0]

                clipped_theta[parameter] = approximate_wavelength
                #print("diff", abs(approximate_wavelength - value), abs(wavelength_tolerance - (approximate_wavelength - value)))
                #if abs(wavelength_tolerance - (approximate_wavelength - value)) < 0.03:
                #    clipped_theta[parameter] = approximate_wavelength

            elif parameter.startswith("sigma_"):
                clipped_theta[parameter] = np.clip(value, 0.03, 1-0.03)

            elif parameter.startswith("fwhm_"):
                clipped_theta[parameter] = np.clip(value, 0.03, 2-0.03)

            elif parameter.startswith("ld_") or parameter == "Po":
                clipped_theta[parameter] = np.clip(value, -1 + 0.03, 1-0.03)

            elif parameter == "Vo" and 0.03 > value:
                clipped_theta[parameter] = 0.03

        return clipped_theta

    def infer(self, data, opt_theta, walkers=50, burn=400, sample=100, threads=24,
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

        walkers = walkers if walkers > 0 else 2*len(opt_theta)

        p0 = emcee.utils.sample_ball(opt_theta, [1e-3]*len(opt_theta),
            size=walkers)

        p0_orig = p0.copy()

        # Do reflections so that we have a reasonable initial starting point.
        limits = {
            "Po": [0, 1],
            "Vo": [0, None],
            "ld": [-1, 1],
            "shape": [0, 1],
            "sigma": [0, 1]
        }
        for i, parameter in enumerate(self.parameters):
            limit = limits.get(parameter.split("_")[0], False)
            if limit:
                p0[:, i] = utils.reflect_about(p0[:, i], limit)

        sampler = emcee.EnsembleSampler(walkers, len(self.parameters),
            _inferencer, args=(self, data), threads=threads, **kwargs)

        mean_acceptance_fractions = np.zeros((burn + sample))
        for i, (pos, lnprob, rstate) in enumerate(sampler.sample(p0, \
            iterations=burn)):
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


    def initial_guess(self, data, **kwargs):
        return [channel.initial_guess(spectrum, **kwargs) \
            for channel, spectrum in zip(self.channels, data)]



    def optimise(self, data, initial_theta=None,
        threads=-1, op_line_kwargs=None):
        """
        Optimise the model parameters given the data.

        """

        # TODO CHANGE INIITAL_THETA TO SINGLE DICT AND NOT LIST OF DICTS
        if initial_theta is None:
            initial_theta = self.initial_guess(data)

        results = [c.optimise(d) for c, d in zip(self.channels, data)]

        table = self.integrate_profiles(results)
        result = self.optimise_stellar_parameters(table)

        raise a

        line_kwargs = { "xtol": 1e-8, "ftol": 1e-8, "maxfun": 10e4, "maxiter": 10e4,
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
        results = []
        for i in range(max_global_iterations):
            for j, channel in enumerate(self.channels):

                result = channel.optimise()
                results.append(result)
                
                # Change the channel continuum coefficient parameter names

                # Change the redshift parameter name if necessary

        raise a

        config = self.config["balance"]
        
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
        

    def integrate_profiles(self, optimal_theta, xlimits=1):
        """
        Integrate measured absorption profiles to find equivalent widths for all
        atomic absorption lines in the current model.

        :param optimal_theta:
            The optimised model parameters for each channel.

        :type optimal_theta:
            list of dicts

        :returns:
            An array of atomic transitions and integrated equivalent widths.
        """

        num_atomic_lines = len(self.config["balance"]["atomic_lines"])
        tabular_results = np.zeros((num_atomic_lines, 3))

        # Fill 'er up.
        tabular_results[:, :2] = np.array(self.config["balance"]["atomic_lines"])
        tabular_results[:, 2] = np.nan

        for optimal_channel_theta in optimal_theta:

            for parameter, value in optimal_channel_theta.iteritems():
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

                    # Integrate the profile.
                    depth = value
                    shape = optimal_channel_theta["shape_{0}".format(transition_index)]
                    sigma = optimal_channel_theta["sigma_{0}".format(transition_index)]
                    y = lambda x: depth * profiles.voigt(wavelength, sigma, shape, x)
                    equivalent_width, integration_error_est = integrate.quad(y,
                        wavelength - xlimits, wavelength + xlimits)

                    # Save it (in mA)
                    tabular_results[transition_index, -1] = equivalent_width * 1000.

        return tabular_results


    def optimise_stellar_parameters(self, equivalent_width_table, 
        initial_stellar_parameters=None, model_outliers=True):
        """
        Optimise stellar parameters (Teff, logg, [Fe/H], xi) given some 
        measured equivalent widths.
        """

        observed = equivalent_width_table[:, 2]
        atomic_data = np.array(self.config["balance"]["atomic_lines"])
        initial_stellar_parameters = [5800.0, 2.0, -1.0, 1.05]

        def func(theta):
            temperature, logg, metallicity, xi = theta
            returned_data = np.array([each[0] for each in si.equivalent_width(
                temperature, logg, metallicity, xi, atomic_data[:, 0], 
                atomic_data[:, 1], threads=4, full_output=True)])

            returned_data[0 >= returned_data[:, 2], 2] = np.nan

            finite = np.isfinite(returned_data[:, 2])
            #slopes = stats.linregress(x=returned_data[finite, 1],
            #    y=returned_data[finite, 2])

            #print(slopes)
            #raise a
            print(returned_data[:, 0] - atomic_data[:, 0])
            expected = returned_data[:, 2]
            difference = (observed - expected)
            total = difference[np.isfinite(difference)]
            print(theta, (total**2).sum())
            return (total**2).sum()

        op_kwargs = {
            "maxfun": 10e3,
            "maxiter": 10e4,
            "xtol": 0.01,
            "ftol": 1e-6,
            "disp": False,
            "full_output": True # Necessary for introspection and provenance.
        }
        t_init = time()
        #op_theta, op_fopt, op_niter, op_nfunc, op_warnflag
        result = op.fmin(func,
            initial_stellar_parameters, **op_kwargs)
        #self._opt_warn_message(op_warnflag, op_niter, op_nfunc)

        raise a
        t_elapsed = time() - t_init
        logger.info("Optimisation took {0:.2f} seconds".format(t_elapsed))

        raise a





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
            parameters.append("z")

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
            
        synthesis_ranges = [(each[0], each[-1]) for each in dispersion]
    
        # Synthesise the model spectra first (in parallel where applicable) and
        # then apply the cheap transformations in serial.
        model_spectra = si.synthesise(theta["teff"], theta["logg"], theta["[M/H]"],
            theta["xi"], synthesis_ranges, **synth_kwargs)

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
            if len(dispersion[i]) > 2:
                expected = np.interp(dispersion[i], model_spectrum[:, 0],
                    model_spectrum[:, 1], left=np.nan, right=np.nan)
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

            model_spectra[i] = model_spectrum

        #if return_continuum:
        #    return (model_spectra, continuum_spectra)
        return model_spectra


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
            ("doppler_sigma_", "0.05")
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


    def _log_probability(self, theta, data, synth_kwargs=None):
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

        :param synth_kwargs: [optional]
            Keyword arguments to pass directly to :func:`si.synthesise`.

        :type synth_kwargs:
            dict

        :returns:
            The logarithmic probability of the parameters ``theta``.

        :rtype:
            float
        """
        
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
            [(wls, wls, wls) for wls in [undersample_rate*np.median(np.diff(s.disp)) \
                for s in data]])
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


    def optimise(self, data, initial_theta=None, synth_kwargs=None, maxfun=500,
        maxiter=500, xtol=0.10, ftol=10e10, full_output=False):
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
        op_theta, op_fopt, op_niter, op_nfunc, op_warnflag = op.fmin(
            lambda theta, data, skw: -self._log_probability(theta, data, skw),
            [initial_theta.get(parameter) for parameter in self.parameters],
            args=(data, synth_kwargs), **op_kwargs)
        self._opt_warn_message(op_warnflag, op_niter, op_nfunc)

        t_elapsed = time() - t_init
        logger.info("Optimisation took {0:.2f} seconds".format(t_elapsed))

        op_theta_dict = dict(zip(self.parameters, op_theta))

        if full_output:
            return (op_theta_dict, op_fopt, op_niter, op_nfunc, op_warnflag)
        return op_theta_dict


    def infer(self, data, optimised_theta=None, p0=None, walkers=20, burn=200, 
        sample=100, threads=-1, **kwargs):
        """
        Infer the model parameters given the data.

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

        mean_acceptance_fractions = np.zeros((burn + sample))
        sampler = emcee.EnsembleSampler(walkers, len(self.parameters),
            _inferencer, args=(self, data), threads=threads, **kwargs)
        
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
       

        for i, (pos, lnprob, rstate) in enumerate(sampler.sample(p0, iterations=burn)):
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
                    np.max(sampler.lnprobability[:, j])))

            if mean_acceptance_fractions[i + j + 1] in (0, 1):
                raise RuntimeError("mean acceptance fraction is {0:.0f}!".format(
                    mean_acceptance_fractions[i + j + 1]))

        # Concatenate the existing chain and lnprobability with the posteriors
        chain = np.concatenate([chain, sampler.chain], axis=1)
        lnprobability = np.concatenate([lnprobability, sampler.lnprobability], axis=1)

        # Get the quantiles.
        posteriors = dict(zip(self.parameters,map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
            zip(*np.percentile(sampler.chain.reshape(-1, len(self.parameters)),
                [16, 50, 84], axis=0)))))

        info = {
            "chain": chain,
            "lnprobability": lnprobability,
            "mean_acceptance_fractions": mean_acceptance_fractions
        }

        return (posteriors, sampler, info)
