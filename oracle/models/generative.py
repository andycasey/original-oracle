# coding: utf-8

from __future__ import absolute_import, print_function

""" Generative Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import logging
from multiprocessing import cpu_count
from functools import partial
from time import time

import emcee
import numpy as np
from scipy import optimize as op, stats
from scipy.ndimage import gaussian_filter1d

from oracle import si, specutils, utils
from oracle.models.model import Model, _log_prior_eval_environment

logger = logging.getLogger("oracle")


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
            or config["model"]["instrumental_broadening"])
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
        if config["model"]["instrumental_broadening"]:
            # [TODO] You idiot.
            parameters.extend(["instrumental_resolution_{0}".format(i) \
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


    def __call__(self, data, synth_kwargs=None, **theta):
        """
        Generate data for the given :math:`\\theta`.

        :param data:
            The observed spectra.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects

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
            logging.warn("Generate call is missing some theta parameters: {0}"\
                .format(", ".join(missing_parameters)))

        if synth_kwargs is None:
            synth_kwargs = {}

        # Update with defaults
        for key, value in self._get_speedy_synth_kwargs(data).iteritems():
            synth_kwargs.setdefault(key, value)
        
        # Synthesis ranges should contain only regions that are not masked out
        synthesis_ranges = utils.invert_mask(self.config["mask"], data=data,
            padding=0.05)
    
        # Synthesise the model spectra first (in parallel where applicable) and
        # then apply the cheap transformations in serial.

        # Extend wavelength_steps to match length of synthesis_ranges
        if "wavelength_steps" in synth_kwargs:
            # Check
            if isinstance(synth_kwargs["wavelength_steps"][0], (tuple, list)) \
            and len(synth_kwargs["wavelength_steps"]) != len(synthesis_ranges):
                synth_kwargs["wavelength_steps"] = [synth_kwargs["wavelength_steps"][0]] * len(synthesis_ranges)

        # Note: Pass keyword arguments so that the cacher works.
        synthesised_spectra = si.synthesise(theta["teff"], theta["logg"], 
            theta["[M/H]"], theta["xi"], wavelengths=synthesis_ranges, 
            line_list_filename=self.config["GenerativeModel"]["line_list_filename"],
            **synth_kwargs)

        model_spectra = []
        for i, observed in enumerate(data):
            model_contributions = np.zeros(len(observed.disp))
            model_fluxes = np.ones(len(observed.disp))

            pixel_size = np.diff(observed.disp).mean()
            for synthesised_spectrum in synthesised_spectra:

                if  np.any(observed.disp[-1] > synthesised_spectrum[:, 0]) \
                and np.any(observed.disp[0]  < synthesised_spectrum[:, 0]):

                    spectral_portion = synthesised_spectrum.copy()

                    # Apply instrumental broadening
                    resolution = theta.get("instrumental_resolution_{0}".format(i), 0)
                    if resolution > 0:
                        pixel_kernel = (spectral_portion[:,0].mean()/resolution)\
                            /(2.3548200450309493*pixel_size)

                        spectral_portion[:, 1] = gaussian_filter1d(
                            spectral_portion[:, 1], pixel_kernel)

                    # Apply any redshift
                    z = theta["z"] if "z" in theta else theta.get("z_{0}".format(i), 0)
                    spectral_portion[:, 0] *= (1. + z)

                    # Put on the observed pixels
                    indices = observed.disp.searchsorted(spectral_portion[[0, -1], 0])
                    if np.ptp(indices) > 0:
                        model_contributions[indices[0]:indices[-1]] += 1
                        disp = observed.disp[indices[0]:indices[-1]]

                        model_fluxes[indices[0]:indices[1]] *= \
                            np.interp(disp, spectral_portion[:, 0],
                                spectral_portion[:, 1], left=1., right=1.)

            # Set nans for the places we did not synthesise
            model_fluxes[model_contributions == 0] = np.nan

            # Apply continuum
            j, coefficients = 0, []
            while "c_{0}_{1}".format(i, j) in theta:
                coefficients.append(theta["c_{0}_{1}".format(i, j)])
                j += 1

            if len(coefficients) > 0:
                model_fluxes[:, 1] *= np.polyval(coefficients, observed.disp)

            model_spectrum = np.vstack([observed.disp, model_fluxes]).T
            model_spectra.append(model_spectrum)

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
        # TODO No!!
        default_rule_stubs = collections.OrderedDict([
            ("instrumental_resolution_", "28000")
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
            if 0 > theta_dict.get("instrumental_resolution_{0}".format(i), 1):
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
            expected = self(data, synth_kwargs=synth_kwargs, **theta_dict)

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
        
        op_kwargs = {
            #"maxiter": maxiter,
            #"epsilon": np.array([10, 0.1, 0.1, 0.1]),
            #"norm": -np.inf,
            #"gtol": 1e-8
            t
        }
        """
        #op.fmin_bfgs(f, x0, fprime=None, args=(), gtol=1e-05, norm=inf, epsilon=1.4901161193847656e-08, maxiter=None, full_output=0, disp=1, retall=0, callback=None)
        

        #['teff', 'logg', '[M/H]', 'xi'
        #teff, vt, logg, feh
        def min_func(theta, data, skw):
            nlp = -self._log_probability(theta, data, skw)
            #return nlp
            return nlp

        result = op.fmin(
            lambda t, d, skw: -self._log_probability(t, d, skw),
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
            threads = cpu_count()

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