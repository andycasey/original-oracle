# coding: utf-8

""" Model a stellar absorption profile """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
import scipy.optimize as op

from scipy.special import wofz
from time import time

import utils
import specutils

logger = logging.getLogger("oracle")


def gaussian(mu, sigma, x):
    """
    Evaluates a Gaussian profile with ``mu`` and ``sigma`` at all values of ``x``.

    :param mu:
        The profile mean.

    :type mu:
        float

    :param sigma:
        The standard deviation of the profile.

    :type sigma:
        float

    :param x:
        The values to calculate the Gaussian profile at.

    :type x:
        :class:`numpy.array`

    :returns:
        An array with values for the calculated absorption profile.

    :rtype:
        :class:`numpy.array`
    """

    return np.exp(-(x - mu)**2 / (2*sigma**2))


def lorentzian(mu, scale, x):
    """
    Evaluates a Lorentzian absorption profile at all `x` values.

    :param mu:
        The centroid of the line.

    :type mu:
        float

    :param scale:
        The scale parameter (also known as gamma).

    :type scale:
        float

    :param x:
        An array of x-points to calculate the profile at.

    :type x:
        :class:`numpy.array`

    :returns:
        The calculated profile points at ``x``.

    :rtype:
        :class:`numpy.array`
    """

    assert scale >= 0

    # 1./np.pi = 0.3183098861837907
    y = scale/((x - mu)**2 + scale**2)
    return (y - np.min(y))/np.ptp(y)
    


def voigt(mu, fwhm, shape, x):
    """
    Evaluates a Voigt absorption profile across the `x`
    values with a given local continuum.

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    :param mu:
        The centroid of the line.

    :param fwhm:
        The full-width half-maximum of the Gaussian profile.

    :param shape:
        The shape parameter.

    :param x:
        An array of x-points.
    """
    n = len(x) if not isinstance(x, float) else 1

    profile = 1. / wofz(np.zeros((n)) + 1j * np.sqrt(np.log(2.0)) * shape).real
    return profile * wofz(2*np.sqrt(np.log(2.0)) * (x - mu)/fwhm \
        + 1j * np.sqrt(np.log(2.0))*shape).real


def pseudo_voigt(mu, sigma, fraction, x):

    z = ((x - mu)/sigma)**2
    gaussian = np.exp(-np.log(2)*z)
    lorentzian = 1.0/(1.0 + z)
    return fraction * lorentzian + (1.0 - fraction) * gaussian 


class AbsorptionProfile(object):

    def __init__(self, wavelength, profile="voigt", mask=None):
        """
        Model an absorption profile in a spectrum.

        :param wavelength:
            Approximate wavelength of the absorption feature.

        :type wavelength:
            float

        :param profile: [optional]
            The type of profile to use. Available profiles are Gaussian or Voigt.

        :type profile:
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
            If the ``profile`` specified is not Voigt or Gaussian, or if the 
            ``wavelength_tolerance``, ``wavelength_contribution`` values are negative.
        """

        self.profile = profile.lower()
        if self.profile not in ("gaussian", "voigt"):
            raise ValueError("profile must be either Gaussian or Voigt")

        #if 0 > wavelength_tolerance:
        #    raise ValueError("wavelength tolerance must be a positive quantity")

        #if 0 > wavelength_contribution:
        #    raise ValueError("wavelength contribution must be a positive quantity")

        self.mask = mask if mask is not None else 1.
        self.outliers = False
        self.approx_wavelength = wavelength
        #self.wavelength_tolerance = wavelength_tolerance
        #self.wavelength_contribution = wavelength_contribution        
        return None


    @property
    def parameters(self):
        """ Return the model parameters. """

        parameters = ["ld"]
        if self.profile == "voigt":
            parameters.extend(["sigma", "shape"])
        elif self.profile == "gaussian":
            parameters.append("sigma")
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
        if self.profile == "voigt":
            depth, shape, sigma = [theta[p] for p in ("ld", "shape", "sigma")]
            flux = continuum * (1. - depth * voigt(wavelength, sigma, shape,
                dispersion))
        elif self.profile == "gaussian":
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
        #if (self.wavelength_tolerance > 0 \
        #and abs(self.approx_wavelength - theta_dict["wl"]) > self.wavelength_tolerance) \
        if not (1 > theta_dict.get("Po", 0.5) > 0) \
        or not (10000 > theta_dict.get("Vo", 1) > 0) \
        or not (1 >= theta_dict.get("ld", 0.5) >= 0) \
        or not (0.5 > theta_dict.get("sigma", 0.0) >= 0) \
        or not (1.0 >= theta_dict.get("shape", 0.0) >= 0):
            return -np.inf

        # It doesn't make a *lot* of sense allowing priors for individual parameters
        # in the AbsorptionProfile model, but we can come back to this.
        # [TODO]
        #ln_prior = 0
        #for parameter, distribution in self.config.get("priors", {}).iteritems():
        #    f = eval(distribution, _log_prior_eval_environment)
        #    ln_prior += f(theta_dict[parameter])
        return 0


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

        """
        # Build an extra mask if only nearby points contribute to the model.
        
        if self.wavelength_contribution > 0:
            wavelength = theta_dict.get("wl", self.approx_wavelength)
            indices = data.disp.searchsorted([
                wavelength - self.wavelength_contribution,
                wavelength + self.wavelength_contribution
            ])
            additional_mask = np.array([np.nan]*len(data.disp))
            additional_mask[indices[0]:indices[1]] = 1.
        
        else:
            additional_mask = 1.
        """


        total_mask = self.mask
        chi_sq = (data.flux - expected)**2 * data.ivariance * total_mask
        
        if not self.outliers:
            likelihood = -0.5 * chi_sq
            finite = np.isfinite(chi_sq)

        else:
            Po, Vo, Yo = [theta_dict.get(each) for each in ("Po", "Vo", "Yo")]
            model_likelihood = -0.5 * (chi_sq - np.log(data.ivariance))
            outlier_ivariance = 1.0/(Vo + data.variance)
            outlier_likelihood = -0.5 * ((data.flux - Yo)**2 * outlier_ivariance \
                * total_mask - np.log(outlier_ivariance))
            likelihood = np.logaddexp(
                np.log(1. - Po) + model_likelihood,
                np.log(Po)      + outlier_likelihood)
            finite = np.isfinite(likelihood)

        return likelihood[finite].sum()


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
            "sigma": 0.05,
            "shape": 0.5,
            "Po": 0.5,
            "Vo": 0.01,
            "Yo": np.median(continuum)
        }
        return dict(zip(self.parameters, [theta[p] for p in self.parameters]))


    @classmethod
    def _opt_warn_message(cls, warnflag, niter, nfunc):
        """ Log a warning message after optimisation. """

        if warnflag > 0:
            message = [
                "Che problem?",
                "Maximum number of function evaluations ({0}) made".format(nfunc),
                "Maximum number of iterations ({0}) made".format(niter)
            ]
            logger.warn("{0}. Optimised values may be inaccurate.".format(
                message[warnflag]))
        return None


    def optimise(self, data, continuum, maxfun=10e3, maxiter=10e3, xtol=1e-12,
        ftol=1e-12, full_output=False):
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
        
        # Let's just provide the data around +/- (wavelength_contribution + wavelength_tolerance)
        """
        if self.wavelength_contribution > 0:
            offset = self.wavelength_contribution + self.wavelength_tolerance
            indices = data.disp.searchsorted([
                self.approx_wavelength - offset,
                self.approx_wavelength + offset
            ])
            _data = specutils.Spectrum1D(disp=data.disp[indices[0]:indices[1]],
                flux=data.flux[indices[0]:indices[1]],
                variance=data.variance[indices[0]:indices[1]])

            if not isinstance(continuum, (int, float)) \
            and len(continuum) == len(data.disp):
                continuum = continuum[indices[0]:indices[1]]

        else:
            _data = data
        """


        t_init = time()
        op_theta, op_fopt, op_niter, op_nfunc, op_warnflag = op.fmin(
            lambda t, d, c: -self._log_probability(t, d, c),
            [initial_theta[p] for p in self.parameters],
            args=(data, continuum), **op_kwargs)
        self._opt_warn_message(op_warnflag, op_niter, op_nfunc)

        t_elapsed = time() - t_init
        logger.info("Optimisation took {0:.2f} seconds".format(t_elapsed))

        if full_output:
            return (op_theta, op_fopt, op_niter, op_funcalls, op_warnflag)
        return op_theta