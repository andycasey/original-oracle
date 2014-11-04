# coding: utf-8

from __future__ import absolute_import, print_function

""" An abstract model class for stellar spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import os
import yaml
import numpy as np
from functools import partial
from scipy import stats

__all__ = ["Model"]

logger = logging.getLogger("oracle")

_prior_env = { 
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

        # Default things that we should have.
        self.config.setdefault("mask", [])
        self.config["mask"] = np.array(self.config["mask"])

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
    def _opt_warn_message(cls, warnflag, niter, nfunc, desc=None):
        """ Log a warning message after optimisation. """

        if warnflag > 0:

            desc = desc if desc is not None else ""
            message = [
                "Che problem?",
                "Maximum number of function evaluations ({0}) made{1}".format(
                    nfunc, desc),
                "Maximum number of iterations ({0}) made{1}".format(niter, desc)
            ]
            logger.warn("{0}. Optimised values may be inaccurate.".format(
                message[warnflag]))
        return None


    def evaluate_lnprior(self, parameter, value):
        """
        Evaluate the logarithmic prior probability as specified in the input
        configuration file for a given parameter at a given theta value.

        :param parameter:
            The model parameter to evaluate.

        :type parameter:
            str
        
        :param value:
            The value of the input theta value to evaluate the prior probability
            for.

        :type value:
            float
        """

        f = eval(self.config.get("priors", {}).get(parameter, "lambda x: 0"),
            _prior_env)
        return f(value)



    def mask(self, dispersion, z=0., fill_value=np.nan, use_cached=False,
        mask_regions=None):
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

        # If kwargs contains a 


        if z == 0 and use_cached: 
            try:
                # Speed hack: see if we have a copy of a rest-frame mask.
                return self._rest_mask
            except AttributeError:
                None

        use_mask_regions = self.config.get("mask", mask_regions)
        if use_mask_regions is not None:
            mask = np.ones(len(dispersion))
            use_mask_regions = np.array(use_mask_regions)
            for start, end in use_mask_regions * (1. + z):
                indices = dispersion.searchsorted([start, end])
                mask[indices[0]:indices[1] + 1] = fill_value

            if z == 0:
                # Speed hack: save a copy of the rest-frame mask.
                self._rest_mask = mask
            return mask

        else:
            return 1.