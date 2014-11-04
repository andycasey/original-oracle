# coding: utf-8


""" An abstract model class for stellar spectra """


from __future__ import absolute_import, print_function

__all__ = ["Model"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import os
import json
import yaml
import numpy as np
from hashlib import md5
from functools import partial
from scipy import stats

logger = logging.getLogger("oracle")

from oracle import utils
from oracle.models import validate

class Model(object):

    # Default configuration
    config = {}

    def __init__(self, configuration, validate=True):
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
            self.config = utils.update_recursively(self.config, configuration)

        else:
            if not os.path.exists(configuration):
                # Probably a string configuration.
                try:
                    supplied_configuration = yaml.load(configuration)

                except:
                    raise IOError("configuration file does not exist or the"\
                        " YAML string provided does not describe a valid "\
                        "dictionary")
                else:
                    # We expect a dictionary.
                    if not isinstance(configuration, dict):
                        raise IOError("configuration file does not exist or the"\
                            " YAML string provided does not describe a valid "\
                            "dictionary")

                    self.config = utils.update_recursively(self.config,
                        supplied_configuration)
            else:
                with open(configuration, "r") as fp:
                    supplied_configuration = yaml.load(fp)

                self.config = utils.update_recursively(self.config,
                    supplied_configuration)

        if validate:
            return validation.validate_configuration(self.config)


    def __str__(self):
        return unicode(self).encode("utf-8")


    def __unicode__(self):
        num_channels = len(self.channels)
        num_models = len(self.grid_points) * num_channels
        num_pixels = sum([len(d) * num_models for d in self.dispersion.values()])
        
        return u"{module}.Model({num_models} {is_cached} models; "\
            "{num_total_parameters} parameters: {num_extra} "\
            "additional parameters, {num_grid_parameters} grid parameters: "\
            "{parameters}; {num_channels} channels: {channels}; ~{num_pixels} "\
            "pixels)".format(module=self.__module__, num_models=num_models, 
            num_channels=num_channels, channels=', '.join(self.channels), 
            num_pixels=utils.human_readable_digit(num_pixels),
            num_total_parameters=len(self.parameters), 
            is_cached=["", "cached"][self.cached],
            num_extra=len(self.parameters) - len(self.grid_points.dtype.names), 
            num_grid_parameters=len(self.grid_points.dtype.names),
            parameters=', '.join(self.grid_points.dtype.names))


    def __repr__(self):
        return u"<{0}.Model object with hash {1} at {2}>".format(self.__module__,
            self.hash[:10], hex(id(self)))


    @property
    def hash(self):
        """ Return a MD5 hash of the JSON-dumped model configuration. """ 
        return md5(json.dumps(self.config).encode("utf-8")).hexdigest()




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