# coding: utf-8

""" Functions to handle non-LTE treatment. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import numpy as np
import scipy.interpolate

def interpolate_departure_coefficients(teff, logg, feh, xi, stellar_parameter_grid,
	departure_coefficient_grid, method="linear", fill_value=np.nan, rescale=False):
	"""
	Interpolate non-LTE departure coefficients for some set of stellar parameters
	within the ``stellar_parameter_grid``.

	:param teff:
		The effective temperature of the model atmosphere (Kelvin).

	:type teff:
		float

	:param logg:
		Surface gravity of the model atmosphere.

	:type logg:
		float

	:param feh:
		Metallicity of the model atmosphere.

	:type feh:
		float

	:param xi:
		Microturbulence in the model atmosphere (km/s).

	:type xi:
		float

	:param stellar_parameter_grid:
		The filename of a memory-mapped file containing a grid of stellar
		parameters.

	:type stellar_parameter_grid:
		str

	:param departure_coefficient_grid:
		The filename of a memory-mapped file containing a grid of departure
		coefficients for a given atom.

	:param method: [optional]
		Method of interpolation. One of: nearest, linear, cubic (1-D), cubic (2-D).

	:type method:
		str

	:param fill_value: [optional]
		Value used to fill in for requested points outside of the convex hull of
		the input points.

	:raises ValueError:
		If no deperature coefficients could be interpolated for the given set of
		stellar parameters.

	:returns:
		An array containing interpolated departure coefficients.

	:rtype:
		:class:`numpy.ndarray`
	"""

	memmap_kwargs = {"mode": "r", "dtype": np.double}
	stellar_parameters = np.memmap(stellar_parameter_grid, **memmap_kwargs).reshape(-1, 4)
	departure_coefficients = np.memmap(departure_coefficient_grid,
		**memmap_kwargs).reshape(-1, 80, 87)

	point = np.array([teff, logg, feh, xi])
	# [TODO] Protect Qhull from falling over.
	interpolated_departure_coefficients = scipy.interpolate.LinearNDInterp(
		stellar_parameters, departure_coefficients, point, method="linear",
		fill_value=fill_value)

	if not np.any(np.isfinite(interpolated_departure_coefficients)):
		raise ValueError("""no departure coefficients could be interpolated from
			{0} for ({1}, {2}, {3}, {4})""".format(departure_coefficient_grid,
				teff, logg, feh, xi))

	return interpolated_departure_coefficients