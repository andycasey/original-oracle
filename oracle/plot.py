# coding: utf-8

""" Visualise results """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import matplotlib.pyplot as plt
import numpy as np

import specutils

def comparison(observed_spectra, model, theta, figsize=None,
    observed_color=u"k", model_color=u"b", mask_color=u"r", mask_alpha=0.5):
    """
    Produce a comparison plot showing the observed and model spectra.

    :param observed_spectra:
        A single observed spectrum, or list of observed spectra.

    :type observed_spectra:
        :class:`specutils.Spectrum1D` object or a list of :class:`specutils.Spectrum1D`
        objects

    :param model:
        The model class.

    :type model:
        :models.Model:

    :param theta:
        The :math:`\Theta` values to use to calculate model spectra for 
        comparison.
    """

    if isinstance(observed_spectra, specutils.Spectrum1D):
        observed_spectra = [observed_spectra]

    K = len(observed_spectra)
    if figsize is None:
        figsize = (25, 3 * len(observed_spectra))

    model_spectra = model(dispersions=[s.disp for s in observed_spectra],
        **theta)

    fig, axes = plt.subplots(K, figsize=figsize)
    axes = [axes] if K == 1 else axes

    mask = np.array(model._configuration.get("mask", []))
    # Redshift all mask wavelengths where necessary
    # [TODO] Need to allow for different redshifts in each channel.
    mask *= 1. + theta.get("z", 0)

    for ax, observed_spectrum, model_spectrum \
    in zip(axes, observed_spectra, model_spectra):

        # Plot the spectra
        ax.plot(model_spectrum[:, 0], model_spectrum[:, 1], model_color)
        ax.plot(observed_spectrum.disp, observed_spectrum.flux, observed_color)
        
        # Show the mask
        obs_start, obs_end = observed_spectrum.disp[0], observed_spectrum.disp[-1]
        for start, end in mask:
            if obs_end >= start and start >= obs_start \
            or obs_end >= end and end >= obs_start:
                # Show the mask in this axes.
                ax.axvspan(start, end, facecolor=mask_color, alpha=mask_alpha,
                    edgecolor='none')

        ax.set_xlim(obs_start, obs_end)
        ax.set_xlabel("Wavelength, $\lambda$ ($\AA$)")
        ax.set_ylabel("Flux, $F_\lambda$")

    return fig

