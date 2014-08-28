# coding: utf-8

""" Visualise results """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import matplotlib.pyplot as plt
import numpy as np

import specutils

def spectrum_comparison(data, model, theta, figsize=None, plot_uncertainties=False,
    observed_color=u"k", model_color=u"b", mask_color=u"r", mask_alpha=0.5):
    """
    Produce a comparison plot showing the observed and model spectra.

    :param data:
        A single observed spectrum, or list of observed spectra.

    :type data:
        :class:`specutils.Spectrum1D` object or a list of :class:`specutils.Spectrum1D`
        objects

    :param model:
        The model class.

    :type model:
        :models.Model:

    :param theta:
        The :math:`\Theta` values to use to calculate model spectra for 
        comparison.

    :type theta:
        dict

    :param figsize: [optional]
        A 2 length tuple (width, height) of the figure size in inches.

    :type figsize:
        tuple

    :param plot_uncertainties: [optional]
        Show the uncertainties in the observed data.

    :type plot_uncertainties:
        bool

    [TODO]: other docs

    :returns:
        A spectrum comparison figure.

    :rtype:
        :class:`matplotlib.Figure`
    """

    if isinstance(data, specutils.Spectrum1D):
        data = [data]

    K = len(data)
    if figsize is None:
        figsize = (25, 4 * len(data))

    # Use speedy synth kwargs if the model allows:
    try:
        synth_kwargs = model._get_speedy_synth_kwargs(data)

    except:
        synth_kwargs = {}

    model_spectra = model(dispersion=[s.disp for s in data], 
        synth_kwargs=synth_kwargs, **theta)

    if not isinstance(model_spectra, list):
        model_spectra = [model_spectra]

    fig, axes = plt.subplots(K, figsize=figsize)
    axes = [axes] if K == 1 else axes

    mask = np.array(model.config.get("mask", []))
    # Redshift all mask wavelengths where necessary
    # [TODO] Need to allow for different redshifts in each channel.
    mask *= 1. + theta.get("z", 0)

    for ax, observed_spectrum, model_spectrum \
    in zip(axes, data, model_spectra):

        # Plot the spectra
        if model_spectrum is not None:
            ax.plot(model_spectrum[:, 0], model_spectrum[:, 1], model_color)
        if plot_uncertainties:
            ax.errorbar(observed_spectrum.disp, observed_spectrum.flux,
                yerr=observed_spectrum.variance**0.5, fmt=None, ecolor=observed_color)
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
        ax.set_ylabel("Flux, $F_\lambda$")
    ax.set_xlabel("Wavelength, $\lambda$ ($\AA$)")    

    return fig

