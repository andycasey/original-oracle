# coding: utf-8

""" Visualise results """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
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


def projection(sampler, model, data, n=100, extents=None, fig=None, figsize=None):
    """
    Project the maximum likelihood values and sampled posterior points as spectra.

    :param sampler:
        The sampler employed.

    :type sampler:
        :class:`emcee.EnsembleSampler`

    :param model:
        The model employed.

    :type model:
        :class:`sick.models.Model`

    :param data:
        The observed spectra.

    :type data:
        iterable of :class:`sick.specutils.Spectrum1D` objects

    :param extents: [optional]
        The wavelength extents to plot for each channel in the form of [(min_chan_1,
        max_chan_1), ..., (min_chan_N, max_chan_N)]
    
    :type extents:
        tuple or None

    :param fig: [optional]
        Overplot onto the provided figure object.

    :type fig:
        :class:`matplotlib.Figure` or None
    
    :param figsize: [optional]
        The figure size (x-dimension, y-dimension) in inches.

    :type figsize:
        tuple or None

    :raises ValueError:
        If a ``fig`` is provided with the incorrect number of axes.

    :raise TypeError:
        If the ``data`` are not provided in the correct type.

    :returns:
        The projection figure.

    :rtype:
        :class:`maplotlib.Figure`
    """

    if not isinstance(data, (tuple, list)) or \
    any([not isinstance(each, specutils.Spectrum1D) for each in data]):
        raise TypeError("Data must be a list-type of Spectrum1D objects.")

    K = len(data)

    factor = 3.0
    lbdim = 0.5 * factor
    trdim = 0.2 * factor
    whspace = 0.10
    width = 8.
    height = factor*K + factor * (K - 1.) * whspace
    dimy = lbdim + height + trdim
    dimx = lbdim + width + trdim

    if figsize is None:
        figsize = (dimx, dimy)
    if fig is None:
        fig, axes = plt.subplots(K, 1, figsize=figsize)

    else:
        try:
            axes = np.array(fig.axes).reshape((1, K))
        except:
            raise ValueError("Provided figure has {0} axes, but data has "
                "parameters K={1}".format(len(fig.axes), K))

    # Find the most probable sampled theta and compute spectra for it
    max_lnprob_index = np.argmax(sampler.lnprobability.flatten())
    max_lnprob_theta = sampler.flatchain[max_lnprob_index]
    max_lnprob_fluxes = [model(dispersion=[spectrum.disp for spectrum in data], 
        **dict(zip(model.parameters, max_lnprob_theta)))[:,1]]

    if n > 0:
        # Draw samples from sampler.chain and compute spectra for them
        sampled_fluxes = []
        n_samples = len(sampler.flatchain)

        for i in range(n):
            sampled_theta = dict(zip(
                model.parameters,
                sampler.flatchain[np.random.randint(0, n_samples)]
            ))

            try:
                sampler_flux = [model(dispersion=[s.disp for s in data],
                    **sampled_theta)[:,1]]
            except:
                continue
            else:
                sampled_fluxes.append(sampler_flux)
    
    if len(data) == 1:
        axes = [axes]

    for k, (ax, max_lnprob_flux, observed_spectrum) in enumerate(zip(axes, max_lnprob_fluxes, data)):

        # Draw the random samples from the chain
        if n > 0:
            for sampled_flux in sampled_fluxes:
                ax.plot(observed_spectrum.disp, sampled_flux[k], color="#666666")

        # Draw the ML spectra
        ax.plot(observed_spectrum.disp, max_lnprob_flux, color="r", lw=2)

        # Plot the data
        ax.plot(observed_spectrum.disp, observed_spectrum.flux, color="k")

        # By default only show common overlap between the model and spectral data
        if extents is None:
            finite_data = np.isfinite(observed_spectrum.flux)
            finite_model = np.isfinite(max_lnprob_flux)

            x_extent = [
                np.max([observed_spectrum.disp[indices][0]  for indices in (finite_model, finite_data)]),
                np.min([observed_spectrum.disp[indices][-1] for indices in (finite_model, finite_data)]),
            ]

            indices = observed_spectrum.disp.searchsorted(x_extent)
            finite_observed_flux = observed_spectrum.flux[indices[0]:indices[1]]
            y_extent = [
                0.9 * np.min(finite_observed_flux[np.isfinite(finite_observed_flux)]),
                1.1 * np.max(finite_observed_flux[np.isfinite(finite_observed_flux)])
            ]
            ax.set_xlim(x_extent)
            ax.set_ylim(y_extent)

        else:
            ax.set_xlim(extents[k][0])
            ax.set_ylim(extents[k][1])

        # Labels and ticks
        if not (k < K - 1):
            ax.set_xlabel("Wavelength, $\lambda$ ($\AA$)")

        ax.set_ylabel("Flux, $F_\lambda$")
        ax.yaxis.set_label_coords(-0.05, 0.5)

        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        [l.set_rotation(45) for l in ax.get_yticklabels()]

    return fig

