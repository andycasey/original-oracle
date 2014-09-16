# coding: utf-8

""" Visualise results """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging

import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy import stats

import acor
import line
import specutils

logger = logging.getLogger(__name__)

def acceptance_fractions(mean_acceptance_fractions, burn=None, **kwargs):
    """
    Plot the mean acceptance fractions as a function of MCMC step.

    :param mean_acceptance_fractions:
        Mean acceptance fractions at each MCMC step.

    :type mean_acceptance_fractions:
        :class:`numpy.array`

    :param burn:
        The number of MCMC steps discarded as burn in.

    :type burn:
        int

    :returns:
        A figure showing the mean acceptance fractions as a function of MCMC step.

    :rtype:
        :class:`matplotlib.Figure`
    """

    plot_kwargs = {"c": "k", "lw": 2}
    plot_kwargs.update(kwargs)

    fig, ax = plt.subplots()
    N = len(mean_acceptance_fractions)
    ax.plot(np.arange(1, 1 + N), mean_acceptance_fractions, **plot_kwargs)
    if burn is not None and (N >= burn > 0):
        ax.axvline(burn, ":")

    ax.set_xlabel("Step")
    ax.set_ylabel("<f_a>")

    return fig


def spectrum_comparison(data, model, theta=None, model_spectra=None, figsize=None,
    observed_color=u"k", model_color=u"#4682b4", mask_color=u"r", mask_alpha=0.1):
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

    :param theta: [optional]
        The :math:`\Theta` values to use to calculate model spectra for 
        comparison. Either ``theta`` or ``model_spectra`` are required.

    :type theta:
        dict

    :param model_spectra: [optional]
        The model spectra to show. Either ``theta`` or ``model_spectra``
        are required.

    :type model_spectra:
        list of :class:`specutils.Spectrum1D` objects

    :param figsize: [optional]
        A 2 length tuple (width, height) of the figure size in inches.

    :type figsize:
        tuple

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

    if model_spectra is None:
        model_spectra = model(dispersion=[s.disp for s in data], 
            synth_kwargs=synth_kwargs, **theta)

    if not isinstance(model_spectra, list):
        model_spectra = [model_spectra]

    fig, axes = plt.subplots(K, figsize=figsize)
    axes = [axes] if K == 1 else axes

    mask = np.array(model.config.get("mask", []))
    # Redshift all mask wavelengths where necessary
    # [TODO] Need to allow for different redshifts in each channel.
    
    for ax, observed_spectrum, model_spectrum \
    in zip(axes, data, model_spectra):

        # Plot the spectra
        if model_spectrum is not None:
            ax.plot(model_spectrum[:, 0], model_spectrum[:, 1], model_color)

        ax.fill_between(observed_spectrum.disp, 
            observed_spectrum.flux - observed_spectrum.variance**0.5,
            observed_spectrum.flux + observed_spectrum.variance**0.5,
            facecolor="#eeeeee", edgecolor="#666666", zorder=-1)
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


def projection(sampler, model, data, n=100, extents=None, fig=None, figsize=None,
    mask_color="r", mask_alpha=0.1):
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
    width = max([len(each.disp) for each in data])/150.
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
                ax.plot(observed_spectrum.disp, sampled_flux[k], color=u"#4682b4", alpha=0.1)

        # Draw the ML spectra
        ax.plot(observed_spectrum.disp, max_lnprob_flux, color=u"#4682b4", lw=2)

        # Plot the data
        ax.plot(observed_spectrum.disp, observed_spectrum.flux, color="k")

        ax.fill_between(observed_spectrum.disp,
            observed_spectrum.flux - observed_spectrum.variance**0.5,
            observed_spectrum.flux + observed_spectrum.variance**0.5, 
            facecolor='#eeeeee', edgecolor="#666666", zorder=-1)
        
        # Show the mask
        mask = np.array(model.config.get("mask", []))
        obs_start, obs_end = observed_spectrum.disp[0], observed_spectrum.disp[-1]
        for start, end in mask:
            if obs_end >= start and start >= obs_start \
            or obs_end >= end and end >= obs_start:
                # Show the mask in this axes.
                ax.axvspan(start, end, facecolor=mask_color, alpha=mask_alpha,
                    edgecolor='none')


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


def autocorrelation(xs, burn_in, labels=None, fig=None):
    """
    Create a plot showing the autocorrelation in each parameter.

    :param xs:
        The sampled values. This should be a three dimensional array of size
        ``(n_walkers, n_steps, n_parameters)``

    :type xs:
        :class:`numpy.array`

    :param burn_in: [optional]
        The number of steps used for burn-in.

    :type burn_in:
        int

    :param labels: [optional]
        The labels for each parameter.

    :type labels:
        tuple of str

    :param fig: [optional]
        Figure class to use for the plotting.

    :type fig:
        :class:`matplotlib.Figure`

    :returns:
        A figure showing the autocorrelation in each parameter at every MCMC step.

    :rtype:
        :class:`matplotlib.Figure`
    """

    n_walkers, n_steps, K = xs.shape

    factor = 2.0
    lbdim = 0.5 * factor
    trdim = 0.2 * factor
    whspace = 0.10
    width = 15.
    height = factor*K + factor * (K - 1.) * whspace
    dimy = lbdim + height + trdim
    dimx = lbdim + width + trdim

    if fig is None:
        fig, axes = plt.subplots(K, 1, figsize=(dimx, dimy))

    else:
        try:
            axes = np.array(fig.axes).reshape((1, K))
        except:
            raise ValueError("Provided figure has {0} axes, but data has "
                "parameters K={1}".format(len(fig.axes), K))

    lm = lbdim / dimx
    bm = lbdim / dimy
    trm = (lbdim + height) / dimy
    fig.subplots_adjust(left=lm, bottom=bm, right=trm, top=trm,
        wspace=whspace, hspace=whspace)

    for k, ax in enumerate(axes):

        ax.plot(acor.function(np.mean(xs[:, burn_in:, k], axis=0)), color="k")

        if burn_in is not None:
            ax.axvline(burn_in, color="k", linestyle=":")

        ax.set_xlim(0, n_steps)
        if k < K - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Step")

        ax.yaxis.set_major_locator(MaxNLocator(4))
        [l.set_rotation(45) for l in ax.get_yticklabels()]
        if labels is not None:
            ax.set_ylabel(labels[k])
            ax.yaxis.set_label_coords(-0.05, 0.5)

    return fig


def chains(xs, labels=None, truths=None, truth_color=u"#4682b4", burn_in=None,
    alpha=0.5, fig=None):
    """
    Create a plot showing the walker values for each parameter at every step.

    Args:
        xs (array_like) : The samples. This should be a 3D array of size 
            (n_walkers, n_steps, n_parameters)

        labels (iterable, optional) : A list of names for the parameters.

        truths (iterable, optional) : A list of reference values to indicate on
            the plots.

        truth_color (str, optional) : A `matplotlib` style color for the `truths`
            markers.

        burn_in (int, optional) : A reference step to indicate on the plots.

        alpha (float between [0, 1], optional) : Transparency of individual walker
            lines.

        fig (`matplotlib.Figure`, optional) : Overplot onto the provided figure object.
    
    Returns:
        A `matplotlib.Figure` object.
    """

    n_walkers, n_steps, K = xs.shape

    if labels is not None:
        assert len(labels) == K

    if truths is not None:
        assert len(truths) == K

    factor = 2.0
    lbdim = 0.5 * factor
    trdim = 0.2 * factor
    whspace = 0.10
    width = 15.
    height = factor*K + factor * (K - 1.) * whspace
    dimy = lbdim + height + trdim
    dimx = lbdim + width + trdim

    if fig is None:
        fig, axes = plt.subplots(K, 1, figsize=(dimx, dimy))

    else:
        try:
            axes = np.array(fig.axes).reshape((1, K))
        except:
            raise ValueError("Provided figure has {0} axes, but data has "
                "parameters K={1}".format(len(fig.axes), K))

    lm = lbdim / dimx
    bm = lbdim / dimy
    trm = (lbdim + height) / dimy
    fig.subplots_adjust(left=lm, bottom=bm, right=trm, top=trm,
        wspace=whspace, hspace=whspace)

    for k, ax in enumerate(axes):

        for walker in range(n_walkers):
            ax.plot(xs[walker, :, k], color="k", alpha=alpha)

        if burn_in is not None:
            ax.axvline(burn_in, color="k", linestyle=":")

        if truths is not None:
            ax.axhline(truths[k], color=truth_color, lw=2)

        ax.set_xlim(0, n_steps)
        if k < K - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Step")

        ax.yaxis.set_major_locator(MaxNLocator(4))
        [l.set_rotation(45) for l in ax.get_yticklabels()]
        if labels is not None:
            ax.set_ylabel(labels[k])
            ax.yaxis.set_label_coords(-0.05, 0.5)

    return fig


def balance(atomic_data_table, title=None):
    """
    Plot the derived abundances as a function of excitation potential and line
    strength, as typically done in classical analysis approaches.

    :param atomic_data_table:
        A record array table containing the wavelength, species, lower excitation
        potential, oscillator strength (loggf), equivalent width, and abundance
        of the atomic transitions used for stellar parameter determination.

    """

    fig, axes = plt.subplots(2)
    excitation_ax, line_strength_ax = axes

    if len(set(atomic_data_table["species"].astype(int))) > 1:
        logger.warn("Multiple elements found in atomic data table. These will "\
            "all be plotted as the same symbol for the moment.")

    # Seperate neutral and ionised species.
    neutral = (atomic_data_table["species"] % 1) == 0
    ionised = ~neutral

    # Plot the excitation potential axes
    try:
        uncertainties = np.any(np.isfinite(np.vstack([
            atomic_data_table["u_pos_abundance"], atomic_data_table["u_neg_abundance"]])))
    except ValueError:
        uncertainties = False

    if uncertainties:
        # Plot uncertainties
        excitation_ax.errorbar(atomic_data_table["excitation_potential"],
            atomic_data_table["abundance"],
            yerr=(np.abs(atomic_data_table["u_neg_abundance"]), atomic_data_table["u_pos_abundance"]),
            fmt=None, ecolor="k")

    excitation_ax.scatter(atomic_data_table["excitation_potential"][neutral],
        atomic_data_table["abundance"][neutral], facecolor="k", zorder=10)
    excitation_ax.scatter(atomic_data_table["excitation_potential"][ionised],
        atomic_data_table["abundance"][ionised], facecolor="b", zorder=10)

    # Measure slopes by linear regression [TODO] and show them
    y_uncertainty = np.nanmax(np.abs(np.vstack([
        atomic_data_table["u_pos_abundance"], atomic_data_table["u_neg_abundance"]])), axis=0)
    assert len(y_uncertainty) == len(atomic_data_table)

    m, b = line.fit(
        x=atomic_data_table["excitation_potential"][neutral],
        y=atomic_data_table["abundance"][neutral],
        y_uncertainty=y_uncertainty[neutral], full_output=True)[:2]

    x_limits = np.array(excitation_ax.get_xlim())
    y_limits = excitation_ax.get_ylim()
    excitation_ax.plot(x_limits, [np.mean(atomic_data_table["abundance"])] * 2, c="#666666")
    excitation_ax.plot(x_limits, m * x_limits + b, ":", c="k", zorder=-1)
    excitation_ax.set_xlim(x_limits)
    excitation_ax.set_ylim(y_limits)

    excitation_ax.set_xlabel("Lower Excitation Potential (eV)")
    excitation_ax.set_ylabel("$\\log_{\\epsilon}({\\rm Fe})$")

    # [TODO]
    # Quote the slope on the axes.
    logger.info("Slope on excitation balance plot: {0:.4f}".format(m))

    # Plot the line strength axes
    reduced_equivalent_width = np.log(atomic_data_table["equivalent_width"]/atomic_data_table["wavelength"])
    if uncertainties:
        x_pos_uncertainties = np.log(
            (atomic_data_table["equivalent_width"] + atomic_data_table["u_pos_equivalent_width"]) \
            /atomic_data_table["wavelength"]) - reduced_equivalent_width
        x_neg_uncertainties = np.abs(np.log(
            (atomic_data_table["equivalent_width"] + atomic_data_table["u_neg_equivalent_width"]) \
            /atomic_data_table["wavelength"]) - reduced_equivalent_width)

        line_strength_ax.errorbar(reduced_equivalent_width,
            atomic_data_table["abundance"],
            xerr=(x_neg_uncertainties, x_pos_uncertainties),
            yerr=(np.abs(atomic_data_table["u_neg_abundance"]), atomic_data_table["u_pos_abundance"]),
            fmt=None, ecolor="k")

    line_strength_ax.scatter(
        reduced_equivalent_width[neutral], atomic_data_table["abundance"][neutral],
        facecolor="k", zorder=10)
    line_strength_ax.scatter(
        reduced_equivalent_width[ionised], atomic_data_table["abundance"][ionised],
        facecolor="b", zorder=10)

    # Measure slopes by linear regression [TODO] and show them
    m, b = line.fit(
        x=reduced_equivalent_width[neutral],
        y=atomic_data_table["abundance"][neutral],
        y_uncertainty=y_uncertainty[neutral], full_output=True)[:2]

    logger.info("Slope on the reduced equivalent width plot: {0:.4f}".format(m))
    x_limits = np.array(line_strength_ax.get_xlim())
    line_strength_ax.plot(x_limits, [np.mean(atomic_data_table["abundance"])] * 2,
        c="#666666", zorder=-1)
    line_strength_ax.plot(x_limits, m * x_limits + b, ":", c="k", zorder=-1)
    line_strength_ax.set_xlim(x_limits)
    line_strength_ax.set_ylim(y_limits)

    line_strength_ax.set_xlabel("Reduced Equivalent Width")
    line_strength_ax.set_ylabel("$\\log_{\\epsilon}({\\rm Fe})$")

    ionisation_difference = np.nanmean(atomic_data_table["abundance"][neutral]) \
        - np.nanmean(atomic_data_table["abundance"][ionised])
    logger.info("Mean neutral and ionised abundance difference: {0:.3f} dex".format(
        ionisation_difference))

    if title is not None:
        excitation_ax.set_title(title)

    return fig