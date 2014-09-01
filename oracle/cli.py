#!/usr/bin/env python

""" oracle, the suppository of all wisdom """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Necessary for all sub-parsers
import argparse
import logging
import os
from time import time

import numpy as np
import matplotlib.pyplot as plt
import triangle

import plot
import models
import specutils
from utils import latexify

logger = logging.getLogger("oracle")

def _check_for_existing_files(args):
    """ Create a list of files that the argument parser will create then check
    to see if those files already exist. """

    if args.clobber:
        return None

    files_produced = []
    if args.command == "solve":
        files_produced.append("{0}-initial.{1}".format(args.output_prefix,
            args.plot_fmt))

    filename_exists = map(os.path.exists, files_produced)
    if any(filename_exists):
        raise IOError("output filename(s) {0} exist and we have been asked not"\
            " to clobber them.".format(" ".join([f for f in files_produced \
                if os.path.exists(f)])))


def solve_generative(args):
    """ Generative Model Solver """

    # Create the model and load the spectra.
    model = models.GenerativeModel(args.config_filename)
    data = map(specutils.Spectrum.load, args.spectra_filenames)

    # Sort the spectra from blue to red
    data.sort(key=lambda spectrum: spectrum.disp.mean())
  
    # Make some initial guess of the model parameters.
    initial_guess = model.scatter(data, args.num_scatter_points)
    logger.info("Initial guess for model parameters:")
    for parameter, value in initial_guess.iteritems():
        logger.info("\t{0}: {1:.2f}".format(parameter, value))

    image_path = lambda s: os.path.abspath(os.path.expanduser(s.format(
        args.output_prefix, args.plot_fmt)))

    # Create figure showing the initial guess
    if args.plotting:
        path = image_path("{0}-initial.{1}")
        fig = plot.spectrum_comparison(data, model, initial_guess)
        fig.savefig(path)
        plt.close(fig)

        logger.info("Saved model spectrum of initial theta to {0}".format(path))
        
    # Perform the optimisation from the initial guess point.
    optimised_theta = model.optimise(data, initial_guess)

    # Plot a projection showing the optimised theta
    if args.plotting:
        path = image_path("{0}-optimal.{1}")
        fig = plot.spectrum_comparison(data, model, optimised_theta)
        fig.savefig(path)
        plt.close(fig)

        logger.info("Saved model spectrum of optimal theta to {0}".format(path))
    
    # Plot the model parameter values against clock time?
    # [TODO]

    # Inference!
    posteriors, sampler, additional_info = model.infer(data, optimised_theta)
    print(posteriors)
    # Make plots, where necessary.
    if args.plotting:

        # Plot the mean acceptance fractions
        path = image_path("{0}-acceptance.{1}")
        fig = plot.acceptance_fractions(additional_info["mean_acceptance_fractions"])
        fig.savefig(path)
        plt.close(fig)

        # Plot the values of the chains
        path = image_path("{0}-chains.{1}")
        fig = plot.chains(sampler.chain)
        fig.savefig(path)
        plt.close(fig)

        # Make a corner plot [of just the astrophysical parameters]
        # [TODO]
        path = image_path("{0}-corner.{1}")
        index = additional_info["burn"] * additional_info["walkers"]
        fig = triangle.corner(sampler.chain.reshape(-1, len(model.parameters))[index:, :],
            labels=latexify(model.parameters))
        fig.savefig(path)
        plt.close(fig)



    raise a


def solve_classical(args):
    """ Classical Model Solver """

    t_init = time()

    # Create the model and load the spectra.
    data = map(specutils.Spectrum.load, args.spectra_filenames)
    model = models.StellarSpectrum(args.config_filename, data)

    # Sort the spectra from blue to red
    data.sort(key=lambda spectrum: spectrum.disp.mean())

    # Optimise the individual channel parameters
    optimised_model_parameters = model.optimise(data)

    image_path = lambda s: os.path.abspath(os.path.expanduser(s.format(
        args.output_prefix, args.plot_fmt)))

    # Plot the spectrum.
    if args.plotting:
        path = image_path("{0}-optimised.{1}")

        model_spectra = [c(s.disp, **theta) for c, s, theta in zip(model.channels, data, optimised_model_parameters)]
        fig = plot.spectrum_comparison(data, model, model_spectra=model_spectra)
        fig.savefig(path)
        plt.close(fig)

    atomic_data = model.integrate_profiles(optimised_model_parameters)
    op_x, op_atomic_data = model.optimise_stellar_parameters(atomic_data,
        full_output=True)

    logger.info("Full analysis took {0:.2f} seconds".format(time() - t_init))

    raise a
    posterior, sampler, info = model.infer(data, optimised_parameters)




def solve(args):

    # Before doing heaps of analysis, look to see if we intend to create some
    # files and if those files already exist -- and we have been told not to 
    # clobber them -- then we should raise an exception now before doing any
    # real work
    _check_for_existing_files(args)

    if args.analysis_type == "generative":
        solve_generative(args)

    elif args.analysis_type == "classical":
        solve_classical(args)

    else:
        raise NotImplementedError("how did you get here?")






def main():
    """ Parse arguments and execute a particular subparser. """

    parser = argparse.ArgumentParser(description="oracle",
        epilog="Contact Andy Casey <arc@ast.cam.ac.uk> with any questions, and "\
        "open a GitHub issue to report a bug or request a feature.")
    
    # oracle solve generative config.yaml <spectra>
    # (solve does optimise and infer)
    # oracle optimise generative config.yaml <spectra>
    # oracle infer generative config.yaml <spectra>


    # Create subparsers
    subparsers = parser.add_subparsers(title="command", dest="command",
        description="Specify the action to perform.")

    # Create a parent subparser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-v", "--verbose", dest="verbose",
        action="store_true", default=False,
        help="Vebose mode. Logger will print debugging messages.")

    parent_parser.add_argument("--clobber", dest="clobber", action="store_true",
        default=False,
        help="Overwrite existing files if they already exist.")

    parent_parser.add_argument("--debug", dest="debug", action="store_true",
        default=False,
        help="Debug mode. Any suppressed exception during runtime will be re-raised.")


    # Create parser for the solve command
    solve_parser = subparsers.add_parser("solve", parents=[parent_parser],
        help="Compute posterior probability distributions for the model "\
        "parameters, given the data.")

    solve_parser.add_argument("analysis_type", type=str,
        choices=("generative", "classical"),
        help="The kind of analysis to perform. Available options are: classical"\
        " or generative (default: %(default)s).")

    solve_parser.add_argument("config_filename", type=str,
        help="The configuration filename in YAML- or JSON-style formatting.")

    solve_parser.add_argument("spectra_filenames", nargs="+",
        help="Filenames of (observed) spectroscopic data.")

    solve_parser.add_argument("-s", "--num_scatter_points", dest="num_scatter_points",
        action="store", default=1, help="Number of scatter points to sample before"\
        " starting numerical optimisation (default: %(default)s).", type=int)

    solve_parser.add_argument("--no-plots", dest="plotting", action="store_false",
        default=True, help="Disable plotting.")

    solve_parser.add_argument("--output-prefix", "-o", dest="output_prefix",
        default="oracle", help="The filename prefix to use for all output files.")

    solve_parser.add_argument("--plot-format", "-pf", dest="plot_fmt",
        action="store", type=str, default="png", help="Format for output plots "\
        "(default: %(default)s). Available formats are (case insensitive):"
        " PDF, JPG, PNG, EPS")
    solve_parser.set_defaults(func=solve)

    # Parse arguments and specify logging level
    args = parser.parse_args()
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    return args.func(args)

if __name__ == "__main__":
    main()
