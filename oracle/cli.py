#!/usr/bin/env python

""" operation put spectroscopists out of a job """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Necessary for all sub-parsers
import argparse
import logging

# Necessary for some sub-parsers
import cPickle as pickle
import json
import multiprocessing
import os
import yaml
from time import time

import numpy as np
import pyfits

import models
import specutils


def solve(args):

    # Load the configuration file into a model.
    model = models.GenerativeModel(args.config_filename)
    logger.info("Model parameters: {0}".format(", ".join(model.parameters)))

    # Load the spectra.
    spectra = map(specutils.Spectrum.load, args.spectra_filenames)

    # Make some initial guess of the model parameters based on the prior
    # distributions and some suitable approximations.
    # NB: This is performed automatically by star.optimise if we supply nothing,
    #     but we would like to keep track of what the initial guess was.
    initial_guess = model.initial_guess(spectra)
    logger.info("Initial guess for model parameters:")
    for parameter, value in initial_guess.iteritems():
        logger.info("{0}: {1:.2f}".format(parameter, value))

    # Create figure showing initial guess?
    
    # Perform the optimisation from the initial guess point.
    optimised_theta = star.optimise(spectra, initial_guess)

    # Plot a projection showing the optimised theta

    # Plot the model parameter values against clock time?

    # Inference!
    




def main():
    """ Parse arguments and execute a particular subparser. """

    parser = argparse.ArgumentParser(description="unnamed",
        epilog="Email Andy Casey <arc@ast.cam.ac.uk> with any questions.")
    
    #COMMAND <solve> config.yaml <spectra>

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
    solve_parser.add_argument("config_filename", type=str,
        help="The configuration filename in YAML- or JSON-style formatting.")
    solve_parser.add_argument("spectra_filenames", nargs="+",
        help="Filenames of (observed) spectroscopic data.")
    solve_parser.set_defaults(func=solve)

    # Parse arguments and specify logging level
    args = parser.parse_args()
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    return args.func(args)

if __name__ == "__main__":
    main()
