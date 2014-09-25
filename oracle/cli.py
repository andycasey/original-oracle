#!/usr/bin/env python

""" oracle, the suppository of all wisdom """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Necessary for all sub-parsers
import argparse
import logging
import multiprocessing as mp
import os
from time import time

import numpy as np
import matplotlib.pyplot as plt
import triangle

 
import oracle

image_path = lambda s, args: os.path.abspath(os.path.expanduser(s.format(
    args.output_prefix, args.plot_fmt)))

logger = logging.getLogger("oracle")


def _check_plot_format(plotting, plot_fmt):
    """ If plotting is enabled, this checks that the requested plot format is
    actually available on this system. """

    if plotting:
        fig = plt.figure()
        available = map(str.lower, fig.canvas.get_supported_filetypes().keys())
        plt.close(fig)

        if plot_fmt.lower() not in available:
            raise ValueError("Plotting format {0} is not available on this "\
                "system. Available formats are: {1}".format(plot_fmt.lower(),
                    ", ".join(available)))
    return True


def _check_for_existing_files(args):
    """ Create a list of files that the argument parser will create then check
    to see if those files already exist. """

    # TODO
    logger.warn("This may overwrite existing plots!")
    return True

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
    model = oracle.models.GenerativeModel(args.config_filename)
    data = map(oracle.specutils.Spectrum.load, args.spectra_filenames)

    # Sort the spectra from blue to red
    data.sort(key=lambda spectrum: spectrum.disp.mean())
  
    # Make some initial guess of the model parameters.
    initial_guess = model.scatter(data, args.num_scatter_points)
    logger.info("Initial guess for model parameters:")
    for parameter, value in initial_guess.iteritems():
        logger.info("\t{0}: {1:.2f}".format(parameter, value))

    # Create figure showing the initial guess
    if args.plotting:
        path = image_path("{0}-initial.{1}", args)
        fig = oracle.plot.spectrum_comparison(data, model, initial_guess)
        fig.savefig(path)
        plt.close(fig)

        logger.info("Saved model spectrum of initial theta to {0}".format(path))
        
    # Perform the optimisation from the initial guess point.
    optimised_theta = model.optimise(data, initial_guess)

    # Plot a projection showing the optimised theta
    if args.plotting:
        path = image_path("{0}-optimal.{1}", args)
        fig = oracle.plot.spectrum_comparison(data, model, optimised_theta)
        fig.savefig(path)
        plt.close(fig)

        logger.info("Saved model spectrum of optimal theta to {0}".format(path))
    
    # Plot the model parameter values against clock time?
    # [TODO]

    # Inference!
    posteriors, sampler, additional_info = model.infer(data, optimised_theta,
        walkers=100, burn=150, sample=50)
    print(posteriors)
    # Make plots, where necessary.
    if args.plotting:

        # Plot the mean acceptance fractions
        path = image_path("{0}-acceptance.{1}", args)
        fig = oracle.plot.acceptance_fractions(additional_info["mean_acceptance_fractions"])
        fig.savefig(path)
        plt.close(fig)

        # Plot the values of the chains
        path = image_path("{0}-chains.{1}", args)
        fig = oracle.plot.chains(sampler.chain)
        fig.savefig(path)
        plt.close(fig)

        # Make a corner plot [of just the astrophysical parameters]
        # [TODO]
        path = image_path("{0}-corner.{1}", args)
        index = additional_info["burn"] * additional_info["walkers"]
        fig = triangle.corner(sampler.chain.reshape(-1, len(model.parameters))[index:, :],
            labels=oracle.utils.latexify(model.parameters))
        fig.savefig(path)
        plt.close(fig)



    raise a


def solve_theremin(data, args):
    """ The-not-ere-min solver """

    t_init = time()
    model = oracle.models.ThereminModel(args.config_filename)

    # We cannot plot in forked processes.
    plotting = args.plotting if args.threads > 1 else False

    # Optimise the model parameters and plot the transition fits at every 10th
    # iteration of stellar parameters
    parameters, state, converged, transitions, spectra, sampled_parameters, \
        parameter_states = model.optimise(data, plotting=plotting, 
            plot_transition_frequency=10, plot_filename_prefix=args.output_prefix,
            full_output=True)

    # Plot the sampled progress, and whether convergence was achieved?

    if converged:
        # Do something with the results?
        logger.info("Convergence achieved!")

    else:
        print("Did not converge. Final parameters tried were {0} with state {1}"\
            .format(parameters, state))

    return (parameters, state, converged, transitions, spectra, 
        sampled_parameters, parameter_states)



def solve_classical(args):
    """ Classical Model Solver """

    t_init = time()

    # Create the model and load the spectra.
    data = map(oracle.specutils.Spectrum.load, args.spectra_filenames)
    model = oracle.models.StellarSpectrum(args.config_filename, data)

    # Sort the spectra from blue to red
    data.sort(key=lambda spectrum: spectrum.disp.mean())

    # Optimise the individual channel parameters
    optimised_model_parameters = model.optimise(data)

    # Plot the spectrum.
    if args.plotting:
        path = image_path("{0}-optimised.{1}", args)

        model_spectra = [c(s.disp, **theta) for c, s, theta in \
            zip(model.channels, data, optimised_model_parameters)]
        fig = oracle.plot.spectrum_comparison(data, model, model_spectra=model_spectra)
        fig.savefig(path)
        logger.info("Saved figure to {0}".format(path))
        plt.close(fig)


    # Integrate the line profiles and optimise stellar parameters
    atomic_data = model.integrate_profiles(optimised_model_parameters)
    op_x, op_atomic_data, op_converged = model.optimise_stellar_parameters(
        atomic_data, full_output=True, 
        ftol=model.config["classical"].get("tolerance", 4e-3),
        maxiter=model.config["classical"].get("maxiter", 3))

    if not np.isfinite(op_x).all():
        logger.info("Unable to determine stellar parameters")
        logger.info("Full analysis {1}took {0:.2f} seconds".format(time() - t_init,
            ["", "(including plotting) "][args.plotting]))
        return None

    if args.plotting:
        path = image_path("{0}-balance.{1}", args)

        fig = oracle.plot.balance(op_atomic_data)
        fig.savefig(path)
        fig.axes[0].set_title(["Not converged", "Converged"][op_converged])
        logger.info("Saved figure to {0}".format(path))
        plt.close(fig)

    logger.info("Full analysis {1}took {0:.2f} seconds".format(time() - t_init,
        ["", "(including plotting) "][args.plotting]))

    #posterior, sampler, info = model.infer(data, optimised_parameters)
    logger.info("Fin.")


def _wrapper(func, star, arg):
    try:
        return apply(func, (star, arg))
    except:
        logger.info("Wrapper failed!")
    else:
        return None

def solve(args):

    # Check availability of plotting format
    _check_plot_format(args.plotting, args.plot_fmt)

    # Before doing heaps of analysis, look to see if we intend to create some
    # files and if those files already exist -- and we have been told not to 
    # clobber them -- then we should raise an exception now before doing any
    # real work
    _check_for_existing_files(args)

    # If the -r flag is given then we are solving for many stars and (probably)
    # doing things in parallel. Let's deal with the spectrum loading here.
    if args.read_from_filename:
        if len(args.spectra_filenames) > 1:
            raise ValueError("only provide a single filename when using the -r flag")

        with open(args.spectra_filenames[0], "r") as fp:
            spectra_filenames_list = map(str.strip, fp.readlines())

        logger.info("Found spectra for {0} stars from {1}".format(
            len(spectra_filenames_list), args.spectra_filenames[0]))

        stars = []
        for filename_list in spectra_filenames_list:
            data = map(oracle.specutils.Spectrum.load, filename_list.split())
            data.sort(key=lambda spectrum: spectrum.disp.mean())
            stars.append(data)

    else:

        logger.info("Single star assumed. Loading {0} files: {1}".format(
            len(args.spectra_filenames), ", ".join(args.spectra_filenames)))
        # This implicitly assumes that each filename is a single extension 1D file
        data = map(oracle.specutils.Spectrum.load, args.spectra_filenames)
        data.sort(key=lambda spectrum: spectrum.disp.mean())
        stars = [data]

    solver = {
        "theremin": solve_theremin,
        "classical": solve_classical,
        "generative": solve_generative
    }[args.analysis_type]

    threads = args.threads if args.threads > 0 else mp.cpu_count()
    if threads > 1:

        # Parallelise
        logger.info("Pooling {0} threads".format(threads))
        pool = mp.Pool(threads)

        if args.plotting:
            logger.warn("Cannot create plots in threaded processes. Disabling "\
                " plots during analysis. Not all plots will be created.")
        
        processes = [pool.apply_async(_wrapper, args=(solver, star, args)) \
            for star in stars]

        # OK, grab the results
        for i, process in enumerate(processes):
            try:
                result = process.get()

            except:
                logger.exception("Failed to retrieve result for star {0}:".format(i))
                if args.debug:
                    raise
                else:
                    continue

            else:
                # Save the outputs?
                logger.info("Finished with star {0}. Not sure what to do with "\
                    "the outputs yet.".format(i))

        # Close the pool
        pool.close()
        pool.join()

    else:
        # Single thread
        logger.info("Performing serial calculations")
        for i, star in enumerate(stars):
            try:
                result = apply(solver, (data, args))

            except:
                logger.exception("Failed to solve for star {0}".format(i))
                if args.debug:
                    raise
                else:
                    continue

            else:
                logger.info("Finished with star {0}. Not sure what to do with "\
                    "the outputs yet.".format(i))

    


def infer(args):
    raise NotImplementedError("let's focus on optimisation first")




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
        help="Numerically optimise the model parameters, given the data.")

    solve_parser.add_argument("analysis_type", type=str,
        choices=("generative", "classical", "theremin"),
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
        action="store", type=str, default="pdf", help="Format for output plots "\
        "(default: %(default)s)")

    solve_parser.add_argument("--from-filename", "-r", dest="read_from_filename",
        action="store_true", default=False, help="Read input spectra from a filename.")

    solve_parser.add_argument("--threads", "-t", dest="threads", action="store",
        type=int, default=1, help="Number of parallel threads to use.")

    solve_parser.set_defaults(func=solve)


    # Create parser for the infer command
    infer_parser = subparsers.add_parser("infer", parents=[parent_parser],
        help="Infer posterior probability distributions for the model "\
        "parameters, given the data.")

    infer_parser.add_argument("analysis_type", type=str,
        choices=("generative", "classical"),
        help="The kind of analysis to perform. Available options are: classical"\
        " or generative (default: %(default)s).")

    infer_parser.add_argument("config_filename", type=str,
        help="The configuration filename in YAML- or JSON-style formatting.")

    infer_parser.add_argument("spectra_filenames", nargs="+",
        help="Filenames of (observed) spectroscopic data.")

    infer_parser.add_argument("-s", "--num_scatter_points", dest="num_scatter_points",
        action="store", default=1, help="Number of scatter points to sample before"\
        " starting numerical optimisation (default: %(default)s).", type=int)

    infer_parser.add_argument("--no-plots", dest="plotting", action="store_false",
        default=True, help="Disable plotting.")

    infer_parser.add_argument("--output-prefix", "-o", dest="output_prefix",
        default="oracle", help="The filename prefix to use for all output files.")

    infer_parser.add_argument("--plot-format", "-pf", dest="plot_fmt",
        action="store", type=str, default="pdf", help="Format for output plots "\
        "(default: %(default)s). Available formats are (case insensitive):"
        " PDF, JPG, PNG, EPS")
    infer_parser.set_defaults(func=infer)


    # Parse arguments and specify logging level
    args = parser.parse_args()
    handler = logging.FileHandler("{0}.log".format(args.output_prefix), "a")
    formatter = logging.Formatter("%(asctime)s [%(levelname)-7s] %(message)s")
    handler.setFormatter(formatter)

    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    logger.addHandler(handler)

    return args.func(args)

if __name__ == "__main__":
    main()
