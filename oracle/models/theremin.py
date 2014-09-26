# coding: utf-8

from __future__ import absolute_import, division, print_function

""" Theremin Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import itertools # why does this make me so happy
import logging
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
from time import time

import numpy as np
import numpy.lib.recfunctions
import scipy.integrate
from scipy import optimize as op, interpolate
from scipy.ndimage import gaussian_filter1d

from oracle import si, specutils, utils
from oracle.plot import transition as plot_transition, balance as plot_balance
from oracle.models import Model, line

logger = logging.getLogger("oracle")


class ConvergenceAchieved(BaseException):
    pass

class ThereminModel(Model):

    _default_constraints = []

    def __init__(self, configuration):
        """
        A class to probabilistically model stellar spectra. This class performs
        on-the-fly synthesis and fits individual absorption lines in order to
        perform an excitation and ionisation balance, while accounting for
        blends.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, or a string containing the configuration. Only YAML-style
            formats are accepted.

        :type configuration:
            str
        """

        super(ThereminModel, self).__init__(configuration)
        return None


    @property
    def parameters(self):
        """
        Return the parameters of the model.
        """

        return ("teff", "logg", "[M/H]", "xi")

    def __call__(self, **theta):

        # For this set of stellar parameters, what are the expected line
        # abundances?
        # This just needs the atomic lines used for SP determination + some
        # solar abundance reference.

        return None

        
    def initial_guess(self):
        """
        Generate an initial guess of the model parameters.
        """

        environment = dict(zip(["locals", "globals", "__name__", "__file__",
            "__builtins__"], [None] * 5))
        environment["np"] = np

        default_rules = collections.OrderedDict([
            ("teff", "np.random.uniform(4000, 6500)"),
            ("logg", "np.random.uniform(1, 5)"),
            ("[M/H]", "np.random.uniform(-2, 0)"),
            ("xi", "(1.28 + 3.3e-4 * (teff - 6000) - 0.64 * (logg - 4.5)) "\
                "if logg > 3.5 else 2.70 - 0.509 * logg")
        ])

        parameter_guess = {}
        for parameter, rule in default_rules.iteritems():

            local_environment = environment.copy()
            local_environment.update(parameter_guess)
            parameter_guess[parameter] = rule if isinstance(rule, (int, float)) \
                else eval(rule, local_environment)

        return [parameter_guess[p] for p in self.parameters]


    def _apply_constraints(self, abundance_data, constraints=None):

        # A default constraint:
        acceptable = np.isfinite(abundance_data["abundance"])
        if constraints is None or len(constraints) == 0:
            return acceptable

        # Create the environment
        environment = dict(zip(["locals", "globals", "__name__", "__file__",
            "__builtins__"], [None] * 5))
        environment["np"] = np
        for column in abundance_data.dtype.names:
            environment[column] = abundance_data[column]

        # Evaluate the constraints
        for constraint in constraints:
            acceptable *= eval(constraint, environment)

        return acceptable


    def _excitation_ionisation_state(self, transitions, metallicity, elements="all"):
        """
        Calculate the excitation and ionisation state for the transitions
        provided.

        :param transitions:
            A record array containing atomic data of all the transitions.

        :type transitions:
            :class:`numpy.core.records.recordarray`

        :param metallicity:
            The model atmosphere metallicity to calculate the excitation and
            ionisation state for.

        :type metallicity:
            float

        :param elements: [optional]
            Filter by particular elements such that only the listed elements 
            will be used for the excitation and ionisation state.

        :type elements:
            tuple
        """

        # species indices
        if isinstance(elements, (str, unicode)) and elements.lower() == "all":
            indices = np.arange(len(transitions))

        else:
            if not isinstance(elements, (list, tuple)):
                raise TypeError("elements must be a list or tuple of strings")

            indices = []
            for element in elements:
                indices.extend(np.where(transitions["atomic_number"] == \
                    utils.atomic_number(element)))
            indices = np.array(list(set(sum(map(list, indices), []))))

        # Excitation
        exc_slope, exc_offset = line.fit(
            x=transitions["excitation_potential"][indices],
            y=transitions["abundance"][indices], full_output=True)[:2]
        
        # Line strength
        lst_slope, lst_offset = line.fit(
            x=np.log(transitions["equivalent_width"]/transitions["wavelength"])[indices],
            y=transitions["abundance"][indices], full_output=True)[:2]

        # Mean ionisation abundance
        neutral = (transitions["ionised"][indices] == 1)
        ionised = (transitions["ionised"][indices] == 2)
        if not any(ionised):
            logger.warn("No acceptable ionised lines!")
            return invalid_response

        ionisation_state = np.median(transitions["abundance"][neutral]) \
            - np.median(transitions["abundance"][ionised])

        # Metallicity state
        metallicity_state = np.median(transitions["abundance"][neutral]) \
            - metallicity

        return np.array([exc_slope, lst_slope, ionisation_state, metallicity_state])


    def optimise(self, data, initial_theta=None, constraints=None,
        state_tolerance=[1e-3, 1e-3, 1e-2, 1e-2], 
        parameter_tolerance=[5, 0.01, 0.01, 0.01], convergence_rule="or",
        op_kwds=None, plotting=True,
        full_output=True, **kwargs):
        """
        Optimise the model parameters given some data.

        :param data:
            The observed stellar spectra.

        :type data:
            list of :class:`specutils.Spectrum1D` objects

        :param constraints: [optional]
            Data quality constraints to apply to the transitions used for stellar
            parameter determination. If None, then it will default to what is
            specified in the configuration file under ThereminModel -> line_constraints.

        :type constraints:
            list
        """

        t_init = time()
        plot_filename_prefix = kwargs.pop("plot_filename_prefix", "opt")
        # A frequency of zero means it won't plot the transitions at each
        # iteration of stellar parameters
        plot_transition_frequency = max([0, kwargs.pop("plot_transition_frequency", 0)])

        invalid_response = np.array([np.inf]*4)
        if initial_theta is None:
            initial_theta = self.initial_guess()

        if constraints is None:
            constraints = self.config["ThereminModel"].get("line_constraints",
                self._default_constraints)
            
        # Some dumb-user (e.g., self-) checking
        convergence_rule = convergence_rule.lower()
        if convergence_rule not in ("or", "and", "state", "parameter"):
            raise ValueError("tolerance rule must be 'or', 'and', 'state', or "\
                "'parameter'")
        assert np.all(np.array(state_tolerance) > 0), "Tolerance values must be"\
            " positive and absolute."
        assert np.all(np.array(parameter_tolerance) > 0), "Tolerance values must"\
            " be positive and absolute."
        assert len(initial_theta) == len(self.parameters)

        message = {
            "state": "state < |{0}|".format(state_tolerance),
            "parameter": "d(parameter) < |{0}|".format(parameter_tolerance),
        }
        message.update({
            "or": "{0} or {1}".format(*message.values()),
            "and": "{0} and {1}".format(*message.values())
        })
        logger.info("Performing excitation and ionisation balance. Requirements"\
            " for convergence are:\n{0}".format(message[convergence_rule]))
        logger.info("Initial stellar parameters: Teff = {0:.0f} K, log(g) = "\
            "{1:.3f}, [M/H] = {2:.3f}, xi = {3:.3f} km/s".format(*initial_theta))

        # At each stellar parameter test we need to re-fit all the lines
        # So we need another model class.
        model = SpectrumModel(self.config, data)

        global sampled_parameters, parameter_states, start_count
        start_count = [0]
        parameter_states = []
        sampled_parameters = []
        @utils.rounder(0, 3, 3, 3)
        @utils.lru_cache(maxsize=25, typed=False)
        def do_balance(*args):

            # Check for a full output
            args = list(args)
            callback = args.pop(-1) if len(args) == len(self.parameters) + 1 else False

            t_a = time()
            global sampled_parameters, parameter_states, start_count

            iteration = len(sampled_parameters) + 1
            logger.info("Performing balance for {0} (iteration #{1})".format(
                args, iteration))

            # Optimise the model parameters given some set of stellar parameters
            try:
                transitions, spectra = model.optimise(*args, full_output=True)
                logger.debug("Finished one iteration")

            except si.SIException:
                # Looks like these parameters are unphysical. Log it, and try 
                # solving again from a new starting point.
                logger.exception("Problem in optimising line abundances for "\
                    "stellar parameters {0}. Trying with new starting point."\
                    .format(args))
                start_count[0] += 1
                raise

            # Apply constraints
            logger.debug("Applying constraints")
            ok = self._apply_constraints(transitions, constraints)

            # Calculate first order derivatives
            logger.debug("Calculating state")
            state = self._excitation_ionisation_state(transitions[ok], args[2])

            # Should we make some plots?
            logger.debug("Plotting..")
            if plotting:
                fig = plot_balance(transitions[ok], title="$T_{{\\rm eff}}$ = "\
                    "{0:.0f} K, log(g) = {1:.3f}, [M/H] = {2:.3f} $\\xi$ = "\
                    "{3:.3f} km/s".format(*args))
                fig.savefig("{0}-iter{1}-state.png".format(plot_filename_prefix,
                    iteration))
                logger.debug("Closing figure")
                plt.close(fig)

                # Have we reached the right plotting transition frequency
                # (by iteration), or is it the first iteration?
                if (plot_transition_frequency > 0 and \
                    (iteration % plot_transition_frequency) == 0) \
                or iteration == 1:
                    # Plot the transitions 

                    logger.info("Plotting transition fits...")
                    for i, (transition, spectrum) in enumerate(zip(transitions, spectra)):

                        if spectrum is None:
                            logger.warn("Skipping plot for transition at {0} because"\
                                " it could not be synthesised.".format(
                                    transition["wavelength"]))
                            continue

                        # Which data index is the transition in?
                        data_index = model._get_data_index(i)

                        path = "{0}-iter{1}-fit-{2:.2f}.png".format(
                            plot_filename_prefix, iteration, transition["wavelength"])
                        fig = plot_transition(data[data_index], self, spectrum, 
                            transition["wavelength"], title="$\chi^2 = {0:.2f}$"\
                            ", R={7:.0f}, [{1} {2:.0f} @ {3:.2f}/H] = {4:.2f}, "\
                            "$EW = {5:.2f} m\AA$ ({6})".format(
                                transition["chi_sq"], 
                                utils.element(transition["atomic_number"]),
                                transition["ionised"],
                                transition["wavelength"],
                                transition["abundance"],
                                transition["equivalent_width"],
                                ["REJECTED", "ACCEPTED"][ok[i]],
                                transition["instrumental_resolution"]))
                        fig.savefig(path)
                        logger.info("Created image {0}".format(path))
                        plt.close(fig)

            logger.info("State for {0} is {1} and took {2:.2f} seconds".format(
                args, state, time() - t_a))

            # Callback calls don't require convergence checks.
            if callback:
                return (state, transitions, spectra)

            # Check for convergence
            state_convergence = np.all(np.less_equal(np.abs(state), state_tolerance))
            if iteration > 1:
                parameter_convergence = np.all(np.less_equal(
                    np.abs(np.array(args) - sampled_parameters[-1]),
                    parameter_tolerance))
            else:
                parameter_convergence = False

            logger.info("State convergence {2}achieved: {0} {3} |{1}|".format(
                state_tolerance, state, ["not ", ""][state_convergence],
                "<>"[state_convergence])) # isn't python awesome?

            if iteration > 1:
                logger.info("Parameter convergence {2}achieved: {0} {3} |{1}|".format(
                    np.abs(parameter_tolerance), np.array(args) - sampled_parameters[-1],
                    ["not ", ""][parameter_convergence], "<>"[parameter_convergence]))

            # Dat logic.
            if (convergence_rule == "state" and state_convergence) \
            or (convergence_rule == "parameter" and parameter_convergence) \
            or (convergence_rule == "or" \
                and (state_convergence or parameter_convergence)) \
            or (convergence_rule == "and" \
                and (state_convergence and parameter_convergence)):
                converged = True
            else:
                converged = False

            # Save the sampled parameters and state
            sampled_parameters.append(args)
            parameter_states.append(state)

            if converged:
                raise ConvergenceAchieved(args, state)

            return state

        op_kwds = {
            "col_deriv": 1,
            "epsfcn": 0,
            "xtol": 1e-16,
            "maxfev": 50,
            "fprime": utils.stellar_jacobian,
            "full_output": True # Necessary for introspection and provenance.
        }
        if op_kwds is not None:
            op_kwds.update(op_kwds)

        while True:
            try:
                opt, infodict, ier, mesg = op.fsolve(
                    lambda a: do_balance(*tuple(a)), initial_theta, **op_kwds)

            except ConvergenceAchieved as (opt, state):
                # Woo hoo!
                # Should we plot transitions?
                logger.info("Convergence achieved. Optimised stellar parameters are"\
                    " Teff = {0:.0f} K, logg = {1:.3f}, [M/H] = {2:.3f}, xi = {3:.3f}"\
                    " km/s".format(*opt))
                logger.info("State was: {0}".format(state))

                # Since we don't have infodict, etc, let's just make them up.
                converged, ier, infodict = True, 0, {}
                mesg = "Convergence achieved after {0} iterations".format(len(sampled_parameters))
                if full_output:
                    opt_with_full_output = list(opt) + [True]
                    state, transitions, spectra = do_balance(*opt_with_full_output)
                break

            except si.SIException:
                # We should catch these and re-raise during debugging.
                if 3 > start_count[0]:
                    initial_theta = self.initial_guess()
                    continue

                else:
                    raise

            else:
                # Should we try again?
                logger.warn("Convergence not achieved: {0}".format(mesg.replace("\n ", "")))

                converged = False
                if full_output:
                    opt_with_full_output = list(opt) + [True]
                state, transitions, spectra = do_balance(*opt_with_full_output)
                break
        
        info = {
            "niterations": len(sampled_parameters) - 1
        }

        if full_output:
            return (opt, state, converged, transitions, spectra, sampled_parameters,
                parameter_states)
        return (opt, state, converged)


    def faster_optimise(self, data, initial_theta=None, constraints=None, update_blended_interval=1,
        plot_transitions=True, plot_filename_prefix="transition-"):
        """
        Optimise the model parameters given some data.

        :param data:
            The observed stellar spectra.

        :type data:
            list of :class:`specutils.Spectrum1D` objects

        :param constraints: [optional]
            Data quality constraints to apply to the transitions to use for
            stellar parameter determination. If no constraints are given, then
            the default ones are employed:

            > default_constraints = [
                "equivalent_width > 5", # mA
                "120 >= equivalent_width",
                "abundance > -5",
                "5 > abundance",
                "chi_sq < 10",
                "equivalent_width/wavelength < np.exp(1e-3)
            ]

        :type constraints:
            list

        :param plot_transitions: [optional]
            Plot the transition figures.

        :type plot_transitions:
            bool

        """
        t_started = time()

        invalid_response = np.array([np.inf]*4)
        if initial_theta is None:
            initial_theta = self.initial_guess()

        if constraints is None:
            constraints = [
                "equivalent_width > 5", # mA 
                "150 >= equivalent_width", # mA
                "warnflag == 0",
                "np.abs(wavelength - 6563) > 100",
                "np.abs(wavelength - 4861) > 100"
 #               "equivalent_width/wavelength < np.exp(1e-3)"
            ]

        logger.info("Initial stellar parameters: Teff = {0:.0f} K, log(g) = "\
            "{2:.3f}, [M/H] = {3:.3f}, xi = {1:.3f} km/s".format(*initial_theta))

        # At each stellar parameter test we need to re-fit all the lines
        # So we need another model class.
        model = SpectrumModel(self.config, data)
        # Since we use a LRU cache, we can do this out here with no slowdown
        #("teff", "xi", "logg", "[M/H]")
        teff, xi, logg, metallicity = initial_theta
        initial_transitions, synthetic_spectra = model.optimise(
            teff, logg, metallicity, xi, full_output=True)
        ok = self._apply_constraints(initial_transitions, constraints)

        # Should we plot the transitions?
        if plot_transitions:
            logger.info("Plotting transition fits...")

            for i, (transition, synthetic_spectrum) \
            in enumerate(zip(initial_transitions, synthetic_spectra)):

                # Which data index is the transition in?
                data_index = model._get_data_index(i)

                path = "{0}fit-{1:.2f}.png".format(plot_filename_prefix,
                    transition["wavelength"])
                fig = plot_transition(data[data_index], self, synthetic_spectrum, 
                    transition["wavelength"], title="$\chi^2 = {0:.2f}$, R={7:.0f}, [{1} "\
                    "{2:.0f} @ {3:.2f}/H] = {4:.2f}, $EW = {5:.2f} m\AA$ ({6})".format(
                        transition["chi_sq"], 
                        utils.element(transition["atomic_number"]),
                        transition["ionised"],
                        transition["wavelength"],
                        transition["abundance"],
                        transition["equivalent_width"],
                        ["REJECTED", "ACCEPTED"][ok[i]], transition["instrumental_resolution"]))
                fig.savefig(path)
                logger.info("Created image {0}".format(path))
                plt.close(fig)

        fig = plot_balance(initial_transitions[ok])
        fig.savefig("balance.png")
        plt.close(fig)
        
        global use_initial_transitions
        use_initial_transitions = [True]

        from oracle import plot

        @utils.rounder(0, 3, 3, 3)
        @utils.lru_cache(maxsize=25, typed=False)
        def fast_balance(effective_temperature, xi, surface_gravity, metallicity):

            t_init = time()
            logger.info("Trying Teff = {0:.0f} K, xi = {1:.3f} km/s, logg = "\
                "{2:.3f}, [M/H] = {3:.3f}".format(effective_temperature, xi,
                    surface_gravity, metallicity))

            transitions = initial_transitions.copy()
            global use_initial_transitions
            if not use_initial_transitions[0]:
                for i, transition in enumerate(transitions):
                    try:
                        transitions[i]["abundance"] = si.find_abundance(
                            transition, effective_temperature, surface_gravity,
                            metallicity, xi, transition["equivalent_width"])

                    except si.SIException:
                        return invalid_response

                    except: 
                        logger.exception("Returning invalid response for parameters"\
                            " {0} K, {1:.2f} km/s, logg = {2:.3f}, [M/H] = {3:.3f}".format(
                                effective_temperature, xi, surface_gravity, metallicity))
                        return invalid_response
            else:
                use_initial_transitions = [False]

            ok = self._apply_constraints(transitions, constraints)
            fig = plot_balance(transitions[ok])
            fig.savefig("balance.png")
            plt.close(fig)

            state = self._excitation_ionisation_state(transitions[ok], metallicity)
            if np.all(np.less_equal(np.abs(state), [1e-3, 1e-3, 1e-3, 1e-3])):
                raise ConvergenceAchieved(state)

            logger.info("State for {0} is {1} and took {2:.2f} seconds".format(
                [effective_temperature, xi, surface_gravity, metallicity], state,
                time() - t_init))
            return state

        op_kwargs = {
            "col_deriv": 1,
            "epsfcn": 0,
            "xtol": 1e-16,
            "maxfev": 20,
            "fprime": utils.stellar_jacobian,
            "full_output": True # Necessary for introspection and provenance.
        }
        t_init = time()
        i, converged = 0, False
        from_theta = initial_theta
        try:
            opt, infodict, ier, mesg = op.fsolve(
                lambda args: fast_balance(*tuple(args)), initial_theta, **op_kwargs)

        except ConvergenceAchieved:
            # woo hoo!
            converged = True

            teff, xi, logg, metallicity = initial_theta
            opt_transitions, synthetic_spectra = model.optimise(
                opt[0], opt[2], opt[3], opt[1], full_output=True)
            ok = self._apply_constraints(initial_transitions, constraints)

            # Should we plot the transitions?
            if plot_transitions:
                logger.info("Plotting transition fits...")

                for i, (transition, synthetic_spectrum) \
                in enumerate(zip(opt_transitions, synthetic_spectra)):

                    # Which data index is the transition in?
                    data_index = model._get_data_index(i)

                    path = "{0}opt-{1:.2f}.png".format(plot_filename_prefix,
                        transition["wavelength"])
                    fig = plot_transition(data[data_index], self, synthetic_spectrum, 
                        transition["wavelength"], title="$\chi^2 = {0:.2f}$, [{1} "\
                        "{2:.0f} @ {3:.2f}/H] = {4:.2f}, $EW = {5:.2f} m\AA$ ({6})".format(
                            transition["chi_sq"], 
                            utils.element(transition["atomic_number"]),
                            transition["ionised"],
                            transition["wavelength"],
                            transition["abundance"],
                            transition["equivalent_width"],
                            ["REJECTED", "ACCEPTED"][ok[i]]))
                    fig.savefig(path)
                    logger.info("Created image {0}".format(path))
                    plt.close(fig)

            if full_output:
                return (opt, converged, opt_transitions, synthetic_spectra)
       
        logger.info("done in ",time() - t_started)

        return (opt, converged)



class SpectrumModel(Model):

    _abundance_key = "[{0}{1:.0f}@{2:.1f}/H]"
    synth_surrounding = 1.

    def __init__(self, configuration, data):
        """
        A class to represent stellar spectra. This class performs on-the-fly 
        synthesis of specific absorption lines in order to perform an excitation
        and ionisation balance, while accounting for blends.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, or a string containing the configuration. Only YAML-style
            formats are accepted.

        :type configuration:
            str
        """

        super(SpectrumModel, self).__init__(configuration)

        if not isinstance(data, (tuple, list)):
            raise TypeError("data must be a list-type of Spectrum1D objects")

        self.data = data
        self._optimised_theta = None

        if not os.path.exists(self.config["ThereminModel"]["clean_line_list_filename"]):
            raise IOError("cannot find clean line list filename {0}".format(
                self.config["ThereminModel"]["clean_line_list_filename"]))
        if not os.path.exists(self.config["ThereminModel"]["blend_line_list_filename"]):
            raise IOError("cannot find blended line list filename {0}".format(
                self.config["ThereminModel"]["blend_line_list_filename"]))

        self.atomic_lines = si.io.read_line_list(
            self.config["ThereminModel"]["clean_line_list_filename"])

        _ = self.parameters
        return None


    def check_line_lists(self, wavelength_tolerance=0.001):
        """
        Check that the transitions in the clean line list do not appear in the
        complete line list.
        """

        complete_filename = self.config["ThereminModel"]["blend_line_list_filename"]
        complete_line_list = si.io.read_line_list(complete_filename)

        for transition in self.atomic_lines:

            # Ensure the transition does not exist in the complete line list.
            matches = np.less_equal(np.abs(complete_line_list["wavelength"] \
                - transition["wavelength"]), wavelength_tolerance)
            matches *= (complete_line_list["atomic_number"] == transition["atomic_number"])
            matches *= (complete_line_list["ionised"] == transition["ionised"])

            if matches.any():
                raise ValueError("{0} {1:.0f} transition at {2:.3f} Angstroms "\
                    "is in both the clean ({3}) and complete ({4}) line list "\
                    .format(utils.element(transition["atomic_number"]),
                        transition["atomic_number"], transition["wavelength"],
                        self.config["ThereminModel"]["clean_line_list_filename"],
                        complete_filename))
        return True


    @property
    def parameters(self):
        """ Return the model parameters. """

        try:
            return self._parameters
        except AttributeError:
            None

        parameters = []
        transition_mapping = {}
        num_channels = len(self.data)
        for i, transition in enumerate(self.atomic_lines):

            wavelength = transition["wavelength"]
            
            # No point adding the transition if it is outside of our observed
            # spectral range.
            in_observed_range = False
            for j, spectrum in enumerate(self.data):
                if (spectrum.disp[-1] >= wavelength >= spectrum.disp[0]):
                    in_observed_range = True
                    if j not in transition_mapping:
                        transition_mapping[j] = []
                    transition_mapping[j].append(i)

            if not in_observed_range:
                # This transition is not in any of the observed spectral ranges
                # so let's ignore it.
                element = utils.element(transition["atomic_number"])
                ionised = transition["ionised"] 
            
                logger.info("Ignoring {0} {1:.0f} transition @ {2:.1f} as it "\
                    "is not in our observed spectral range".format(element,
                        ionised, wavelength))
                continue

            parameters.append(self._format_abundance_key(transition))

        # Dumb check.
        assert len(parameters) == len(set(parameters)), "At least one atomic "\
            "transition was listed more than once."

        # Any redshift parameters?
        if self.config["model"]["redshift"]:
            # [TODO] Allow individual channel redshifts
            parameters.extend(["z_{0}".format(i) for i in range(len(self.data))])

        # Any continuum parameters?
        if self.config["model"]["continuum"]:
            # Don't create unnecessary parameters. Only create parameters for
            # channels that are defined *and* have data.
            num_channels_with_continuum = min([
                len(self.config["continuum"]["order"]),
                num_channels
            ])

            for i, each in enumerate(self.config["continuum"]["order"]):
                parameters.extend(["c_{0}_{1}".format(i, j) \
                    for j in range(num_channels_with_continuum)])

        # Any doppler broadening?
        if self.config["model"]["instrumental_broadening"]:
            # [TODO] You idiot.
            parameters.extend(["instrumental_resolution_{0}".format(i) \
                for i in range(num_channels)])

        self._transition_mapping = transition_mapping
        self._parameters = parameters
        return parameters



    def __call__(self, effective_temperature, surface_gravity, metallicity, xi,
        **theta):
        """
        Produce some synthetic spectra to compare against the data.
        """

        clean_line_list = self.config["ThereminModel"]["clean_line_list_filename"]
        blend_line_list = self.config["ThereminModel"]["blend_line_list_filename"]

        # For each channel
        expected_fluxes = []
        for i, observed in enumerate(self.data):

            expected_flux = np.ones(len(observed.disp))
            pixel_size = np.mean(np.diff(observed.disp))
            z = theta["z"] if "z" in theta else theta.get("z_{0}".format(i), 0)

            # For each transition
            for j in self._transition_mapping[i]:

                element = utils.element(self.atomic_lines[j]["atomic_number"])
                ionised = self.atomic_lines[j]["ionised"]
                wavelength = self.atomic_lines[j]["wavelength"]
                abundance = theta[self._format_abundance_key(self.atomic_lines[j])]
                
                wavelength_region = (
                    wavelength - self.synth_surrounding,
                    wavelength + self.synth_surrounding
                )

                # Synthesise the blended spectrum
                blended_spectrum = si.synthesise(effective_temperature,
                    surface_gravity, metallicity, xi, wavelength_region,
                    wavelength_steps=(pixel_size, pixel_size, pixel_size),
                    line_list_filename=blend_line_list)

                # Convolve the blended spectrum flux
                # [TODO] This should be done by resolution
                if self.config["model"]["instrumental_broadening"]:
                    resolution = theta["instrumental_resolution_{0}".format(i)]
                    kernel = (wavelength/resolution)/(2.3548200450309493*pixel_size)

                    blended_spectrum[:, 1] = gaussian_filter1d(
                        blended_spectrum[:, 1], kernel)

                # Introduce the blended spectrum *only* at the points where
                # the expected flux is 1 (e.g., definitely no synthesis there)
                indices = (expected_flux == 1.)
                expected_flux[indices] *= np.interp(observed.disp[indices],
                    blended_spectrum[:, 0], blended_spectrum[:, 1],
                    left=1., right=1.)

                # Synthesise the transition.
                transition_spectrum = si.synthesise_transition(effective_temperature,
                    surface_gravity, metallicity, xi, tuple(self.atomic_lines[j]),
                    abundance, self.synth_surrounding, (pixel_size, pixel_size, 
                        pixel_size))

                # Smoothing
                # [TODO] Do by resolution
                if self.config["model"]["instrumental_broadening"]:
                    transition_spectrum[:, 1] = gaussian_filter1d(
                        transition_spectrum[:, 1], kernel)

                # Put onto observed dispersion map
                transition_flux = np.interp(observed.disp,
                    transition_spectrum[:, 0] * (1. + z),
                    transition_spectrum[:, 1],
                    left=1., right=1.)

                # [TODO] Only applicable for weak lines
                expected_flux *= transition_flux

            # Apply continuum to expected fluxes
            if self.config["model"]["continuum"]:
                
                # Apply continuum transformation to the model spectra.
                j, coefficients = 0, []
                while "c_{0}_{1}".format(i, j) in theta:
                    coefficients.append(theta["c_{0}_{1}".format(i, j)])
                    j += 1

                if len(coefficients) > 0:
                    expected_flux *= np.polyval(coefficients, observed.disp)

            expected_fluxes.append(expected_flux)

        return expected_fluxes


    def _get_data_index(self, i):

        # Reference the transition to the self.atomic_lines first
        for data_index, transition_indices in self._transition_mapping.iteritems():
            if i in transition_indices:
                break
        else:
            raise ValueError("transition at {0:.3f} Angstroms not found in any "\
                "data index".format(transition["wavelength"]))
        return data_index


    def d(self, effective_temperature, surface_gravity, metallicity,
        xi, resolution=25000):
        """
        Provide an initial guess of all the SpectrumModel model parameters
        given some stellar parameters.
        """

        d = {}
        synthesis_required = \
            (self.config["model"]["redshift"] or self.config["model"]["continuum"])

        line_list_filename = self.config["ThereminModel"]["blend_line_list_filename"]

        # Synthesise spectrum.
        if synthesis_required:
            expected_fluxes = [self._synthesise(i, effective_temperature, 
                surface_gravity, metallicity, xi, line_list_filename=line_list_filename) \
                for i in range(len(self.data))]

        # Cross-correlate with the data to obtain estimate of radial velocity.
        if self.config["model"]["redshift"]:
            raise NotImplementedError

        # Estimate continuum parameters
        if self.config["model"]["continuum"]:
            raise NotImplementedError

        # Estimate doppler broadening parameters from data spectral resolution
        if self.config["model"]["instrumental_broadening"]:
            for i in range(len(self.data)):
                initial_guess["instrumental_resolution_{0}".format(i)] = resolution

        # Remaining parameters should just be abundances.
        remaining_parameters = set(self.parameters).difference(initial_guess)
        initial_guess.update(dict(zip(remaining_parameters, \
            [metallicity] * len(remaining_parameters))))

        return initial_guess


    def _optimise_abundance(self, data, transition, effective_temperature,
        surface_gravity, metallicity, xi, tolerance=0.1, free_broadening=False,
        acceptable_rv_shift=None, full_output=False, **theta):
        """
        Fits an absorption profile leaving the redshift + continuum fixed.
        """

        element = utils.element(transition["atomic_number"])
        ionised = transition["ionised"]
        wavelength = transition["wavelength"]
        
        data_index = self.data.index(data)
        assert data.disp[-1] >= wavelength >= data.disp[0]

        # Get data around that line
        wavelength_region = [
            wavelength - self.synth_surrounding,
            wavelength + self.synth_surrounding
        ]
        indices = data.disp.searchsorted(wavelength_region)
        indices[-1] += 1 # It's an index thing.
        disp = data.disp.__getslice__(*indices)
        flux = data.flux.__getslice__(*indices).copy()
        ivariance = data.ivariance.__getslice__(*indices)
        
        # Create a blended spectrum.
        pixel_size = np.mean(np.diff(disp))
        
        blending_spectrum = si.synthesise(effective_temperature, surface_gravity,
            metallicity, xi, wavelength_region,
            wavelength_steps=(pixel_size, pixel_size, pixel_size),
            line_list_filename=self.config["ThereminModel"]["blend_line_list_filename"])

        # [TODO] err. Real World(tm) is hard.
        mask = self.mask(disp, 0, use_cached=False)
        z = 0.
        continuum = 1.

        # This is some speedyness:
        with si.instance(transition=transition) as siu:

            # Synthesise extremes and do interpolation?
            initial_abundances = np.array([
                metallicity - 2,
                metallicity - 1,
                metallicity,
                metallicity + 1,
                metallicity + 2,
            ])

            initial_fluxes = []
            for initial_abundance in initial_abundances:
                transition_spectrum, equivalent_width, stdout = siu.synthesise_transition(
                    effective_temperature, surface_gravity, metallicity, xi,
                    initial_abundance, self.synth_surrounding, (pixel_size, pixel_size,
                        pixel_size), full_output=True)
                initial_fluxes.append(transition_spectrum[:,1])
                #if equivalent_width == 0:
                #    initial_fluxes.extend([transition_spectrum[:,1]]*(len(initial_abundances) - len(initial_fluxes)))
                #    break

            # Create an interpolator
            interpolator = interpolate.interp1d(initial_abundances,
                np.array(initial_fluxes).T)
            transition_dispersion = transition_spectrum[:,0].copy()

            def approximate_transition(abundance):
                try:
                    return np.vstack([transition_dispersion, interpolator(abundance)]).T
                except ValueError:
                    transition_spectrum, equivalent_width, stdout = siu.synthesise_transition(
                        effective_temperature, surface_gravity, metallicity, xi, 
                        abundance, self.synth_surrounding, (pixel_size, pixel_size,
                            pixel_size), full_output=True)
                    return transition_spectrum

            
            def chi_sq(args, full_output=False):
                abundance = args[0]
 

                #transition_spectrum = approximate_transition(abundance)

                #transition_spectrum, equivalent_width, stdout = siu.synthesise_transition(
                #    effective_temperature, surface_gravity, metallicity, xi, 
                #    abundance, self.synth_surrounding, (pixel_size, pixel_size,
                #        pixel_size), full_output=True)

                if acceptable_rv_shift is not None:
                    v = args[-1]

                    if not (acceptable_rv_shift[1] >= v >= acceptable_rv_shift[0]):
                        return np.inf
                    transition_spectrum[:, 0] *= (1. + v/299792458.0)

                # Any convolution?
                if not free_broadening:
                    resolution = theta.get("instrumental_resolution_{0}".format(data_index), 0)
                    kernel = (wavelength/resolution)/(2.3548200450309493*pixel_size)

                    if resolution > 0:
                        convolved_flux = gaussian_filter1d(transition_spectrum[:, 1],
                            kernel)

                        background_flux = np.interp(disp,
                            blending_spectrum[:, 0] * (1. + z),
                            gaussian_filter1d(blending_spectrum[:, 1], kernel))

                    else:
                        convolved_flux = transition_spectrum[:, 1].copy()
                        background_flux = np.interp(disp,
                            blending_spectrum[:, 0] * (1. + z),
                            blending_spectrum[:, 1])

                else:
                    resolution = args[1]
                    if 50000 > resolution > 20000:
                        kernel = (wavelength/resolution)/(2.3548200450309493*pixel_size)
                        convolved_flux = gaussian_filter1d(transition_spectrum[:, 1],
                            kernel)

                        background_flux = np.interp(disp,
                            blending_spectrum[:, 0] * (1. + z),
                            gaussian_filter1d(blending_spectrum[:, 1], kernel))

                    else:
                        return np.inf

                # Multiply the two spectra together.
                # Redshift and put it on the observed pixel space.
                # [TODO] This approximation is only suitable for weak regimes.
                expected = background_flux * np.interp(disp, 
                    transition_spectrum[:, 0] * (1. + z),
                    convolved_flux,
                    left=np.nan, right=np.nan)

                chi_sq_i = (flux - expected)**2 * ivariance * mask
                dof = np.isfinite(chi_sq_i).sum()
                chi_sq = chi_sq_i[np.isfinite(chi_sq_i)].sum()
                r_chi_sq = chi_sq/(dof - 1)

                if full_output:
                    if not np.any(transition_spectrum[:, 1] < 1.):
                        line_range = [
                            transition_spectrum[0, 0],
                            transition_spectrum[-1, 0]
                        ]

                    else:
                        line_range = [
                            transition_spectrum[transition_spectrum[:, 1] < 1., 0].min(),
                            transition_spectrum[transition_spectrum[:, 1] < 1., 0].max()
                        ]
                    indices = disp.searchsorted(line_range)
                    line_chi_sq_i = chi_sq_i[indices[0]:indices[1]]
                    line_chi_sq = line_chi_sq_i[np.isfinite(line_chi_sq_i)]
                    r_line_chi_sq = line_chi_sq.sum()/(len(line_chi_sq) - 1)
                    
                    return (chi_sq, r_line_chi_sq, equivalent_width, 
                        np.vstack([disp, expected, background_flux]).T)
                return chi_sq

            # Minimise the function.
            p0 = [metallicity]
            if free_broadening:
                p0.append(28000)

            if acceptable_rv_shift is not None:
                p0.append(0.)

            xopt, fopt, niter, nfunc, warnflag = op.fmin(
                chi_sq, p0, xtol=tolerance, ftol=tolerance, disp=False, full_output=True)
            self._opt_warn_message(warnflag, niter, nfunc, "fitting {0} {1:.0f} "\
                "transition at {2:.1f} Angstroms".format(element, ionised, wavelength))
            

            # Plot transitions again


            #epsilon = [0.01, 1000]
            #xopt, fopt, nfunc, ngrad, warnflag = op.fmin_cg(
            #    chi_sq, p0, gtol=0.001, epsilon=epsilon[::-1], disp=False,
            #    full_output=True)
            #direc, niter = None, None

            #epsilon = [1e-2, 1e-8]
            #xopt = op.fmin_bfgs(chi_sq, p0, epsilon=epsilon)
            
            #p0 = np.array(p0)
            #epsilon = np.array([1e-3, 10e6])
            #xopt, fopt, info = op.fmin_l_bfgs_b(chi_sq, p0, approx_grad=True,
            #    epsilon=epsilon, bounds=[(-5, 1), (0, None)], factr=10.)
            #direc, niter, nfunc, warnflag = None, None, None, None
            #self._opt_warn_message(info["warnflag"], info["nit"], info["funcalls"],
            #    "fitting {0} {1:.0f} transition at {2:.1f} Angstroms".format(
            #        element, ionised, wavelength))

            #xopt, nfeval, rc = op.fmin_tnc(chi_sq, p0, approx_grad=True, 
            #    bounds=[(-5, 1), (0, None)], epsilon=epsilon,
            #    disp=False)
            #fopt, direc, niter, nfunc, warnflag = [None] * 5
            #xopt = op.fmin(chi_sq, p0, xtol=tolerance, ftol=tolerance)

            #xopt = op.basinhopping(chi_sq, p0, niter=10)
            #xopt = xopt["x"]

            #xopt, jmin, T, feval, iters, accept, status = op.anneal(chi_sq, p0)
            
            # Necessary:
            xopt = xopt.reshape(-1)
            
            result = [xopt]
            result.extend(chi_sq(xopt, True)[1:])


        if full_output:
            return (result, fopt, direc, niter, nfunc, warnflag)
        return result


    def _format_abundance_key(self, transition):
        wavelength, atomic_number, ionised = [transition[k] \
            for k in ("wavelength", "atomic_number", "ionised")]
        return self._abundance_key.format(
            utils.element(atomic_number), ionised, wavelength)


    def _blending_spectrum(wavelength_region):

        element = utils.element(transition["atomic_number"])
        ionised = transition["ionised"]
        wavelength = transition["wavelength"]
        
        data_index = self.data.index(data)
        assert data.disp[-1] >= wavelength >= data.disp[0]

        # Get data around that line
        wavelength_region = [
            wavelength - self.synth_surrounding,
            wavelength + self.synth_surrounding
        ]
        indices = data.disp.searchsorted(wavelength_region)
        indices[-1] += 1 # It's an index thing.
        disp = data.disp.__getslice__(*indices)
        flux = data.flux.__getslice__(*indices).copy()
        ivariance = data.ivariance.__getslice__(*indices)
        
        # Create a blended spectrum.
        pixel_size = np.mean(np.diff(disp))
        
        return si.synthesise(effective_temperature, surface_gravity,
            metallicity, xi, wavelength_region,
            wavelength_steps=(pixel_size, pixel_size, pixel_size),
            line_list_filename=self.config["ThereminModel"]["blend_line_list_filename"])


    def _optimise_transition(self, data, transition, effective_temperature,
        surface_gravity, metallicity, xi, full_output=True, **initial_theta):

        invalid_response = np.inf

        # Create a blending spectrum
        wavelength = transition["wavelength"]
        wavelength_region = (
            wavelength - self.synth_surrounding,
            wavelength + self.synth_surrounding
        )

        indices = data.disp.searchsorted(wavelength_region)
        indices[-1] += 1

        disp = data.disp.__getslice__(*indices)
        flux = data.flux.__getslice__(*indices).copy() 
        variance = data.variance.__getslice__(*indices)
        ivariance = data.ivariance.__getslice__(*indices)
        pxs = np.diff(disp).mean()

        blending_spectrum = si.synthesise(effective_temperature, surface_gravity,
            metallicity, xi, wavelength_region, wavelength_steps=(pxs, pxs, pxs),
            line_list_filename=self.config["ThereminModel"]["blend_line_list_filename"])

        # [TODO] err. Real World(tm) is hard.
        mask = self.mask(disp, 0, use_cached=False)
        z = 0.
        continuum = 1.

        def chi_sq(args, full_output=False):

            depth, resolution, wl_offset = args

            if not (1 >= depth >= 0) or not (0.05 > wl_offset > -0.05) \
            or 0 > resolution \
            or not np.isfinite(self.evaluate_lnprior("instrumental_resolution", 
                resolution)):
                return invalid_response

            sigma = (wavelength/resolution)/2.35482
            background_flux = continuum * np.interp(disp,
                blending_spectrum[:, 0] * (1. + z),
                gaussian_filter1d(blending_spectrum[:, 1], sigma/pxs))

            # Only allow +/- 0.25 Angstroms to continue to line fitting
            #contribution = np.ones(len(disp)) * np.nan
            #indices = disp.searchsorted([
            #    wavelength - 0.4, wavelength + 0.4])
            #contribution[indices[0]:indices[1] + 1] = 1.
            
            # [TODO] This approximation is only suitable for weak regimes.
            profile = lambda x: depth * np.exp(-(x - (wavelength + wl_offset) \
                * (1. + z))**2 / (2*sigma**2))
            expected = background_flux * (1.0 - profile(disp)) 

            chi_sq_i = (flux - expected)**2 * ivariance * mask
            
            chi_sq = chi_sq_i[np.isfinite(chi_sq_i)].sum()
            dof = np.isfinite(chi_sq_i).sum()
            r_chi_sq = chi_sq/(dof - 1)

            if full_output:
                equivalent_width = depth * sigma * np.sqrt(2 * np.pi) * 10e2
                return (chi_sq, r_chi_sq, equivalent_width,
                    np.vstack([disp, expected, background_flux]).T)
            return r_chi_sq

        # Minimise the function.
        index = disp.searchsorted(wavelength)
        # TODO: specify initial resolution guess
        p0 = [1.0 - flux[index], 28000, 0]#, 0.1, np.mean(variance)]
        if not np.isfinite(p0[0]):
            p0[0] = 1.0

        xopt, fopt, niter, nfunc, warnflag = op.fmin(
            chi_sq, p0, disp=False, full_output=True)
        self._opt_warn_message(warnflag, niter, nfunc, " fitting {0} {1:.0f} "\
            "transition at {2:.1f} Angstroms".format(
                utils.element(transition["atomic_number"]), transition["ionised"],
                wavelength))
        xopt = xopt.reshape(-1)

        if fopt != invalid_response:
            r_chi_sq, equivalent_width, spectrum = chi_sq(xopt, True)[1:]    
        else:
            r_chi_sq, equivalent_width, spectrum = (np.nan, 0, None)

        result = [xopt, r_chi_sq, equivalent_width, spectrum]

        if full_output:
            return (result, fopt, niter, nfunc, warnflag)
        return result


    @utils.rounder(None, 0, 3, 3, 3)
    @utils.lru_cache(maxsize=42, typed=False)
    def optimise(self, effective_temperature, surface_gravity, metallicity, xi,
        initial_guess=None, free_broadening=True, use_previously_optimised=False,
        full_output=False):
        """
        Given some fixed stellar parameters, optimise the SpectrumModel parameters.

        :param effective_temperature:
            The given stellar effective temperature.

        :type effective_temperature:
            float

        :param surface_gravity:
            The given surface gravity of the model star.

        :type surface_gravity:
            float

        :param metallicity:
            The given metallicity of the model star.

        :type metallicity:
            float

        :param xi:
            The given microturbulence of the model star.

        :type xi:
            float

        :param initial_guess: [optional]
            The initial guess of the SpectrumModel parameters.

        :type initial_guess:
            dict
        """


        if initial_guess is None:
            if use_previously_optimised and self._optimised_theta is not None:
                logger.info("Using previously optimised theta as initial guess:"\
                    "{0}".format(self._optimised_theta))
                initial_theta = self._optimised_theta

            else:
                initial_theta = {}
                logger.warn("Remember that this bit is here because I'm not 100%"\
                    " sure whether it's necessary or not.")


        # Create the final dictionary result
        optimal_theta = {}
        optimal_theta.update(initial_theta)

        # Fit each line/neighbouring individually using the current estimate of
        # redshift, continuum, etc.

        t_init = time()
        rows = []
        spectra = []
        for i, observed_spectrum in enumerate(self.data):

            for k, transition in enumerate((self.atomic_lines[j] \
                for j in self._transition_mapping[i])):
                
                result, fopt, niter, nfunc, warnflag = self._optimise_transition(
                    observed_spectrum, transition, effective_temperature,
                    surface_gravity, metallicity, xi, full_output=True, 
                    free_broadening=free_broadening, **initial_theta)

                xopt, r_chi_sq, equivalent_width, spectrum = result

                # Calculate abundance
                abundance = si.find_abundance(transition, effective_temperature,
                    surface_gravity, metallicity, xi, equivalent_width)

                # Save and announce the results
                key = self._format_abundance_key(transition)
                optimal_theta[key] = xopt[0]
                message = [
                    "(Invalid line: EW = 0)",
                    "EW = {0:.1f} mA".format(equivalent_width)
                ][equivalent_width > 0]
                logger.info("Found {1} = {0:.2f} with chi^2 = {2:.1f} {3}".format(
                    abundance, key, r_chi_sq, message))

                # Save to table
                rows.append([abundance, xopt[0], xopt[1], r_chi_sq, equivalent_width,
                    warnflag])
                spectra.append(spectrum)

        logger.info("Optimised line fitting in {0:.2f} seconds".format(
            time() - t_init))
        # Without instance (+/- 1 Angstroms): 141 seconds [no plotting]
        # With instance, ftol=0.0001 (+/- 1 Angstroms): 129 seconds [no plotting]
        # With instance, ftol=0.01: 94.64 seconds [no plotting]
        # With instance, ftol=0.10 and precision 2 in abundance for LRU_cache: 82 seconds
        # With instance, ftol=0.10, precision 2 in abundance w/ LRU, initial R~28k: 57 seconds

        # Create a record array with the equivalent width, abundance, smoothing?
        # reduced equivalent width.
        logger.debug("Creating record array")
        transition_indices = sum([self._transition_mapping[i] for i in xrange(len(self.data))], [])
        results = numpy.lib.recfunctions.append_fields(
            np.array([self.atomic_lines[i] for i in transition_indices]),
            ("abundance", "line_depth", "instrumental_resolution", "chi_sq", 
                "equivalent_width", "warnflag"), np.array(rows).T, usemask=False)
        logger.debug("Returning results")

        if full_output:
            return (results, spectra)
        return results



def _calculate_cog(transition, teff, logg, metallicity, xi, abundances, **kwargs):

    # ensure abundances are sorted high->low
    kwargs["full_output"] = True
    assert abundances[0] > abundances[-1]
    equivalent_widths = np.zeros(len(abundances))
    with si.instance(transition=transition) as siu:
        for k, abundance in enumerate(abundances):
            try:
                spectrum, equivalent_width, stdout = siu.synthesise_transition(
                    teff, logg, metallicity, xi, abundance, **kwargs)

            except si.SIException:
                equivalent_widths[k] = np.nan

            else:
                equivalent_widths[k] = equivalent_width
                if equivalent_width == 0:
                    break
    return equivalent_widths




def calculate_cog_tables(model, abundances=[-3, -2.5, -2, -1.5, -1, -0.5, 0.0, 0.5, 1.0],
    limits=None, spacings=None, threads=1, **kwargs):
    """
    At each permutation of stellar parameters, calculate the abundances for all
    lines given some measured equivalent widths.
    """

    use_limits = {
        "teff": [4000, 6500],
        "logg": [0, 5],
        "[M/H]": [-3, 1.0],
        "xi": [0.5, 4.0]
    }

    use_spacings = {
        "teff": 250,
        "logg": 0.5,
        "[M/H]": 0.5,
        "xi": 0.5
    }
    
    if limits is not None:
        use_limits.update(limits)

    if spacings is not None:
        use_spacings.update(spacings)

    # Ensure abundances are sorted high->low
    abundances = np.sort(abundances)[::-1]

    # Create all possible combinations
    points = []
    parameters = ("teff", "logg", "[M/H]", "xi")
    for parameter in parameters:
        points.append(np.arange(use_limits[parameter][0], use_limits[parameter][1]
            + use_spacings[parameter], use_spacings[parameter]))
    
    kwargs["full_output"] = True

    # We listify the stellar parameter combinations 
    stellar_parameter_combinations = list(itertools.product(*points))

    # Create an array to store all the equivalent widths
    n_c, n_ew, n_abund = map(len, 
        [stellar_parameter_combinations, model.atomic_lines, abundances])
    equivalent_widths = np.zeros((n_c, n_ew, n_abund))

    for i, stellar_parameters in enumerate(stellar_parameter_combinations):

        logger.info("At point {0:.0f}/{1:.0f}: {2}".format(i+1, n_c, stellar_parameters))
        teff, logg, metallicity, xi = stellar_parameters

        # Calculate the equivalent widths for all transitions at this set of
        # stellar parameters
        if threads > 1:
            processes = []
            pool = mp.Pool(threads)

            for j, transition in enumerate(model.atomic_lines):
                process = pool.apply_async(_calculate_cog,
                    args=(transition, teff, logg, metallicity, xi, abundances))
                processes.append([j, process])

            for j, process in processes:
                row = process.get()[::-1]
                #if np.all(np.isfinite(equivalent_widths)):
                #    assert np.all(np.diff(equivalent_widths[::-1]) >= 0)
    
                equivalent_widths[i, j, :] = row

            pool.close()
            pool.join()

        else:
            for j, transition in enumerate(model.atomic_lines):
                equivalent_widths[i, j, :] = _calculate_cog(
                    transition, teff, logg, metallicity, xi, abundances)[::-1]


    # Set minimum equivalent width as zero
    equivalent_widths = np.clip(equivalent_widths, 0, np.inf)
    return (stellar_parameter_combinations, equivalent_widths)


