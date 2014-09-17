# coding: utf-8

from __future__ import absolute_import, print_function

""" Theremin Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import logging
import os

import numpy as np
from scipy import optimize as op
from scipy.ndimage import gaussian_filter1d

from oracle import si, specutils, utils
from oracle.models import Model

logger = logging.getLogger("oracle")


class ThereminModel(Model):

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
            ("teff", "np.random.uniform(3000, 7000)"),
            ("logg", "np.random.uniform(0, 5)"),
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


    def optimise(self, data):
        """
        Optimise the model parameters given some data.

        :param data:
            The observed stellar spectra.

        :type data:
            list of :class:`specutils.Spectrum1D` objects
        """

        # At each stellar parameter test we need to re-fit all the lines
        # So we need another model class.

        model_spectra = SpectrumModel(self.configuration, data)

        initial_theta = model_spectra.initial_guess()

        def min(teff, logg, mh, xi):

            model_data = model_spectra.optimise(teff, logg, mh, xi)
            expected = self(teff, logg, mh, xi)

            # Use jacobian approximation for 




        return None



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

        assert os.path.exists(self.config["ThereminModel"]["clean_line_list_filename"])
        assert os.path.exists(self.config["ThereminModel"]["blend_line_list_filename"])

        self.atomic_lines = si.io.read_line_list(
            self.config["ThereminModel"]["clean_line_list_filename"])

        _ = self.parameters


        """
        # Check the line list filename is provided, and match against the atomic
        # transitions listed
        line_list_filename = self.config["ThereminModel"]["line_list_filename"]
        logger.info("Checking line list {0}".format(line_list_filename))

        wl_tolerance = 0.001
        all_transitions = si.io.read_line_list(line_list_filename)
        for transition in self.config["ThereminModel"]["atomic_lines"]:

            # Check the transition exists in the line list.
            wavelength, species = transition[:2]
            matches = np.less_equal(np.abs(all_transitions["wavelength"] - wavelength),
                wl_tolerance)
            matches *= (all_transitions["atomic_number"] == int(species))
            matches *= (all_transitions["ionised"] == int((species % 1)*10 + 1))

            if not matches.any():
                raise ValueError("{0} {1:.0f} transition at {2:.3f} Angstroms "\
                    "not in line list {3}".format(utils.element(species),
                        ((species % 1)*10)+1, wavelength, line_list_filename))
        """

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
            element = utils.element(transition["atomic_number"])
            ionised = transition["ionised"] 
            
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
                logger.info("Ignoring {0} {1:.0f} transition @ {2:.1f} as it "\
                    "is not in our observed spectral range".format(element,
                        ionised, wavelength))
                continue

            parameters.append(self._abundance_key.format(element, ionised,
                wavelength))

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
        if self.config["model"]["doppler_broadening"]:
            # [TODO] You idiot.
            parameters.extend(["doppler_sigma_{0}".format(i) \
                for i in range(num_channels)])

        self._transition_mapping = transition_mapping
        self._parameters = parameters
        return parameters


    def _synthesise(self, data_index, effective_temperature, surface_gravity,
        metallicity, xi, line_list_filename=None, **theta):

        observed = self.data[data_index]
        expected_flux = np.ones((len(observed.disp)))

        # Synthesise everything except the atomic lines:

        for index in self._transition_mapping[data_index]:
            wavelength = self.atomic_lines[index]["wavelength"]
            element = utils.element(self.atomic_lines[index]["atomic_number"])
            ionised = self.atomic_lines[index]["ionised"] 

            wavelength_region = (
                wavelength - self.synth_surrounding,
                wavelength + self.synth_surrounding
            )
            keyword = self._abundance_key.format(element, ionised, wavelength)
            abundance = theta.get(keyword, metallicity)

            synthesised_spectrum = si.synthesise(effective_temperature,
                surface_gravity, metallicity, xi, wavelength_region,
                abundances={utils.element(species): abundance},
                line_list_filename=line_list_filename)

            # Interpolate the synthesised spectrum onto the observed pixels
            wl_indices = observed.disp.searchsorted(wavelength_region)
            resampled_flux = np.interp(
                observed.disp[wl_indices[0]:wl_indices[1]],
                synthesised_spectrum[:, 0], synthesised_spectrum[:, 1],
                left=1., right=1.)

            # The spectrum enters multiplicatively.
            expected_flux[wl_indices[0]:wl_indices[1]] *= resampled_flux

        return expected_flux


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

            # Synthesise the blended spectrum
            # [TODO] Only do this for the ranges we care about.
            blended_spectrum = si.synthesise(effective_temperature,
                surface_gravity, metallicity, xi, (observed.disp[0], observed.disp[-1]),
                wavelength_steps=(pixel_size, pixel_size, pixel_size),
                line_list_filename=blend_line_list)

            # Convolve the blended spectrum flux
            # [TODO] This should be done by resolution
            if self.config["model"]["doppler_broadening"]:
                sigma = theta["doppler_sigma_{0}".format(i)]
                kernel = (sigma/2.3548200450309493)/pixel_size

                blended_spectrum[:, 1] = gaussian_filter1d(blended_spectrum[:, 1],
                    kernel)

            expected_flux *= np.interp(observed.disp,
                blended_spectrum[:, 0], blended_spectrum[:, 1], 
                left=1., right=1.)

            # Redshift to apply
            z = theta["z"] if "z" in theta else theta.get("z_{0}".format(i), 0)
            for j in self._transition_mapping[i]:

                element = utils.element(self.atomic_lines[j]["atomic_number"])
                ionised = self.atomic_lines[j]["ionised"]
                wavelength = self.atomic_lines[j]["wavelength"]
                abundance = theta[self._abundance_key.format(element, ionised, 
                    wavelength)]
                
                wavelength_region = (
                    wavelength - self.synth_surrounding,
                    wavelength + self.synth_surrounding
                )

                # Synthesise the transition
                # [TODO] Create a singular line list
                transition_spectrum = si.synthesise_transition(effective_temperature,
                    surface_gravity, metallicity, xi, self.atomic_lines[j],
                    surrounding=self.synth_surrounding,
                    abundances={element: abundance},
                    wavelength_steps=(pixel_size, pixel_size, pixel_size))

                # Smoothing
                if self.config["model"]["doppler_broadening"]:
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


    def initial_guess(self, effective_temperature, surface_gravity, metallicity,
        xi, resolution=25000):
        """
        Provide an initial guess of all the SpectrumModel model parameters
        given some stellar parameters.
        """

        initial_guess = {}
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
        if self.config["model"]["doppler_broadening"]:
            for i in range(len(self.data)):
                initial_guess["doppler_sigma_{0}".format(i)] = \
                    np.mean(self.data[i].disp)/resolution

        # Remaining parameters should just be abundances.
        remaining_parameters = set(self.parameters).difference(initial_guess)
        initial_guess.update(dict(zip(remaining_parameters, \
            [metallicity] * len(remaining_parameters))))

        return initial_guess


    def _optimise_transition(self, observed, blending_spectrum, effective_temperature,
        surface_gravity, metallicity, xi, transition_index, full_output=False,
        **theta):
        """
        Fits an absorption profile leaving the redshift + continuum fixed.
        """

        element = utils.element(self.atomic_lines[transition_index]["atomic_number"])
        ionised = self.atomic_lines[transition_index]["ionised"]
        wavelength = self.atomic_lines[transition_index]["wavelength"]
        clean_line_list_filename = self.config["ThereminModel"]["clean_line_list_filename"]

        data_index = self.data.index(observed)
        assert observed.disp[-1] >= wavelength >= observed.disp[0]

        # Get data around that line
        wavelength_region = [
            wavelength - self.synth_surrounding,
            wavelength + self.synth_surrounding
        ]
        indices = observed.disp.searchsorted(wavelength_region)
        indices[-1] += 1 # It's an index thing.
        disp = observed.disp.__getslice__(*indices)
        flux = observed.flux.__getslice__(*indices).copy()
        ivariance = observed.ivariance.__getslice__(*indices)
        
        # Create a mask for this redshift, and apply it to the flux values
        z = theta["z"] if "z" in theta else theta.get("z_{0}".format(data_index), 0)
        flux *= self.mask(disp, z, use_cached=False)
        sigma = theta.get("doppler_sigma_{0}".format(data_index), 0)
        pixel_size = np.mean(np.diff(disp))
        kernel = (sigma/2.3548200450309493)/pixel_size
        
        # Continuum
        if self.config["model"]["continuum"]:
            j, coefficients = 0, []
            while "c_{0}_{1}".format(data_index, j) in theta:
                coefficients.append(theta["c_{0}_{1}".format(data_index, j)])
                j += 1

            if len(coefficients) > 0:
                continuum = np.polyval(coefficients, disp)
        else:
            continuum = 1.

        
        def chi_sq(args, full_output=False):

            abundance = args[0]
            transition_spectrum = si.synthesise_transition(effective_temperature,
                surface_gravity, metallicity, xi, self.atomic_lines[transition_index],
                surrounding=self.synth_surrounding, abundances={element: abundance},
                wavelength_steps=(pixel_size, pixel_size, pixel_size))

            blending_flux = np.interp(transition_spectrum[:, 0],
                blending_spectrum[:, 0], blending_spectrum[:, 1],
                left=np.nan, right=np.nan)
            
            # [TODO] This approximation is only suitable for weak regimes.
            transition_spectrum[:, 1] *= blending_flux

            # Any convolution?
            if kernel > 0:
                transition_spectrum[:, 1] = gaussian_filter1d(
                    transition_spectrum[:, 1], kernel)

            # Put it on the observed pixel space.
            expected = np.interp(disp, transition_spectrum[:, 0] * (1. + z),
                transition_spectrum[:, 1], left=np.nan, right=np.nan)

            # Continuum
            expected *= continuum

            chi_sq = np.nansum((flux - expected)**2 * ivariance)
            if full_output:
                return (chi_sq, expected)
            print(wavelength, args, chi_sq)
            return chi_sq

        # Minimise the function.
        xopt, fopt, direc, niter, nfunc, warnflag = op.fmin_powell(
            chi_sq, [metallicity], xtol=0.01, disp=False, full_output=True)
        self._opt_warn_message(warnflag, niter, nfunc, "fitting {0} {1:.0f} "\
            "transition at {2:.1f} Angstroms".format(element, ionised, wavelength))

        xopt = xopt.reshape(-1)
        
        if full_output:
            return (xopt, fopt, direc, niter, nfunc, warnflag)
        return xopt


    def optimise(self, effective_temperature, surface_gravity, metallicity, xi,
        initial_guess=None, niter=3, use_previously_optimised=True):
        """
        Given some fixed stellar parameters, optimise the model parameters.
        """

        if initial_guess is None:
            if use_previously_optimised and self._optimised_theta is not None:
                logger.info("Using previously optimised theta as initial guess:"\
                    "{0}".format(self._optimised_theta))
                initial_theta = self._optimised_theta

            else:
                initial_theta = self.initial_guess(effective_temperature,
                    surface_gravity, metallicity, xi)

        # Create the final dictionary result
        optimal_theta = {}
        optimal_theta.update(initial_theta)

        # Synthesise all the channels as they are, taking continuum, etc into
        # account
        line_list_filename = self.config["ThereminModel"]["blend_line_list_filename"]

        blending_spectra = [si.synthesise(effective_temperature, surface_gravity,
            metallicity, xi, (s.disp[0], s.disp[-1]), line_list_filename=line_list_filename) \
                for s in self.data]


        # Fit each line/neighbouring individually using the current estimate of
        # redshift, continuum, etc.
        for i, (observed_spectrum, blending_spectrum) \
        in enumerate(zip(self.data, blending_spectra)):

            sigma = initial_theta.get("doppler_sigma_{0}".format(i), 0)
            if sigma > 0:
                # [TODO] do by resolution
                pixel_size = np.mean(np.diff(observed_spectrum.disp))
                kernel = (sigma/2.3548200450309493)/pixel_size
            
                blending_spectrum[:, 1] = gaussian_filter1d(
                    blending_spectrum[:, 1], kernel)

            transition_indices = self._transition_mapping[i]
            for j in transition_indices:

                xopt, fopt, direc, niter, nfunc, warnflag = self._optimise_transition(
                    observed_spectrum, blending_spectrum, effective_temperature,
                    surface_gravity, metallicity, xi, j, full_output=True, **initial_theta)

                wavelength = self.atomic_lines[j]["wavelength"]
                element = utils.element(self.atomic_lines[j]["atomic_number"])
                ionised = self.atomic_lines[j]["ionised"] 

                print("Found {0:.2f} for {1}".format(xopt[0], 
                    self._abundance_key.format(element, ionised, wavelength)))

                optimal_theta[self._abundance_key.format(element, ionised,
                    wavelength)] = xopt[0]


        # Fit everything simultaneously
        initial_fluxes = self(effective_temperature, surface_gravity, metallicity,
            xi, **initial_theta)

        optimal_fluxes = self(effective_temperature, surface_gravity, metallicity,
            xi, **optimal_theta)

        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(self.data))
        for ax, observed, initial, optimal, blended, \
        in zip(axes, self.data, initial_fluxes, optimal_fluxes, blending_spectra):
            ax.plot(observed.disp, observed.flux, 'k')
            ax.plot(observed.disp, initial, 'r', label='initial')
            ax.plot(observed.disp, optimal, 'g', label='final')
            ax.plot(blended[:, 0], blended[:,1], c="#666666", label='blended')

        ax.legend()
        for i, transition in enumerate(self.atomic_lines):
            if i not in sum(self._transition_mapping.values(), []): continue
            
            # Which thing are we in?
            for j, indices in self._transition_mapping.iteritems():
                if i in indices: break

            axes[j].axvline(transition["wavelength"], 0, 1.2, c="b")


        # Fit all lines + redshift + continuum + doppler simultaneously
        raise a


        return None



