# coding: utf-8

from __future__ import absolute_import, division, print_function

""" Theremin Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import logging
import os
from time import time

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
        if self.config["model"]["instrumental_broadening"]:
            for i in range(len(self.data)):
                initial_guess["instrumental_resolution_{0}".format(i)] = resolution

        # Remaining parameters should just be abundances.
        remaining_parameters = set(self.parameters).difference(initial_guess)
        initial_guess.update(dict(zip(remaining_parameters, \
            [metallicity] * len(remaining_parameters))))

        return initial_guess


    def _optimise_abundance(self, data, transition, effective_temperature,
        surface_gravity, metallicity, xi, xtol=0.01, free_broadening=False,
        full_output=False, **theta):
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

        def chi_sq(args, full_output=False):

            abundance = args[0]
            transition_spectrum, equivalent_width, stdout = si.synthesise_transition(
                effective_temperature, surface_gravity, metallicity, xi, 
                tuple(transition), abundance, self.synth_surrounding,
                (pixel_size, pixel_size, pixel_size),
                full_output=True)

            # Any convolution?
            if not free_broadening:
                resolution = theta.get("instrumental_resolution_{0}".format(data_index), 0)
                kernel = (wavelength/resolution)/(2.3548200450309493*pixel_size)

                if kernel > 0:
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

            line_range = [
                transition_spectrum[transition_spectrum[:, 1] < 1., 0].min(),
                transition_spectrum[transition_spectrum[:, 1] < 1., 0].max()
            ]
            indices = disp.searchsorted(line_range)
            line_chi_sq_i = chi_sq_i[indices[0]:indices[1]]
            line_chi_sq = line_chi_sq_i[np.isfinite(line_chi_sq_i)]
            r_line_chi_sq = line_chi_sq.sum()/(len(line_chi_sq) - 1)

            
            if full_output:
                return (chi_sq, r_line_chi_sq, equivalent_width, np.vstack([disp, expected, background_flux]).T)
            print(wavelength, args, chi_sq, r_line_chi_sq)
            return chi_sq

        # Minimise the function.
        p0 = [metallicity]
        if free_broadening:
            p0.append(25000)
        xopt, fopt, direc, niter, nfunc, warnflag = op.fmin_powell(
            chi_sq, p0, xtol=xtol, disp=False, full_output=True)
        self._opt_warn_message(warnflag, niter, nfunc, "fitting {0} {1:.0f} "\
            "transition at {2:.1f} Angstroms".format(element, ionised, wavelength))

        xopt = xopt.reshape(-1)
        
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(disp, flux, "-.", c="k")
        ax.plot(disp,flux * self.mask(disp, z, use_cached=False), 'k')
        chi, r_chi, eqw, bar = chi_sq(xopt, True)
        ax.set_title("$\chi^2 = {0:.2f}$".format(r_chi))
        ax.plot(disp, bar[:,1], 'b')
        ax.plot(disp,bar[:,2], ":", c="b")
        ax.set_xlim(wavelength_region)
        ax.axvline(wavelength, c="g")
        fig.savefig("fit-{0:.2f}.png".format(transition["wavelength"]))
        #if wavelength > 4732:
        #    raise a
        plt.close("all")

        result = [xopt]
        result.extend(chi_sq(xopt, True)[1:])

        if full_output:
            return (result, fopt, direc, niter, nfunc, warnflag)
        return result


    def _format_abundance_key(self, transition):

        wavelength = transition["wavelength"]
        element = utils.element(transition["atomic_number"])
        ionised = transition["ionised"] 

        return self._abundance_key.format(element, ionised, wavelength)


    def optimise(self, effective_temperature, surface_gravity, metallicity, xi,
        initial_guess=None, free_broadening=True, use_previously_optimised=True):
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

        :param use_previously_optimised: [optional]
            If this SpectrumModel class has been optimised previously, use the
            previously optimised values as an initial guess. This parameter is
            ignored if ``initial_guess`` is anything but ``None``.

        :type use_previously_optimised:
            bool
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

        # Fit each line/neighbouring individually using the current estimate of
        # redshift, continuum, etc.
        t_init = time()

        line_results = []
        wavelengths = []
        for i, observed_spectrum in enumerate(self.data):

            for transition in (self.atomic_lines[j] for j in self._transition_mapping[i]):

                result, fopt, direc, niter, nfunc, warnflag = self._optimise_abundance(
                    observed_spectrum, transition, effective_temperature,
                    surface_gravity, metallicity, xi, full_output=True, 
                    free_broadening=free_broadening, **initial_theta)

                # result contains:
                # (optimised_values, equivalent_width, spectrum)
                xopt, r_chi_sq, equivalent_width, spectrum = result
            
                wavelengths.append(transition["wavelength"])
                line_results.append([xopt, r_chi_sq, equivalent_width])

                theta_key = self._format_abundance_key(transition)
                optimal_theta[theta_key] = xopt[0]

                print("Found {0:.2f} for {1}".format(result[0][0], theta_key))

        print("Finished the line fitting in {0:.2f} seconds".format(time() - t_init))

        # Fit everything simultaneously
        initial_fluxes = self(effective_temperature, surface_gravity, metallicity,
            xi, **initial_theta)

        optimal_fluxes = self(effective_temperature, surface_gravity, metallicity,
            xi, **optimal_theta)

        def chi_sq(theta, full_output=False):

            theta_dict = dict(zip(self.parameters, theta))
            expected = self(effective_temperature, surface_gravity, metallicity,
                xi, **theta_dict)

            chi_sq = 0
            for o, e in zip(self.data, expected):

                # TODO
                mask = self.mask(o.disp, 0., use_cached=True)
                difference = (o.flux - e)**2 * o.ivariance * mask
                chi_sq += difference[np.isfinite(difference)].sum()

            print(theta, chi_sq)
            if full_output:
                return (chi_sq, expected)
            return chi_sq

        #result, fopt, direc, niter, nfunc, warnflag = op.fmin_powell(
        #    chi_sq, [optimal_theta.get(p) for p in self.parameters])
        
        #final_theta = dict(zip(self.parameters, result))
        #final_fluxes = self(effective_temperature, surface_gravity, metallicity,
        #    xi, **final_theta)


        import matplotlib.pyplot as plt

        if free_broadening:
            fig, ax = plt.subplots()
            ax.scatter(wavelengths, [each[0][1] for each in line_results])
            ax.set_xlabel("wavelength")
            ax.set_ylabel("resolution")
            fig.savefig("resolutions.png")

        transition_indices = np.array(sum([self._transition_mapping[i] for i in xrange(len(self.data))], []))

        fig, axes = plt.subplots(2)
        ordered_transitions = self.atomic_lines[transition_indices]
        wavelengths = np.array(wavelengths)
        abundances = np.array([each[0][0] for each in line_results])
        chi_sqs = np.array([each[1] for each in line_results])
        equivalent_widths = np.array([each[2] for each in line_results])

        ok = (equivalent_widths > 0)

        scat = axes[0].scatter(ordered_transitions["excitation_potential"][ok], abundances[ok], c=chi_sqs[ok], cmap='YlOrRd')
        scat2 = axes[1].scatter(np.log(equivalent_widths/wavelengths)[ok], abundances[ok], c=chi_sqs[ok], cmap='YlOrRd')
        cbar = plt.colorbar(scat, cmap='YlOrRd')
        cbar.set_label("$\chi^2$")
        fig.savefig("excitation.png")


        fig, axes = plt.subplots(len(self.data))
        for i, (ax, observed, initial, optimal) \
        in enumerate(zip(axes, self.data, initial_fluxes, optimal_fluxes)):
            ax.plot(observed.disp, observed.flux, 'k')
            ax.plot(observed.disp, initial, 'r', label='initial')
            ax.plot(observed.disp, optimal, 'g', label='final')
            #ax.plot(observed.disp, final, 'b', label='real final')
            #ax.plot(blended[:, 0], blended[:,1], c="#666666", label='blended')

            #[ax.plot(each[:,0], each[:,1], 'g') for each in opt_spectra[i]]
            ax.set_xlim(observed.disp[0], observed.disp[-1])

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



