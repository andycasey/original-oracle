# coding: utf-8

""" Interface with Spectrum Investigation """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import multiprocessing as mp
import numpy as np
import os
import re
import shutil
from random import choice
from signal import alarm, signal, SIGALRM, SIGKILL
from string import ascii_letters
from subprocess import PIPE, Popen
from textwrap import dedent

from utils import atomic_number

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("si")

class SIException(Exception):
    pass

class instance(object):
    
    # It is assumed that the setup.py installer has placed si_lineform on your $PATH
    _executable = "si_lineform"
    _acceptable_return_codes = (0, )

    def __init__(self, twd_base_dir="/tmp", prefix="si", chars=10, debug=False):
        """
        A context manager for interfacing with SI.

        :param twd_base_dir: [optional]
            Base directory for temporary working directory required by SI.

        :type twd_base_dir: str

        :param prefix: [optional]
            Filename prefix to use for temporary files.

        :type prefix: str

        :param chars: [optional]
            Number of random characters to use for temporary file names.

        :type chars: int

        :param debug: [optional]
            If SI raises an exception, keep temporary files and re-raise the
            exception.

        :type debug: bool
        """

        if prefix is None:
            prefix = ""

        self.debug = debug
        self.chars = chars
        self.prefix = prefix
        self.twd_base_dir = twd_base_dir
        if not os.path.exists(self.twd_base_dir):
            os.mkdir(self.twd_base_dir)


    def __enter__(self):
        # Create a temporary working directory.
        self.twd = os.path.join(self.twd_base_dir, self.prefix \
            + "".join([choice(ascii_letters) for _ in xrange(self.chars)]))
        while os.path.exists(self.twd):
            self.twd = os.path.join(self.twd_base_dir, self.prefix \
                + "".join([choice(ascii_letters) for _ in xrange(self.chars)]))
        os.mkdir(self.twd)
        logging.debug("Temporary working directory: {0}".format(self.twd))

        # Link the line list.
        siu_line_list = os.path.abspath(os.path.join(os.path.dirname(os.path.expanduser(__file__)), 
            "../si/linedata/master_line.dat")) 
        os.symlink(siu_line_list, os.path.join(self.twd, "linedata.dat"))

        return self


    def execute(self, filename="fort.10", timeout=60, shell=False, env=None):
        """
        Execute SI with the input filename.
        
        :param filename: [optional]
            The input filename to execute with SI.

        :type filename: str

        :param timeout: [optional]
            Number of seconds to wait for output before timing out.

        :type timeout: int

        :param shell: [optional]
            Execute the command in shell.

        :type shell: bool

        :param env: [optional]
            Enviroment variables to pass through to the SI environment.

        :type env: dict
        """

        logger.debug("Executing input file: {0}".format(
            os.path.join(self.twd, filename)))

        class Alarm(Exception):
            pass

        def alarm_handler(signum, frame):
            raise Alarm

        default_env = {}
        default_env.update(os.environ.copy())
        # We will make this relative to __file__ until we can remove
        # all instances of SIU_MAIN in the fortran code itself.
        default_env["SIU_MAIN"] = os.path.abspath(
            os.path.join(os.path.dirname(os.path.expanduser(__file__)), "../si"))
        
        if env is not None:
            default_env.update(env)

        p = Popen([self._executable],
            shell=shell, bufsize=2056, cwd=self.twd, stdin=PIPE, stdout=PIPE,
            stderr=PIPE, env=default_env, close_fds=True)

        if timeout != -1:
            signal(SIGALRM, alarm_handler)
            alarm(timeout)

        try:
            pipe_input = "\n" if -6 in self._acceptable_return_codes else ""
            
            stdout, stderr = p.communicate(input=pipe_input)
            logger.debug("SI stdout (code {0}): {1}".format(p.returncode, stdout))
            if timeout != -1:
                alarm(0)

        except Alarm:
            # Process might have died before getting to this line so wrap it to
            # avoid "OSError: no such process"
            try:
                os.kill(p.pid, SIGKILL)
            except OSError:
                pass
            return (-9, "", "The process was killed due to timeout.")

        if p.returncode not in self._acceptable_return_codes:
            logger.warn("SI returned the following error (code {0:d}):".format(
                p.returncode))
            logger.warn(stdout)

            raise SIException(stderr)

        return (p.returncode, stdout, stderr)


    def equivalent_width(self, teff, logg, metallicity, xi, rest_wavelength,
        species, wavelength_steps=(0.10, 0.005, 1.5), abundances=None, lte=True,
        full_output=False):
        """
        Return an equivalent width for a transition from ``species`` that occurs
        near or at ``rest_wavelength`` for the given stellar parameters and
        abundances.

        :param teff:
            Effective temperature of the model atmosphere (Kelvin).

        :type teff:
            float

        :param logg:
            Surface gravity of the model atmosphere.

        :type logg:
            float

        :param metallicity:
            Mean metallicity of the model atmosphere.

        :type metallicity:
            float

        :param xi:
            Microturbulence in the model atmosphere (km/s).

        :type xi:
            float

        :param rest_wavelength:
            The rest wavelength of the transition.

        :type rest_wavelength:
            float

        :param species:
            The atomic number and ionisation state of the transition (e.g., for
            Mg I use 12.01).
    
        :type species:
            float
    
        :param wavelength_steps: [optional]
            The optimal, minimum and maximum wavelength step to synthesise
            (Angstroms).

        :type wavelength_steps:
            tuple

        :param abundances: [optional]
            The abundances (values) of different elements (keys).

        :type abundances:
            dict

        :param lte: [optional]
            Employ the approximation of local thermodynamic equilibrium.

        :type lte:
            bool

        :param full_output: [optional]
            Return the synthesised spectra and the SI output in addition to the
            equivalent width. If ``True``, the return format is ``equivalent_width``,
            ``synthetic_spectra``, ``stdout``.

        :type full_output:
            bool

        :returns:
            The integrated equivalent width of the transition.

        :rtype: 
            float
        """

        # Remove any old output files so that we don't accidentally confuse
        # ourselves.
        filenames = ("fort.14", "fort.16")
        for filename in filenames:
            try:
                os.remove(os.path.join(self.twd, filename))
            except OSError:
                continue

        # Write the new input file/
        with open(os.path.join(self.twd, "fort.10"), "w+") as fp:

            # SI requires the minimum step width to be in mA, and the maximum
            # step width to be in Angstroms. Because that makes perfect sense.
            fp.write(dedent("""
            {teff:.0f}
            {logg:.3f}
            {metallicity:+.3f}
            {xi:.3f}
            {rest_wavelength:.3f} {species:.3f}
            {min_wl_step:.1f} {max_wl_step:.3f}
            {si_bit:.0f}
            {opt_wl_step:.3f}
            """.format(teff=teff, logg=logg, metallicity=metallicity, xi=xi,
                rest_wavelength=rest_wavelength, species=species,
                min_wl_step=wavelength_steps[1] * 1000.,
                max_wl_step=wavelength_steps[2],
                opt_wl_step=wavelength_steps[0],
                si_bit=8421423 if lte else 26247215)).lstrip())

            if abundances is None:
                abundances = {}

            fp.write("{0:.0f}\n".format(len(abundances)))
            for element, abundance in abundances.iteritems():
                fp.write("{atomic_number:.0f} {abundance:.3f}\n".format(
                    atomic_number=atomic_number(element),
                    abundance=abundance))
        # Execute it.
        returncode, stdout, stderr = self.execute()

        try:
            equivalent_width = np.loadtxt(os.path.join(self.twd, "fort.16"),
                usecols=(4, )).flatten()[0]

        except IOError:
            logger.warn("No equivalent width found in {0}:".format(
                os.path.join(self.twd, "fort.16")))
            logger.warn("SI output (code {0}):\n{1}\n{2}".format(
                returncode, stdout, stderr))
            raise SIException("no equivalent width found in {0}".format(
                os.path.join(self.twd, "fort.16")))

        if full_output:
            # Grab the synthetic spectra too.
            try:
                synthetic_spectra = np.loadtxt(os.path.join(self.twd, "fort.14"))

            except IOError:
                logger.warn("No synthetic spectra found in {0}:".format(
                    os.path.join(self.twd, "fort.14")))
                logger.warn("SI output (code {0}):\n{1}\n{2}".format(
                    returncode, stdout, stderr))
                raise SIException("no synthetic spectra found in {0}".format(
                    os.path.join(self.twd, "fort.14")))

            if len(synthetic_spectra) == 0:
                logger.warn("No synthetic spectra found in {0}:".format(
                    os.path.join(self.twd, "fort.14")))
                logger.warn("SI output (code {0}):\n{1}\n{2}".format(
                    returncode, stdout, stderr))
                raise SIException("no synthetic spectra found in {0}".format(
                    os.path.join(self.twd, "fort.14")))

            return (equivalent_width, synthetic_spectra, stdout)
        return equivalent_width


    def synthesise(self, teff, logg, metallicity, xi, wavelengths,
        wavelength_steps=(0.10, 0.005, 1.5), abundances=None, full_output=False):
        """
        Synthesise spectra for a given effective temperature, surface gravity,
        metallicity, microturbulence, and abundances.

        :param teff:
            Effective temperature of the model atmosphere (Kelvin).

        :type teff:
            float

        :param logg:
            Surface gravity of the model atmosphere.

        :type logg:
            float

        :param metallicity:
            Mean metallicity of the model atmosphere.

        :type metallicity:
            float

        :param xi:
            Microturbulence in the model atmosphere (km/s).

        :type xi:
            float

        :param wavelength_steps: [optional]
            The optimal, minimum and maximum wavelength step to synthesise
            (Angstroms).

        :type wavelength_steps:
            tuple

        :param abundances: [optional]
            The abundances (values) of different elements (keys).

        :type abundances:
            dict

        :param full_output: [optional]
            Return the SI output as well as the synthesised spectra.

        :type full_output:
            bool

        :raises:
            SIException

        :returns:
            The synthetic spectra as an array of wavelengths and intensities.

        :rtype: :class:`numpy.array`
        """

        # Remove any old output files so that we don't accidentally confuse
        # ourselves.
        try:
            os.remove(os.path.join(self.twd, "fort.14"))
        except OSError:
            pass

        # Write the new input file/
        with open(os.path.join(self.twd, "fort.10"), "w+") as fp:

            # SI requires the minimum step width to be in mA, and the maximum
            # step width to be in Angstroms. Because that makes perfect sense.
            fp.write(dedent("""
            {teff:.0f}
            {logg:.3f}
            {metallicity:+.3f}
            {xi:.3f}
            {wl_start:.2f} {wl_end:.2f}
            {min_wl_step:.2f} {max_wl_step:.4f}
            32811
            {opt_wl_step:.4f}
            """.format(teff=teff, logg=logg, metallicity=metallicity, xi=xi,
                wl_start=min(wavelengths), wl_end=max(wavelengths),
                min_wl_step=wavelength_steps[1] * 1000.,
                max_wl_step=wavelength_steps[2],
                opt_wl_step=wavelength_steps[0])).lstrip())

            if wavelengths[0] > wavelengths[1]:
                logging.warn("Synthesis region requested is out of order:\
                    {0:.0f} > {1:.0f}".format(wavelengths[0], wavelengths[1]))

            if abundances is None:
                abundances = {}

            fp.write("{0:.0f}\n".format(len(abundances)))
            for element, abundance in abundances.iteritems():
                fp.write("{atomic_number:.0f} {abundance:.3f}\n".format(
                    atomic_number=atomic_number(element),
                    abundance=abundance))
        # Execute it.
        returncode, stdout, stderr = self.execute()

        try:
            spectrum = np.loadtxt(os.path.join(self.twd, "fort.14"))

        except IOError:
            logger.warn("No synthetic spectra found in {0}:".format(
                os.path.join(self.twd, "fort.14")))
            logger.warn("SI output (code {0}):\n{1}\n{2}".format(
                returncode, stdout, stderr))
            raise SIException("no fluxes synthesised by SI in {0}".format(
                os.path.join(self.twd, "fort.14")))

        if len(spectrum) == 0:
            logger.warn("No synthetic spectra found in {0}:".format(
                os.path.join(self.twd, "fort.14")))
            logger.warn("SI output (code {0}):\n{1}\n{2}".format(
                returncode, stdout, stderr))
            raise SIException("no synthetic spectra found in {0}".format(
                    os.path.join(self.twd, "fort.14")))

        if full_output:
            return (spectrum, stdout)
        return spectrum


    def __exit__(self, exit_type, value, traceback):
        # Remove the temporary working directory and any files in it.
        if exit_type not in (IOError, SIException) and not self.debug:
            shutil.rmtree(self.twd)
        else:
            logger.info("Temporary directory {0} has been kept to allow debugging"
                .format(self.twd))
        return False


def _synthesise_wrapper(*args, **kwargs):
    """ Wrapper for synthesiser so we can parallelise it. """
    with instance() as si:
        spectrum = si.synthesise(*args, **kwargs)
    return spectrum


def _equivalent_width_wrapper(*args, **kwargs):
    """ Wrapper for the equivalent width function so we can parallelise it. """
    with instance() as si:
        equivalent_width = si.equivalent_width(*args, **kwargs)
    return equivalent_width


def equivalent_width(teff, logg, metallicity, xi, rest_wavelengths, species,
    wavelength_steps=(0.10, 0.005, 1.5), abundances=None, lte=True,
    full_output=False, threads=1):
    """
    Return an equivalent width for a transition from ``species`` that occurs
    near or at ``rest_wavelength`` for the given stellar parameters and
    abundances.

    :param teff:
        Effective temperature of the model atmosphere (Kelvin).

    :type teff:
        float

    :param logg:
        Surface gravity of the model atmosphere.

    :type logg:
        float

    :param metallicity:
        Mean metallicity of the model atmosphere.

    :type metallicity:
        float

    :param xi:
        Microturbulence in the model atmosphere (km/s).

    :type xi:
        float

    :param rest_wavelengths:
        The rest wavelength of the transitions.

    :type rest_wavelengths:
        float or list-type of floats

    :param species:
        The atomic number and ionisation state of the transition (e.g., for
        Mg I use 12.01).

    :type species:
        float or list-type of floats

    :param wavelength_steps: [optional]
        The optimal, minimum and maximum wavelength step to synthesise
        (Angstroms).

    :type wavelength_steps:
        tuple

    :param abundances: [optional]
        The abundances (values) of different elements (keys).

    :type abundances:
        dict

    :param lte: [optional]
        Employ the approximation of local thermodynamic equilibrium.

    :type lte:
        bool

    :param full_output: [optional]
        Return the synthesised spectra and the SI output in addition to the
        equivalent width.

    :type full_output:
        bool

    :param threads: [optional]
        The maximum number of threads to use. Specifying -1 will set the number
        of threads as equal to the number of available cores.

    :type threads:
        int

    :returns:
        The integrated equivalent width of the transition.

    :rtype: 
        float or list-type of floats
    """

    # Check if we only have one spectral region to syntheise.
    single_call = isinstance(rest_wavelengths, (int, float))
    if single_call:
        species = [species]
        rest_wavelengths = [rest_wavelengths]
    else:
        assert len(species) == len(rest_wavelengths), """Species and rest
            wavelengths must be of the same length"""

    # Check number of threads to use.
    threads = [threads, mp.cpu_count()][threads < 0]

    equivalent_widths = []
    if threads == 1:
        # Serial killer.
        with instance() as si:
            for specie, rest_wavelength in zip(species, rest_wavelengths):
                equivalent_widths.append(si.equivalent_width(teff, logg,
                    metallicity, xi, rest_wavelength, specie, wavelength_steps,
                    abundances, lte, full_output))
    else:
        
        pool = mp.Pool(threads)
        results = []
        for specie, rest_wavelength in zip(species, rest_wavelengths):
            results.append(pool.apply_async(_equivalent_width_wrapper,
                args=(teff, logg, metallicity, xi, rest_wavelength, specie,
                    wavelength_steps, abundances, lte, full_output)))

        equivalent_widths = [result.get() for result in results]

        # Winter is coming.
        pool.close()
        pool.join()

    return equivalent_widths[0] if single_call else equivalent_widths


def synthesise(teff, logg, metallicity, xi, wavelengths,
    wavelength_steps=(0.10, 0.005, 1.5), abundances=None, full_output=False,
    threads=1, chunk=True):
    """
    Synthesise spectra for a given effective temperature, surface gravity,
    metallicity, microturbulence, and abundances.

    :param teff:
        Effective temperature of the model atmosphere (Kelvin).

    :type teff:
        float

    :param logg:
        Surface gravity of the model atmosphere.

    :type logg:
        float

    :param metallicity:
        Mean metallicity of the model atmosphere.

    :type metallicity:
        float

    :param xi:
        Microturbulence in the model atmosphere (km/s).

    :type xi:
        float

    :param wavelength_steps: [optional]
        The optimal, minimum and maximum wavelength step to synthesise
        (Angstroms).

    :type wavelength_step_limits:
        tuple

    :param abundances: [optional]
        The abundances (values) of different elements (keys).

    :type abundances:
        dict

    :param full_output: [optional]
        Return the SI output as well as the synthesised spectra.

    :type full_output:
        bool

    :param threads: [optional]
        The maximum number of parallel threads to use for synthesis.

    :type threads:
        int

    :param chunk: [optional]
        Chunk the spectra into smaller portions for parallelisation.

    :type chunk:
        bool

    :raises:
        SIException

    :returns:
        The synthetic spectra as an array of wavelengths and intensities.

    :rtype: :class:`numpy.array`
    """

    # Check if we only have one spectral region to syntheise.
    single_spectrum = (len(wavelengths) == 2 and isinstance(wavelengths[0], (int, float)))
    if single_spectrum:
        wavelengths = [wavelengths]
        wavelength_steps = [wavelength_steps]

    else:
        # Check if we have been given many wavelength ranges but only one set of
        # wavelength steps.
        if isinstance(wavelength_steps[0], (int, float)):
            wavelength_steps = [wavelength_steps] * len(wavelengths)
    
    # Check number of threads to use.
    threads = [threads, mp.cpu_count()][threads < 0]

    # Synthesise all wavelength regions.
    spectra = []
    stdouts = []
    if threads == 1:
        # Serial killer.
        with instance() as si:
            for wavelength_range, wavelength_step in zip(wavelengths, wavelength_steps):
                spectra.append(si.synthesise(teff, logg, metallicity, xi,
                    wavelength_range, wavelength_step, abundances, False))
    else:
        # Thread dat shit.
        if chunk:
            # How many times should every channel be split?
            # (Since we pool it, best to over-chunk than under-chunk)
            num_chunks = int(np.ceil(float(threads)/len(wavelengths)))
            chunks = []
            chunked_wavelength_steps = []
            for wavelength_range, wavelength_step in zip(wavelengths, wavelength_steps):
                spacing = np.linspace(wavelength_range[0], wavelength_range[1],
                    num_chunks + 1)
                chunks.extend([spacing[j:j+2] for j in range(num_chunks)])
                chunked_wavelength_steps.extend([wavelength_step] * num_chunks)

            wavelengths = chunks
            wavelength_steps = chunked_wavelength_steps

        pool = mp.Pool(threads)
        results = []
        for wavelength_range, wavelength_step in zip(wavelengths, wavelength_steps):
            results.append(pool.apply_async(_synthesise_wrapper, args=(teff, logg,
                metallicity, xi, wavelength_range, wavelength_step, abundances,
                True)))

        # Join the chunks back together. 
        chunked_spectrum = []
        for i, result in enumerate(results):
            spectrum, stdout = result.get()
            stdouts.append(stdout)

            if len(spectrum) == 0:
                raise SIException("no fluxes found in spectrum chunk")

            if chunk:
                # Ignore last pixel of the bluer spectrum, except for the last one.
                chunked_spectrum.append(spectrum[:-1, :] if len(results) > i + 1 \
                    else spectrum)
                if (i + 1) % num_chunks == 0:
                    spectra.append(np.vstack(chunked_spectrum))
                    chunked_spectrum = []
            else:
                spectra.append(spectrum)

        # Winter is coming.
        pool.close()
        pool.join()

    if full_output:
        return (spectra[0], stdouts[0]) if single_spectrum else (spectra, stdouts)
    return spectra[0] if single_spectrum else spectra
    
