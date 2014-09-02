# coding: utf-8

""" A Pythonic Interface to MOOG(SILENT) """

__author__ = "Andy Casey <andy@astrowizici.st>"


# Standard library
import logging
import os
import re
import shutil
from operator import itemgetter
from random import choice
from signal import alarm, signal, SIGALRM, SIGKILL
from string import ascii_letters
from subprocess import PIPE, Popen
from textwrap import dedent

# Third party
import numpy as np
from atmospheres import interpolator as atmospheres

from utils import atomic_number

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

__all__ = ["instance"]


#def equivalent_width(teff, logg, metallicity, xi, rest_wavelengths, species,
#    wavelength_steps=(0.10, 0.005, 1.5), abundances=None, lte=True,
#    full_output=False, threads=1):

class MOOGException(BaseException):
    pass


class instance(object):
    """ A context manager for dealing with MOOG """

    _executable = "MOOGSILENT"
    _acceptable_return_codes = (0, )

    def __init__(self, twd_base_dir="/tmp/", chars=10, debug=False):
        """
        Initialisation class allows the user to specify a base temporary
        working directory.

        :param twd_base_dir: [optional]
            Base directory for temporary working directories.

        :type twd_base_dir:
            str

        :param chars: [optional]
            Number of characters to use for temporary directory names.

        :type chars:
            int

        :param debug: [optional]
            Enable debugging. If MOOG falls over and debug is set to ``True``,
            then the temporary working directory will not be removed.
        """

        self.debug, self.chars, self.twd_base_dir = debug, chars, twd_base_dir
        if not os.path.exists(self.twd_base_dir):
            os.mkdir(self.twd_base_dir)

    def __enter__(self):
        """ Create temporary working directory """

        self.twd = os.path.join(self.twd_base_dir, 
            "".join([choice(ascii_letters) for _ in xrange(self.chars)]))
        while os.path.exists(self.twd):
            self.twd = os.path.join(self.twd_base_dir, 
                "".join([choice(ascii_letters) for _ in xrange(self.chars)]))
        
        os.mkdir(self.twd)
        if len(self.twd) > 40:
            warnings.warn("MOOG has trouble dealing with absolute paths greater"\
                " than 40 characters long. Consider a shorter absolute path for"\
                " your temporary working directory.")
        return self


    def execute(self, filename="batch.par", timeout=30, shell=False, env=None):    
        """
        Execute a MOOG input file with a timeout, after which it will be forcibly
        killed.

        :param filename: [optional]
            The filename to execute with MOOG.

        :type filename:
            str

        :param timeout: [optional]
            The number of seconds to wait before timing out.

        :type timeout:
            int

        :param shell: [optional]
            Execute MOOG in a shell environment.

        :type shell:
            bool

        :param env: [optional]
            Environment variables to pass to the execution environment.

        :type env:
            dict
        """

        logger.info("Executing input file: {0}".format(filename))

        class Alarm(Exception):
            pass

        def alarm_handler(signum, frame):
            raise Alarm

        if env is None and len(os.path.dirname(self._executable)) > 0:
            env = {"PATH": os.path.dirname(self._executable)}

        p = Popen([os.path.basename(self._executable)], shell=shell, bufsize=2056,
            cwd=self.twd, stdin=PIPE, stdout=PIPE, stderr=PIPE, env=env, close_fds=True)

        if timeout != -1:
            signal(SIGALRM, alarm_handler)
            alarm(timeout)

        try:
            pipe_input = "\n" if -6 in self._acceptable_return_codes else ""
            pipe_input += os.path.basename(filename) + "\n"*100

            stdout, stderr = p.communicate(input=pipe_input)
            if timeout != -1:
                alarm(0)

        except Alarm:

            # process might have died before getting to this line
            # so wrap to avoid OSError: no such process
            try:
                os.kill(p.pid, SIGKILL)
            except OSError:
                pass
            return (-9, "", "The process timed out.")

        if p.returncode not in self._acceptable_return_codes:
            logger.warn("MOOG returned the following message (code: {0:d}:".format(
                p.returncode))
            logger.warn(stdout)
            raise MOOGException(stderr)
            
        return (p.returncode, stdout, stderr)


    def _cp_to_twd(self, filename):
        # [TODO] make this symlink instead.
        if os.path.dirname(filename) != self.twd:
            shutil.copy(filename, self.twd)
            filename = os.path.join(self.twd, os.path.basename(filename))

        elif not os.path.exists(filename):
            raise IOError("filename {0} does not exist".format(filename))

        return filename


    def _format_ew_input(self, measurements, comment=None):
        """
        measurments should be recarray
        """
        
        output = comment.rstrip() if comment is not None else ""
        output += "\n"

        line = "{0:10.3f} {1:9.3f} {2:8.2f} {3:6.2f}                             {4:5.1f}\n"

        if np.all(measurements["loggf"] > 0):
            warnings.warn("The atomic line list contains no lines with positive oscillator "
                "strengths. MOOG will not treat these as logarithmic oscillator strengths!")

        # Sort all the lines first transition, then by wavelength
        measurements = sorted(measurements, key=itemgetter("species", "wavelength"))
        for i, measurement in enumerate(measurements):
            # [TODO]: Ignoring van Der Waal damping coefficients for the moment << implement if they exist!
            output += line.format(*[measurement[col] for col in ["wavelength", "species",
                "excitation_potential", "loggf", "equivalent_width"]])
    
        return output
        

    def _format_abfind_input(self, model_atmosphere_filename, line_list_filename,
        standard_out, summary_out, terminal="x11", atmosphere=1, molecules=1, 
        truedamp=1, lines=1, freeform=0, flux_int=0, damping=0, units=0):

        output = """
        abfind
        terminal '{terminal}'
        standard_out '{standard_out}'
        summary_out '{summary_out}'
        lines_in '{model_atmosphere_filename}'
        model_in '{line_list_filename}'
        atmosphere {atmosphere}
        molecules {molecules}
        lines {lines}
        freeform {freeform}
        flux/int {flux_int}
        damping {damping}
        plot 0
        """.format(**locals())
        return dedent(output).lstrip()


    def _parse_abfind_summary_output(self, filename):
        """ Reads the summary output filename after MOOG's `abfind` has been
        called and returns a numpy record array """

        with open(filename, "r") as fp:
            output = fp.readlines()

        data = []
        columns = ("wavelength", "species", "excitation_potential", "loggf",
            "equivalent_width", "abundance")

        for i, line in enumerate(output):
            if line.startswith("Abundance Results for Species "):
                element, ionization = line.split()[4:6]
                current_species = atomic_number(element) + 0.1 * (len(ionization) - 1)
                
                # Check if we already had this species. If so then MOOG has run >1 iteration.
                if len(data) > 0:
                    exists = np.where(np.array(data)[:, 1] == current_species)

                    if len(exists[0]) > 0:
                        logger.debug("Detecting more than one iteration from MOOG")
                        data = list(np.delete(np.array(data), exists, axis=0))
                continue

            elif re.match("^   [0-9]", line):
                line_data = map(float, line.split())
                # Delete the logRW column
                del line_data[4]
                # Delete the del_avg column
                del line_data[-1] 

                # Insert a species column
                line_data.insert(1, current_species)
                data.append(line_data)
                continue

        return np.core.records.fromarrays(
            np.array(data).T, names=columns, formats=["f8"] * len(columns))


    def abfind(self, teff, logg, metallicity, xi, line_list_filename, clobber=False, 
        **kwargs):
        """ Call `abfind` in MOOG """

        # Interpolate a model atmosphere
        atmosphere_filename = os.path.join(self.twd, "model")
        atmospheres.interpolate(atmosphere_filename, teff, logg, metallicity, +0.4, xi,
            clobber=clobber)

        # Copy line list filename to temporary working directory
        line_list_filename = self._cp_to_twd(line_list_filename)

        # Prepare the input and output filenames
        input_filename, standard_out, summary_out = [os.path.join(self.twd, filename) \
            for filename in ("batch.par", "abfind.std", "abfind.sum")]
        
        # Write the abfind file
        with open(input_filename, "w") as fp:
            fp.write(self._format_abfind_input(line_list_filename,
                atmosphere_filename, standard_out, summary_out, **kwargs))

        # Execute MOOG
        result, stdout, stderr = self.execute()
        return self._parse_abfind_summary_output(summary_out)


    def __exit__(self, exit_type, value, traceback):
        # Remove the temporary working directory and any files in it
        if exit_type not in (IOError, MOOGException) and not self.debug:
            shutil.rmtree(self.twd)

        else:
            logger.info("Temporary directory {0} has been kept to allow debugging".format(self.twd))
        return False


