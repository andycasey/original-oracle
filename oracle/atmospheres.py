# coding: utf-8

""" Stellar atmosphere model interpolator """

__author__ = "Andy Casey <andy@astrowizici.st>"

# Standard libraries
import os
import logging
import re

from collections import OrderedDict
from glob import glob
from gzip import open as gzip_open
from textwrap import dedent

# Third party imports
import numpy as np
import scipy.interpolate
from scipy.io import readsav

# Module imports
from utils import atomic_number

__all__ = [
    "AtmosphereInterpolator", "AtmosphereParser", "CastelliKuruczParser", \
    "AtmoATLASParser", "MARCSParser", "StaggerGridATLASParser", "parsers"
]

logger = logging.getLogger("oracle")

class AtmosphereParser(object):
    """Base class for all AtmosphereParser objects.
    
    Notes
    ----
    Any sub-classes derived from this `AtmosphereParsers` class should have at
    least two methods: `parse_thermal_structure`, and `parse_filename` where 
    each only takes a single input - the model atmosphere filename.
    
    The `parse_thermal_structure` method should interpret the atmosphere filename
    and return an `np.array` containing the deck values.
    
    The `parse_filename` method should interpret the filename path and return a
    list-type containing the temperature, surface gravity, metallicity, and alpha
    enhancement.
    
    Lastly, the class should also hold an attribute of `basename_match`, which
    will be a string match to identify model atmosphere files (e.g. '*.dat')."""
    
    # By default, just do linear interpolation of properties
    # This is a legacy thing.
    logarithmic_structure_properties = None
    structure_properties_to_interpolate = "all"
    pass


class AtmoATLASParser(AtmosphereParser):
    """A class capable of interpreting ATMO 3D stellar atmosphere files that have
    been exported in an ATLAS-style format for use with the `AtmosphereInterpolator`
    class."""

    basename_match = 'atlas.*_atmo'

    @staticmethod 
    def parse_filename(filename):
        """Parses the standard filename for an ATMO 3D model atmosphere filename.

        [teff, logg, feh, alpha]
        """

        filename = os.path.basename(filename)[7:-5]
        # XXXXgYY[m|p]ZZ where X = teff, Y = logg, Z = feh

        teff = float(filename.split('g')[0])
        logg = float(filename.split('g')[1].replace('m', 'p').split('p')[0]) /10.
        feh = float(filename[-3:].replace('m', '-').replace('p', '+')) / 10.

        alpha = 0.4

        return [teff, logg, feh, alpha]


    @staticmethod
    def parse_thermal_structure(filename):
        """Reads in the thermal structure of a <3D> Stagger-grid model that has
        been exported into 'ATLAS-style' format and returns an array of the 
        thermal structure."""

        with open(filename, 'r') as fp:
            contents = fp.readlines()
        
        in_deck, deck = False, []
        
        for line in contents:
            if line.startswith('READ DECK'):
                in_deck = True
                continue
            elif line.startswith('PRADK'):
                break
            
            if in_deck:

                # First two columns have peculiar lengths
                values = [line[:15], line[15:24]]
                cut_line = line[24:]
                
                # After that it's all 10 character columns.
                values.extend([cut_line[i*10:(i + 1)*10] for i in xrange(len(cut_line)/10)])

                deck.append(map(float, values))
            
        return np.array(deck)

    @staticmethod
    def write_atmosphere(*args):
        raise NotImplementedError


class StaggerGridATLASParser(AtmosphereParser):
    """A class capable of interpreting Stagger-grid 3D stellar amostphere files
    that have been exported in an ATLAS-style format for use with the 
    `AtmosphereInterpolator` class."""

    basename_match = 'atlas.t*i'

    @staticmethod
    def parse_filename(filename):
        """Parses the standard filename for a Stagger-grid 3D model atmosphere filename.

        [teff, logg, feh, alpha]
        """

        filename = os.path.basename(filename)[7:-1]
        # XXXXgYY[m|p]ZZ where X = teff, Y = logg, Z = feh

        teff = float(filename.split('g')[0])
        logg = float(filename.split('g')[1].replace('m', 'p').split('p')[0]) / 10.
        feh = float(filename[-3:].replace('m', '-').replace('p', '+')) /10.

        alpha = 0.4

        return [teff, logg, feh, alpha]


    @staticmethod
    def parse_thermal_structure(filename):
        """Reads in the thermal structure of a 3D Stagger-grid model that has
        been exported into 'ATLAS-style' format and returns an array of the 
        thermal structure."""

        with open(filename, 'r') as fp:
            contents = fp.readlines()
        
        in_deck, deck = False, []
        
        for line in contents:
            if line.startswith('READ DECK'):
                in_deck = True
                continue
            elif line.startswith('PRADK'):
                break
            
            if in_deck:

                # First two columns have peculiar lengths
                values = [line[:15], line[15:24]]
                cut_line = line[24:]
                
                # After that it's all 10 character columns.
                values.extend([cut_line[i*10:(i + 1)*10] for i in xrange(len(cut_line)/10)])

                deck.append(map(float, values))
            
        return np.array(deck)

    @staticmethod
    def write_atmosphere(*args):
        raise NotImplementedError


class MARCSParser(AtmosphereParser):
    """A class capable of interpreting the MARCs 2011 model atmosphere
    files for use with the `AtmosphereInterpolator` class. """

    basename_match = "*_*.mod.gz" # Avoid the sun.mod.gz
    logarithmic_structure_properties = ["lgTauR", "lgTau5"]
    structure_properties_to_interpolate = \
        ["k", "lgTauR", "lgTau5", "Depth", "T", "Pe", "Pg", "Prad", "Pturb",
        "KappaRoss", "Density", "Mu", "Vconv", "Fconv/F", "RHOX"]

    @staticmethod
    def parse_thermal_structure(filename, full_output=False):
        """
        Reads in a MARCs model atmosphere and returns the thermal
        structure of the photosphere.

        Parameters
        ----------
        filename : str
            The model atmosphere filename.
        """

        open_fn, open_args = (gzip_open, "r") if filename.endswith(".gz") else (open, "rb")
        with open_fn(filename, open_args) as fp:
            contents = fp.readlines()

        thermal_structure = {}
        num_structures_read, num_depth_points, thermal_structure_starts = 0, np.nan, -1
        for i, line in enumerate(contents):

            if line.endswith("Number of depth points\n"):
                num_depth_points = int(line.strip().split()[0])
                continue

            elif line in ("Model structure\n", "Assorted logarithmic partial pressures\n"):
                num_structures_read, thermal_structure_starts = 0, i + 1
                continue

            if thermal_structure_starts + num_structures_read*(num_depth_points + 1) == i:

                # Header information?
                headers = line.replace(" H I ", "H_I").strip().split()
                for header in headers:
                    thermal_structure[header] = []
                num_structures_read += 1
                continue

            if i > thermal_structure_starts > 0:
                data = map(float, re.sub("E[\+|-]\d{2}", lambda x: x.group(0) + " ", line).replace("******", "  nan  ").strip().split())
                for j, header in enumerate(headers):
                    thermal_structure[header].append(data[j])

        if full_output:
            # Create structured array of all the data
            thermal_structure = np.core.records.fromarrays(thermal_structure.values(),
                names=thermal_structure.keys(), formats=["f8"] * len(thermal_structure))
            return thermal_structure

        else:
            # What data columns should we be interpolating?
            data_columns = thermal_structure.keys()
            # Create structured array of the relevant data
            thermal_structure = np.core.records.fromarrays([thermal_structure[column] for column in data_columns],
                names=data_columns, formats=["f8"] * len(data_columns))
            return thermal_structure


    @staticmethod
    def parse_filename(filename, full_output=False):
        """
        Parses the standard filename for a MARCs atmosphere filename and
        returns the stellar parameters.
        """
        
        filename = os.path.basename(filename).lower()
        #p3600_g+4.0_m0.0_t02_ap_z-0.50_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.krz.gz
        #dimensional_type = "parallel" if filename.startswith("p") else "spherical"

        segments = filename.split("_")
        teff = int(filename[1:5])
        logg = segments[1].lstrip("g")
        feh = segments[5].lstrip("z")
        alpha = segments[6].lstrip("a")
        logg, feh, alpha = map(float, [logg, feh, alpha])

        if full_output:
            raise NotImplementedError

        else:
            return np.array([teff, logg, feh, alpha])


    @staticmethod
    def write_atmosphere(filename, teff, logg, fe_h, xi, thermal_structure, 
        solar_abundances=None, molecules=None, clobber=False, lambda_std=5000):

        if os.path.exists(filename) and not clobber:
            raise IOError("Filename {0} exists and we've been asked not to clobber it.".format(filename))

        num_depth_points = thermal_structure.shape[0]
        contents = dedent(
        """
        KURTYPE
        MARCS:  {teff:5.0f}./    {logg:1.2f}/   {fe_h:1.2f}      mic = {xi:1.4f}                                
                     {num_depth_points:.0f}
        {lambda_std:4.1f}
        """.format(teff=teff, logg=logg, fe_h=fe_h, xi=xi, num_depth_points=num_depth_points, lambda_std=lambda_std)).lstrip()

        # rhox, Teff, Pg, Ne
        # 0.93435814E-02   4380.4 2.935E+02 4.585E+10
        # Note that MOOG will convert Ne --> Pe for us automagically
        for i in xrange(num_depth_points):
            contents += " {0:14.8e} {1:7.1f} {2:9.3e} {3:9.3e}\n".format(
                thermal_structure.RHOX[i], thermal_structure["T"][i], thermal_structure.Pg[i], thermal_structure.Pe[i])

        # Write microturbulence
        contents += "     {0:1.2e}\n".format(xi)
        
        # Abundances
        if solar_abundances is None:
            contents += "NATOMS        0           {0:1.2f}\n".format(fe_h)

        else:
            contents += "NATOMS    {0:4.0f}     {1:1.2f}\n".format(len(solar_abundances), fe_h)

            for i, (element, abundance) in enumerate(solar_abundances.iteritems()):
                contents += "{0:4.0f} {1:3.2f}    ".format(atomic_number(element), abundance)
                if (i > 0 and not i % 8) or (i + 1 == len(solar_abundances.keys())): contents += "\n"

        # Molecules
        if molecules is not None:
            contents += "NMOL      %4i\n" % (len(molecules), )
            for i in xrange(np.ceil(len(molecules) / 6)):
                contents += ''.join(['%10s' % (str(int(molecule)), ) for molecule in molecules[6 * i:6 * (i + 1)]]) + "\n"
        else:
            contents += "NMOL          0\n"

        # Add the trailing newline character
        contents += "\n"

        # Write abundances and molecules
        with open(filename, "w") as fp:
            fp.write(contents)

        return True


class CastelliKuruczParser(AtmosphereParser):
    """A class capable of interpreting Castelli & Kurucz (2003) stellar
    atmosphere files for use with the `AtmosphereInterpolator` class."""
    
    basename_match = '*.dat'
    
    @staticmethod
    def parse_thermal_structure(filename):
        """Reads in a Castell-Kurucz stellar atmosphere model and returns only
        an array of the deck."""
        
        with open(filename, 'r') as fp:
            contents = fp.readlines()
        
        in_deck, deck = False, []
        
        for line in contents:
            if line.startswith('READ DECK'):
                in_deck = True
                continue
            elif line.startswith('PRADK'):
                break
            
            if in_deck: deck.append(map(float, line.split()))
            
        return np.array(deck)
        
    @staticmethod
    def parse_filename(filename):
        """Parses the standard filename for a Castelli-Kurucz atmosphere filename.
        
        [teff, logg, feh, alpha]
        """
        filename = os.path.basename(filename)
        
        feh = float(filename.split('t')[0][1::].replace('m', '-').replace('p', '+').rstrip('a')) /10.
        teff = float(filename.split('t')[1].split('g')[0])
        logg = float(filename.split('g')[1].split('k')[0]) /10.
        if filename[4] == 'a': alpha = 0.4
        else: alpha = 0.0
        
        return [teff, logg, feh, alpha]
    
    @staticmethod
    def write_atmosphere(filename, teff, logg, M_H, vt, thermal_structure, solar_abundances=None,
        molecules=None, clobber=False):
        """
        Writes a KURUCZ-style atmosphere file to disk.

        Parameters
        ----------
        filename : str
            The filename to write the atmosphere file to

        teff : float
            Effective surface temperature of the star.

        logg : float
            Surface gravity of the star.

        M_H : float
            Metallicity of the star, [M/H]

        vt : float
            Level of microturbulence (km/s) 

        thermal_structure : `np.ndarray` or `np.core.records.recarray`
            The thermal structure of the model atmosphere

        abundances : `dict`
            A dictionary containing the solar abundances to employ. This should be in
            the form where the elements are as keys, and the abundances as values. The
            element keys can either be by atomic number (e.g. 26), or a string element
            representation (e.g. 'Fe').

        molecules : `list`-type of `int`s
            Molecules to include during the molecular equilibrium.
        """

        if os.path.exists(filename) and not clobber:
            raise IOError("Output filename {0} already exists and we've been asked not to clobber it.".format(filename))

        if not isinstance(teff, (float, int)):
            raise TypeError("Temperature must be a floating point or integer-type.")

        if not isinstance(logg, (float, int)):
            raise TypeError("Surface gravity must be a floating point or integer-type.")

        if not isinstance(M_H, (float, int)):
            raise TypeError("Model metallicity must be a floating point or integer-type.")

        if not isinstance(vt, (float, int)) or vt < 0:
            raise TypeError("Microturbulence must be a positive floating point or integer-type.")

        if solar_abundances is not None:
            if not isinstance(solar_abundances, dict):
                raise TypeError("Abundances must be a dictionary-type.")

            try:
                map(float, solar_abundances.values())

            except TypeError:
                raise TypeError("Abundances must be floating point-types.")

            # Check that the keys are floats already
            try:
                map(float, solar_abundances.values())

            except TypeError:
                
                abundances_with_atomic_numbers = {}
                for element, abundance in solar_abundances.iteritems():
                    abundances_with_atomic_numbers[atomic_number(element)] = abundance

                # Update the normal abundances dictionary
                solar_abundances = abundances_with_atomic_numbers

        if molecules is not None:
            if not isinstance(molecules, (list, tuple)):
                raise TypeError("Molecules must be a list-type.")

            try:
                map(float, molecules)

            except TypeError:
                raise TypeError("Molecules must be floating point-types.")

        output_string = """
        KURUCZ
                  TEFF   %i.  GRAVITY %2.5f LTE
        NTAU        %i
        """ % (teff, logg, len(thermal_structure), )
        
        output_string = dedent(output_string).lstrip()
        
        for line in thermal_structure:
            # Need to accompany thermal structures of different array sizes
            output_string += " %1.8e " % (line[0], )
            output_string += ("%10.3e" * (len(line[1:])) % tuple(line[1:]) + "\n")

        output_string += "         %1.2f\n" % (vt, )

        if solar_abundances is None:
            output_string += "NATOMS        0     %2.1f\n" % (M_H,)

        else:
            output_string += "NATOMS    %4i     %2.1f\n" % (len(solar_abundances), M_H, )
            for i, (element, abundance) in enumerate(solar_abundances.iteritems()):
                output_string += "%4i %3.2f    " % (atomic_number(element), abundance )
                if (i > 0 and not i % 8) or (i + 1 == len(solar_abundances.keys())): output_string += "\n"


        if molecules is not None:
            output_string += "NMOL      %4i\n" % (len(molecules), )
            for i in xrange(np.ceil(len(molecules) / 6)):
                output_string += ''.join(['%10s' % (str(int(molecule)), ) for molecule in molecules[6 * i:6 * (i + 1)]]) + "\n"
        else:
            output_string += "NMOL          0\n"

        # Add the trailing newline character
        output_string += "\n"

        # Write it out
        with open(filename, 'w') as fp:
            fp.write(output_string)

        return True



class StaggerGridInterpolator:
    """Interpolates between stellar atmospheres on a column mass scale """

    def __init__(self, path):

        if not os.path.exists(path):
            raise OSError("path '{0}' does not exist".format(path))

        self.grid = readsav(path)["mmd"]

        # Remove the sun from the grid points?
        #self.grid = 


    def zaz_get_intpolmod(self, teff, logg, feh, verbose=True):
        """ Find the cube surrounding a point """

        teff_attribute = "teffaim" if "teffaim" in self.grid.dtype.names else "teff"

        teffs = np.unique(self.grid[teff_attribute])
        loggs = np.unique(self.grid["logg"])
        fehs = np.unique(self.grid["feh"])

        min_teff, max_teff = min(teffs), max(teffs)
        min_logg, max_logg = min(loggs), max(loggs)
        min_feh, max_feh = min(fehs), max(fehs)

        # TODO - should we calculate this?
        d_teff, d_logg, d_feh = 500, 0.5, 0.5

        if teff < min_teff or teff > max_teff:
            logger.warn("Temperature {teff:.0f} is outside grid boundary: {min:.0f} / {max:.0f}"
                .format(teff=teff, min=min_teff, max=max_teff))
        if logg < min_logg or logg > max_logg:
            logger.warn("Surface gravity {logg:.1f} is outside grid boundary: {min:.1f} / {max:.1f}"
                .format(logg=logg, min=min_logg, max=max_logg))
        if feh < min_feh or feh > max_feh:
            logger.warn("Metallicity {feh:.1f} is outside grid boundary: {min:.1f} / {max:.1f}"
                .format(feh=feh, min=min_feh, max=max_feh))

        # Get the cube
        teff_0 = max(teffs[np.where(teffs <= teff)[0]])
        teff_1 = min(teffs[np.where(teffs > teff)[0]])
        logg_0 = max(loggs[np.where(loggs <= logg)[0]])
        logg_1 = min(loggs[np.where(loggs > logg)[0]])
        feh_0 = max(fehs[np.where(fehs <= feh)[0]])
        feh_1 = min(fehs[np.where(fehs > feh)[0]])

        # May need to edit this to do some checking for whether points are actually returned or not
        cube = [
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_0) & (self.grid.feh == feh_0))[0][0], # 000
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_0) & (self.grid.feh == feh_1))[0][0], # 001
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_0) & (self.grid.feh == feh_0))[0][0], # 010
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_0) & (self.grid.feh == feh_1))[0][0], # 011
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_1) & (self.grid.feh == feh_0))[0][0], # 100 
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_1) & (self.grid.feh == feh_1))[0][0], # 101 <-- IS THIS RIGHT?
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_1) & (self.grid.feh == feh_0))[0][0], # 110
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_1) & (self.grid.feh == feh_1))[0][0], # 111
        ]

        cube = [
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_0) & (self.grid.feh == feh_0))[0][0], # 000
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_0) & (self.grid.feh == feh_1))[0][0], # 001
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_1) & (self.grid.feh == feh_0))[0][0], # 010
            np.where((self.grid[teff_attribute] == teff_0) & (self.grid.logg == logg_1) & (self.grid.feh == feh_1))[0][0], # 011
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_0) & (self.grid.feh == feh_0))[0][0], # 100 
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_1) & (self.grid.feh == feh_1))[0][0], # 101 <-- IS THIS RIGHT?
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_1) & (self.grid.feh == feh_0))[0][0], # 110
            np.where((self.grid[teff_attribute] == teff_1) & (self.grid.logg == logg_1) & (self.grid.feh == feh_1))[0][0], # 111
        ]

        if verbose:
            names = dict(zip(
                ["n000", "n001", "n010", "n011", "n100", "n101", "n110", "n111"],
                [self.grid[index].name for index in cube]
                ))

            message_args = {
                "min_teff": int(teff_0), "max_teff": int(teff_1),
                "min_logg": logg_0, "max_logg": logg_1,
                "min_feh": feh_0, "max_feh": feh_1,
            }
            message_args.update(names)
            verbose_message = \
                """
                Stagger-grid Interpolation Diagram:
                --------------------------------------------------------
                           [Fe/H]: [{min_feh}:{max_feh}]
                               |
                          {n011:s} ____ {n111:s}
                         /     :          / |
                   {n001:s} ____ {n101:s} |
                       |       :     |      |
                       |       :     |      |
                       |       :     |      |
                       |  {n010:s} _|__ {n110:s}  _ logg: [{min_logg}:{max_logg}]
                       |/            |  /
                   {n000:s} ____ {n100:s}
                     /
                   Teff: [{min_teff:d}:{max_teff:d}]
                --------------------------------------------------------
                """

            logger.info(dedent(verbose_message.format(**message_args)))

        return self.grid[cube]



    def zaz_sgi(self, teff, logg, feh, attributes=None, verbose=True):

        if not isinstance(teff, (int, float)):
            raise TypeError("temperature should be an integer or float-like object")
        if not isinstance(logg, (int, float)):
            raise TypeError("logg should be an integer or float-like object")
        if not isinstance(feh, (int, float)):
            raise TypeError("feh should be an integer or float-like object")

        # Get neighbouring models
        cube = self.zaz_get_intpolmod(teff, logg, feh)
        n_depth = len(cube[0].ltaur) # ndep -> n_depth

        # Set up the attributes list
        if attributes is None:
            attributes = ["depth", "nltau", "ltaur", "ltau5", "tt", "rho", "ptot", "pgas",
                "pel", "nel", "kapr"] + ["fconv"]
        else:
            attributes = [attribute for attribite in cube[0].dtype.names if len(cube[0][attribute]) == n_depth]

        # Shift zero depth point
        # NOTE: cube VARIABLE EXPECTED TO BE LENGTH 8, LIKE THAT COOL TESSERACT THING FROM THE AVENGERS
        for i in xrange(8):
            if cube[i].depth[0] != 0:
                cube[i].depth -= cube[i].depth[0]

            # Avoid the infinities
            cube[i].depth += 1e-8

        # Check validity of stellar parameter range
        # all these min_X, max_X, d_X, etc are: minX -> min_X
        teff_attribute = "teffaim" if "teffaim" in cube[0].dtype.names else "teff"
        min_teff, max_teff = min(cube[teff_attribute]), max(cube[teff_attribute])
        min_logg, max_logg = min(cube["logg"]), max(cube["logg"])
        min_feh, max_feh = min(cube["feh"]), max(cube["feh"])

        d_teff, d_logg, d_feh = map(np.ptp, [cube[attribute] for attribute in [teff_attribute, "logg", "feh"]])

        xx = (teff - min_teff)/d_teff if min_teff != max_teff else 0
        yy = (logg - min_logg)/d_logg if min_logg != max_logg else 0
        zz = (feh - min_feh)/d_feh if min_feh != max_feh else 0

        metadata = {
            "name": self.zaz_get_sg_name(teff, logg, feh),
            "teff": teff,
            "logg": logg,
            "feh": feh
            }

        # Start interpolating
        n_attributes = len(attributes)
        # value -> attributes

        structure = {}
        for i, attribute in enumerate(attributes):
            do_logarithmically = self.zaz_get_dolog(attribute)

            values = [point[attribute] for point in cube]
            if do_logarithmically:
                values = np.log10(values)

            interpolated_values = self.zaz_trilinear_interpolate(xx, yy, zz, values)

            if do_logarithmically:
                interpolated_values = 10**interpolated_values

            structure[attribute] = interpolated_values

        # Turn structure into structured array
        structure = np.core.records.fromarrays(structure.values(),
            names=structure.keys(), formats=['f8'] * len(structure))

        # Shift zero point
        structure.depth -= 1e-8
        ltau_attribute = "nltau" if "nltau" in structure.dtype.names else "ltaur"
        tau_23rds_index = np.where(structure[ltau_attribute] >= np.log10(2./3))[0]
        if len(tau_23rds_index):
            structure.depth -= structure.depth[tau_23rds_index]

        else:
            logger.warn("Photospheric depth has not been re-scaled to tau^(2/3)!")

        # Add interpolation information to the metadata
        metadata.update({
            "teff0": min_teff, "teff1": max_teff, "dteff": d_teff, "xx": xx,
            "logg0": min_logg, "logg1": max_logg, "dlogg": d_logg, "yy": yy,
            "feh0": min_feh, "feh1": max_feh, "dfeh": d_feh, "zz": zz,
            "cube": cube.name
            })

        return (structure, metadata)


    def generate_atmosphere(self, output, teff, logg, feh, alpha, xi, solar_abundances=None, molecules=None):

        photospheric_structure  = self.zaz_sgi(teff, logg, feh)




    def export_as_marcs(self, output, structure, teff, logg, feh, alpha, xi, solar_abundances=None, molecules=None,
        tau_min=-4.8, tau_max=4.4, dim=50, rescale=True):

        teff_attribute = "teffaim" if "teffaim" in structure else "teff"

        # Re-scale model to dimensions required between tau_min and tau_max
        if rescale:
            # Do rescale
            raise NotImplementedError

        n_depth = len(structure.ltaur)

        # Shift zero point of depth to tau=1
        #structure = self.shift_depth(structure) ## UNKNOWN

        # Re-scale fconv
        fconv = -structure.fconv/(teff**(4*5.6705e-16))

        # Convert to CGS units
        #structure = convert_structure_cgs(structure) ## UNKNOWN (convert_md_cgs)
        structure.fconv = fconv

        # Get values
        xi = structure.xi if "xi" in structure.dtype.names else 1
        lambda_std = structure["lambda"] if "lambda" in structure.dtype.names else 5000.
        alpha = structure.alpha if "alpha" in structure.dtype.names else 1.5
        beta = structure.beta if "beta" in structure.dtype.names else 0.

        pny, y = 0, 0
        mu = np.zeros(n_depth)
        g_rad = np.zeros(n_depth)

        vconv = structure.vconv if "vconv" in structure.dtype.names else structure.uyrms
        abund = structure.abund if "abund" in structure.dtype.names else get_abund(feh) ## UNKNOWN

        atmosphere_contents = dedent(
            """
              TEFF={teff:8.1f} [K],  LG G={logg:6.2f} [cgs],  XI_TURB={xi:6.2f} [KM/S],{name:24s}     
              Depth points={n_depth:3d},  Lambda_std={lambda_std:8.1f} [AA],  ALPHA={alpha:5.2f}, BETA={beta:5.2f}, PNY={pny:5.2f}, Y={y:6.3f}
            """.format(teff=teff, logg=logg, xi=xi, name="foo", n_depth=n_depth, lambda_std=lambda_std, alpha=alpha, beta=beta, pny=pny, y=y))

        #printf,u,strupcase(abund.el),f='("  ",11a-6)'
        #printf,u,abund.abund,f='(11f6.2)'
        atmosphere_contents += "0   K    TAUROSS   TAU({lambda:d})  GEOM. DEPTH     T        PE          PG         PRAD        PTURB     KAPPAROSS     K\n"
        for i in xrange(n_depth):
            atmosphere_contents += "{0:5d}{1:11.3e}{2:11.3e}{3:13.3e}{4:9.0f}{5:12.3e}{6:12.3e}{7:12.3e}{8:12.3e}{9:12.3e}{10:6d}\n".format(
                i+1, structure.ltaur[i], structure.ltau5[i], structure.depth[i], structure.tt[i], structure.pel[i], structure.ptot[i],
                structure.prad[i], structure.pturb[i], structure.kapr[i], i+1)

        atmosphere_contents += "0  K   TAUROSS   DENSITY     MU         CP          CV       ADGRAD       G_RAD     SOUND VEL.  CONV. VEL. FCONV/F   K\n"
        for i in xrange(n_depth):
            atmosphere_contents += "{0:3d}{1:11.3e}{2:11.3e}{3:8.3f}{4:12.3e}{5:12.3e}{6:12.3e}{7:12.3e}{8:12.3e}{9:12.3e}{10:9.5f}{11:4d}\n".format(
                i+1, structure.ltaur[i], structure.rho[i], mu[i], structure.cp[i], structure.cv[i], strcuture.nad[i], g_rad[i], structure.cs[i],
                vconv[i], structure.fconv[i], i+1)

        with open(output, "w") as fp:
            fp.write(atmosphere_contents)

        print(atmosphere_contents)



    def zaz_trilinear_interpolate(self, x, y, z, data):

        ndep = len(data)
        v000 = data[0]
        v001 = data[1]
        v010 = data[2]
        v011 = data[3]
        v100 = data[4]
        v101 = data[5]
        v110 = data[6]
        v111 = data[7]

        v = np.zeros(ndep)
        for i in xrange(ndep):
            v[i] = v000[i] + x*(-v000[i] + v100[i]) + y*(-v000[i] + v010[i]) + z*(-v000[i] + v001[i]) \
                + x*y*(v000[i] - v100[i] - v010[i] + v110[i]) \
                + x*z*(v000[i] - v100[i] - v001[i] + v101[i]) \
                + y*z*(v000[i] - v010[i] - v001[i] + v011[i]) \
                + x*y*z*(-v000[i] + v100[i] + v010[i] + v001[i] - v011[i] - v101[i] - v110[i] + v111[i])

        return v




    def zaz_get_dolog(self, attribute):
        """ Returns whether an attribute should be interpolated logarithmically
        or linearly. """

        attribute = attribute.lower()
        # Check for some non-logarithmic attributes first so they don't get mixed
        # up in our acceptable_variations
        if attribute in ["pp", "u", "ph"]: return False

        depth_attributes = ["depth"]
        density_attributes = ["rho", "nel"]
        pressure_attributes = ["ptot", "pth", "pgas", "prad", "pel", "pturb"]
        opacity_attributes = ["kapr", "kap5", "kapp", "alphar", "alpha5", "alphap"]
        granule_attributes  = ["area", "diam", "perim"]

        logarithmic_attributes = density_attributes + pressure_attributes \
            + opacity_attributes + granule_attributes + depth_attributes

        if attribute in logarithmic_attributes: return True

        acceptable_variations = ["ln{0}", "med_{0}", "{0}u", "{0}d"]
        for variation in acceptable_variations:
            if attribute in [variation.format(attr) for attr in logarithmic_attributes]: return True

        # By this stage, we have no reason to think logarithmic, right?
        return False



    def zaz_get_sg_name(self, teff, logg, feh, alpha=0):
        """ Returns a shorthand name for a model atmosphere given the input stellar parameters """

        # Create a model name from the stellar parameters
        if teff % 100. == 0: teff_str = "{0:d}".format(teff/100)
        elif teff % 10. == 0: teff_str = "{0:d}".format(teff/10)
        else:
            teff_str = "{0:d}".format(teff) if int(teff) == teff else "{0:.2f}".format(teff)
        separator = "m" if feh < 0 else "p"
        alpha_str = "" if alpha == 0 else "a{0:02df}".format(alpha*10)        
        name = "t{teff_str}g{logg:02d}{separator}{feh:02d}{alpha_str}".format(
            teff_str=teff_str, logg=int(10*logg), feh=int(abs(10*feh)), alpha_str=alpha_str, separator=separator)

        return name






class AtmosphereInterpolator:
    """Interpolate between stellar atmospheres.
    
    When initiating an AtmosphereInterpolator class you must provide the models
    directory where all the files are kept, as well as a parser which is used to
    (1) identify what the stellar parameters are for each file and (2) parse the
    deck within the file.
    
    """
    
    def __init__(self, models_folder, parser, logarithmic_columns=[2, 3]):
        """Initialises an class for interpolating atmospheres from existing
        atmosphere model files.
        
        Parameters
        ----
        models_folder : str
            Path containing all the model atmosphere files to use during
            interpolation. The folder does not have to only include model
            atmosphere files; the basename match is derived from the
            `parser`.
        
        parser : `AtmosphereParser` instance
            An `AtmosphereParser` class capable of parsing the atmosphere
            filename nomenclature, and reading in the deck from each file.
            
            See the `AtmosphereParser` class for more details.
        
        logarithmic_columns : list of ints, default
            A list of integers indicating which stellar parameter columns
            are logarithmic values and should be interpolated appropriately.
            In general the columns are referenced as temperature, surface
            gravity, metallicity and alpha enhancement. Thus, column indices
            1, 2, and 3 are all logarithmic values.
        """
        
        if not isinstance(parser, AtmosphereParser):
            raise TypeError("Expected parser to be an AtmosphereParser class.")
        
        self.folder = models_folder
        self.parser = parser
        
        if hasattr(parser, 'basename_match'):
            basename_match = parser.basename_match
        else:
            basename_match = '*'
            logger.warn("No basename_match attribute found in AtmosphereParser."+
                " We are matching on all filenames in '{0}'.".format(self.folder))
        
        # We need to determine all the actual grid points
        filename_match = os.path.join(models_folder, basename_match)
        filenames = glob(filename_match)
        
        if len(filenames) is 0:
            raise ValueError("No filenames found matching the wildmask '%s'" % (filename_match, ))
            
        for i, filename in enumerate(filenames):
            point = self.parser.parse_filename(filename)
            if i == 0:
                points = np.zeros((len(filenames), len(point)))
            points[i, :] = point
        

        # Before scaling the grid points we should be interpolating [Fe/H], [alpha/Fe], and surface gravity in logarithmic space
        for column in logarithmic_columns:
            points[:, column] = pow(10, points[:, column])
        
        self.logarithmic_columns = logarithmic_columns
        
        self.points = self.scale_gridpoints(points)
        self.filenames = filenames
        
        # What's the minimum number of model atmospheres required to interpolate in this many dimensions?
        if pow(2, self.points.shape[1]) > len(self.filenames):
            logger.warn('Warning: Only %i atmosphere files were found with "%s"' % (len(self.filenames), filename_match, ))
    
    
    def scale_gridpoints(self, points):
        """Scales the grid points to be between zero and unity.
        
        Parameters
        ----
        points : `np.array`
            A grid of points to scale down to between zero and unity.
        """
        
        try: (self.scales, self.offsets)
        except AttributeError: self.get_scales(points)
        
        points = (points - self.offsets)/self.scales
            
        return points
    
    
    def get_scales(self, points):
        """Determines scaling and offset values for the grid points such that
        all dimensions are scaled to between zero and unity.
        
        Parameters
        ----
        points : `np.array`
            The physical points in order to determine a unit scale from.
        """
        
        n = points.shape[1]
        self.offsets = np.zeros(n)
        self.scales = np.ones(n)
        
        for column in xrange(n):
            offset, scale = np.min(points[:, column]), np.max(points[:, column])
            
            self.offsets[column] = offset
            self.scales[column] = scale
        
        return np.vstack([self.offsets, self.scales])
    
    
    def scale_point(self, point):
        """Scales a single grid point to be between zero and unity.
        
        Parameters
        ----
        point : list containing floats
            The grid point to be scaled down to between zero and unity.
        """
        
        try: (self.scales, self.offsets)
        except AttributeError: self.get_scales()
        
        return (point - self.offsets)/self.scales
        
        
    def rescale_point(self, point):
        """Re-scales a single grid point from between zero and unity to its
        real, physical value.
        
        Parameters
        ----
        point : list containing floats
            The grid point (between zero and unity) to scale back to a physical
            value.
        """
        
        try: (self.scales, self.offsets)
        except AttributeError: raise ArithmeticError("No scales or offsets have been applied to the grid -- cannot descale point.")
        
        return point * self.scales + self.offsets
    

    def interpolate(self, filename, teff, logg, fe_h, alpha_fe, xi, solar_abundances=None, 
        molecules=None, clobber=False):
        """
        Interpolates a model photosphere for the stellar parameters provided and writes 
        the photosphere to a MOOG-compatible file.

        Parameters
        ----------
        filename : str
            The output filename for the interpolated model atmosphere.

        teff : float
            The effective temperature to interpolate to.

        logg : float
            The effective surface gravity to interpolate to.

        fe_h : float
            The overall metallicity (e.g., [Fe/H]) to interpolate to.

        alpha_fe : float
            The overall [alpha/Fe] to interpolate to.
        """

        logger.info("Generating model atmosphere with Teff = {0:.0f}, xi = {1:.2f} km/s, logg = {2:.2f}, [Fe/H] = {3:.2f}, [alpha/Fe] = {4:.2f}"
            .format(teff, xi, logg, fe_h, alpha_fe))
        thermal_structure = self.interpolate_thermal_structure(teff, logg, fe_h, alpha_fe)

        # Writing/Reading atmospheres is done by the parser
        return self.parser.write_atmosphere(filename, teff, logg, fe_h, xi, thermal_structure,
            solar_abundances, molecules, clobber)


    def interpolate_thermal_structure(self, *point):
        """
        Linearly interpolates through grid points of stellar atmospheres to
        a given point.
        
        Parameters
        ----
        point : list of floats-like objects
            A list containing the physical values to interpolate to. This is
            generally ordered to include temperature, surface gravity,
            metallicity and alpha enhancement. 
        """
        
        # Convert to array first
        neighbours = 1
        point = np.array(point)
        
        # Any logarithmic columns?
        if len(self.logarithmic_columns) > 0:
            point[self.logarithmic_columns] = pow(10, point[self.logarithmic_columns])
        
        point = self.scale_point(point)
        
        # Check for an exact match
        exact_match = np.where((self.points == point).all(axis=1))[0]        
        if len(exact_match) == 1:
            return self.parser.parse_thermal_structure(self.filenames[exact_match])
            
        elif len(exact_match) > 1:
            raise ValueError("Found multiple exact matches in atmospheric gridpoints.")
            
        # No exact match found - let's interpolate
        # Find the nearest grid points in all dimensions
        indices = []
        diff = self.points - point
        
        for column in xrange(diff.shape[1]):
            if np.all(self.points[:, column] == point[column]):
                continue
            
            values = np.unique(diff[:, column])
            
            negvals = np.sort(values[np.where(values < 0)])[-neighbours:]
            posvals = np.sort(values[np.where(values > 0)])[:neighbours]

            idx = []
            # find matching points first            
            match = np.where(self.points[:, column] == point[column])[0]
            idx.extend(list(match))
            
            if len(match) == 0:
                # Do neighbours
                
                for negval in negvals:
                    idx.extend(list(np.where(diff[:, column] == negval)[0]))
                
                for posval in posvals:
                    idx.extend(list(np.where(diff[:, column] == posval)[0]))
            
            if column == 0:
                indices.extend(idx)
            
            else:
                indices = set(indices).intersection(idx)
        
        # Generate our subgrid of points                    
        indices = list(indices)
        subset_points = np.copy(self.points[indices])
        
        # Load in our subgrid values
        for num, index in enumerate(indices):
            deck = self.parser.parse_thermal_structure(self.filenames[index])

            # Any logarithmic structure properties?
            if self.parser.logarithmic_structure_properties is not None:
                if not hasattr(self, "_warned_log_scaling"):
                    logger.warn("Logarithmically scaling thermal structure properties before interpolation: {0}"
                        .format(self.parser.logarithmic_structure_properties))
                    self._warned_log_scaling = True

                for logarithmic_structure_property in self.parser.logarithmic_structure_properties:
                    # Are we dealing with a structured array?
                    if deck.dtype.names is not None:
                        deck[logarithmic_structure_property] = 10**deck[logarithmic_structure_property]
                    else:
                        deck[:, logarithmic_structure_property] = 10**deck[:, logarithmic_structure_property]
            
            if isinstance(deck, np.core.records.recarray):
                if num == 0:
                    subset_values = np.zeros((len(subset_points), deck.shape[0], len(self.parser.structure_properties_to_interpolate)))

                for j, structure_property in enumerate(self.parser.structure_properties_to_interpolate):
                    subset_values[num, :, j] = deck[structure_property] 
                thermal_structure_shape = subset_values.shape[1:]

            else:
                if num == 0:
                    subset_values = np.zeros((len(subset_points), deck.shape[0], deck.shape[1]))
                
                if subset_values[num, :].shape[1] != deck.shape[1]:
                    logger.warn("Thermal structures in atmosphere grid points is different! %s != %s. Filled missing points with zeros." \
                        % (subset_values[num, :].shape, deck.shape, ))
                    
                    deck_new = np.zeros((subset_values[num, :].shape[0], np.max(subset_values[num, :].shape[1], deck.shape[1])))
                    deck_new[:deck.shape[0], :np.min([subset_values[num, :].shape[1], deck.shape[1]])] = deck
                    deck = deck_new

                subset_values[num, :] = deck
                thermal_structure_shape = deck.shape

        # Protect griddata from QHull errors
        superfluous_columns = []
        for column in xrange(subset_points.shape[1]):
            if len(np.unique(subset_points[:, column])) == 1:
                superfluous_columns.append(column)

        interpolated_point = scipy.delete(point, superfluous_columns, axis=0)        
        subset_points = scipy.delete(subset_points, superfluous_columns, axis=1)
        
        # Protect griddata from dimensionality complexes
        try:
            interpolated_deck = scipy.interpolate.griddata(
                subset_points, subset_values, 
                interpolated_point.reshape(1, len(interpolated_point))).reshape(thermal_structure_shape)
        except:
            raise ValueError("parameter outside of boundaries")
            
        if np.all(np.isnan(interpolated_deck)):
            raise ValueError("QHull interpolation fell over")
        
        if isinstance(deck, np.core.records.recarray):
            interpolated_deck = np.core.records.fromarrays(interpolated_deck.T,
                names=self.parser.structure_properties_to_interpolate,
                formats=["f8"] * len(self.parser.structure_properties_to_interpolate))

        # Any logarithmic structure properties that need rescaling?
        if self.parser.logarithmic_structure_properties is not None:
            for logarithmic_structure_property in self.parser.logarithmic_structure_properties:
                # Are we dealing with a structured array?
                if isinstance(interpolated_deck, np.core.records.recarray):
                    interpolated_deck[logarithmic_structure_property] = np.log10(interpolated_deck[logarithmic_structure_property])
                else:
                    interpolated_deck[:, logarithmic_structure_property] = np.log10(interpolated_deck[:, logarithmic_structure_property])

        return interpolated_deck
    

interpolator = AtmosphereInterpolator(
    os.path.abspath(os.path.join(os.path.dirname(os.path.expanduser(__file__)), 
        "atmospheres/castelli-kurucz-2004")), CastelliKuruczParser())