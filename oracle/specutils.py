# coding: utf-8

""" Utilities for dealing with spectra. """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["Spectrum1D", "Spectrum"]

# Standard library
import logging
import os

# Third-party
import numpy as np
import pyfits

class Spectrum(object):
    """ A general spectrum class. """

    @classmethod
    def load(cls, filename, **kwargs):
        """
        Load a spectrum from a filename.

        :param filename:
            The filename to load.

        :returns:
            A spectrum.

        :rtype: :class:`Spectrum1D`
        """

        # Try as a Spectrum1D class first
        methods = (Spectrum1D.load, )#, OTHER_LOAD_FUNCTIONS

        for method in methods:
            try:
                spectrum = method(filename)
            except:
                continue

            else:
                if isinstance(spectrum, Spectrum1D) and spectrum.variance is None:
                    raise ValueError("no variance array found")

                return spectrum

        raise IOError("could not interpret spectrum in {0}".format(filename))


class Spectrum1D(object):
    """
    This is a temporary class holder for a Spectrum1D object until the
    :class:`astropy.specutils.Spectrum1D` module has sufficiently matured.
    """
    
    def __init__(self, disp, flux, variance=None, headers=None):
        """Initializes a `Spectrum1D` object with the given dispersion and flux
        arrays.
        
        :param disp:
            Dispersion of the spectrum (i.e. the wavelength points).
            
        :type disp:
            :class:`numpy.array`

        :param flux:
            Flux points for each `disp` point.

        :type flux:
            :class:`numpy.array`

        :param variance: 
            The variance in flux points for each dispersion point.

        :type variance:
            :class:`numpy.array`

        :returns:
            A spectrum.
        """

        if len(disp) != len(flux):
            raise ValueError("dispersion and flux must have the same length")

        if len(disp) == 0:
            raise ValueError("dispersion and flux cannot be empty arrays")
        
        self.disp = disp
        self.flux = flux
        if variance is None:
            variance = flux
        self.variance = variance
        # Better to send an extra array (ivariance) around than calculate it at
        # every single likelihood call.
        self.ivariance = 1.0/variance 
        if headers is not None:
            self.headers = headers
        else:
            self.headers = {}
        return None


    def copy(self):
        """ Creates a copy of the object. """

        return self.__class__(self.disp.copy(), self.flux.copy(),
            variance=self.variance, headers=self.headers)
    

    @classmethod
    def load(cls, filename, **kwargs):
        """
        Load a Spectrum1D from a given filename.
        
        :param filename:
            Path of the filename to load. Can be either simple FITS extension
            or an ASCII filename.

        :type filename:
            str
            
        :notes:
            If you are loading from an non-standard ASCII file, you can pass
            kwargs to :func:`numpy.loadtxt` through this function.
        """
        
        if not os.path.exists(filename):
            raise IOError("path {0} does not exist".format(filename))
        
        variance = None

        if filename.endswith(".fits"):
            image = pyfits.open(filename, **kwargs)
            
            header = image[0].header
            
            # Check for a tabular data structure
            if len(image) > 1 and image[0].data is None:

                names = [name.lower() for name in image[1].data.names]
                dispersion_key = "wave" if "wave" in names else "disp"
                
                disp, flux = image[1].data[dispersion_key], image[1].data["flux"]

                if "error" in names or "variance" in names:
                    variance_key = "error" if "error" in names else "variance"
                    variance = image[1].data[variance_key]

            else:

                # According to http://iraf.net/irafdocs/specwcs.php ....
                #li = a.headers["LTM1_1"] * np.arange(a.headers["NAXIS1"]) + a.headers["LTV1"]
                #a.headers["CRVAL1"] + a.headers["CD1_1"] * (li - a.headers["CRPIX1"])

                if np.all([key in header.keys() for key in ("CDELT1", "NAXIS1", "CRVAL1")]):
                    disp = header["CRVAL1"] + np.arange(header["NAXIS1"]) * header["CDELT1"]
            
                if "LTV1" in header.keys():
                    disp -= header["LTV1"] * header["CDELT1"]

                #disp -= header["LTV1"] if header.has_key("LTV1") else 0
                flux = image[0].data
            
                # Check for logarithmic dispersion
                if "CTYPE1" in header.keys() and header["CTYPE1"] == "AWAV-LOG":
                    disp = np.exp(disp)

            # Add the headers in
            headers = {}
            for row in header.items():
                key, value = row
                
                # Check the value is valid
                try:
                    str(value)

                except TypeError:
                    continue

                if len(key) == 0 or len(str(value)) == 0: continue
                
                if key in headers.keys():
                    if not isinstance(headers[key], list):
                        headers[key] = [headers[key]]
                    
                    headers[key].append(value)

                else:
                    headers[key] = value

            for key, value in headers.iteritems():
                if isinstance(value, list):
                    headers[key] = "\n".join(map(str, value))

        else:
            headers = {}
            # Try for variance too first
            try:
                disp, flux, variance = np.loadtxt(filename, unpack=True, **kwargs)
            except:
                disp, flux = np.loadtxt(filename, unpack=True, **kwargs)
            
        return cls(disp, flux, variance=variance, headers=headers)


    def save(self, filename, clobber=True):
        """
        Saves the spectrum to disk.
        
        :param filename:
            The filename to save the spectrum to.

        :type filename:
            str

        :param clobber: [optional]
            Whether to overwite the ``filename`` if it already exists.

        :type clobber:
            bool
        
        :raises IOError:
            If the filename exists and we are not asked to clobber it.
        """
        
        if os.path.exists(filename) and not clobber:
            raise IOError("Filename already exists and we have been asked not \
                to clobber it." % (filename, ))
        
        if not filename.endswith("fits"):
            # ASCII
            
            if self.variance is not None:
                data = np.hstack([
                    self.disp.reshape(-1, 1),
                    self.flux.reshape(-1, 1),
                    self.variance.reshape(-1, 1)
                    ])
            else:
                data = np.hstack([self.disp.reshape(len(self.disp), 1), 
                    self.flux.reshape(len(self.disp), 1)])
            
            np.savetxt(filename, data)
            
        else:
            # FITS
            crpix1, crval1 = 1, self.disp.min()
            
            cdelt1 = np.mean(np.diff(self.disp))
            
            test_disp = (crval1 + np.arange(len(self.disp), dtype=self.disp.dtype)\
                * cdelt1).astype(self.disp.dtype)
            
            if np.max(self.disp - test_disp) > 10e-2 or self.variance is not None:

                # Non-linear dispersion map, or we have variance information too
                # Create a tabular FITS format.

                col_disp = pyfits.Column(name="disp", format="1D", array=self.disp)
                col_flux = pyfits.Column(name="flux", format="1D", array=self.flux)

                if self.variance is not None:
                    col_variance = pyfits.Column(name="variance", format="1D",
                        array=self.variance)

                    table_hdu = pyfits.new_table([col_disp, col_flux, col_variance])

                else:
                    table_hdu = pyfits.new_table([col_disp, col_flux])

                # Create Primary HDU
                hdu = pyfits.PrimaryHDU()

                # Update primary HDU with headers
                for key, value in self.headers.iteritems():
                    if len(key) > 8:
                        # To deal with ESO compatibility
                        hdu.header.update("HIERARCH {0}".format(key), value)
                    
                    try:
                        hdu.header.update(key, value)

                    except ValueError:
                        logger.warn("Could not save header: {0} = {1}".format(
                            key, value))
                    
                # Create HDU list with our tables
                hdulist = pyfits.HDUList([hdu, table_hdu])

                hdulist.writeto(filename, clobber=clobber)

            else:
                # Linear dispersion map.
                # Create a PrimaryHDU file.

                # Ensure we have an array!
                hdu = pyfits.PrimaryHDU(np.array(self.flux))
                
                headers = self.headers.copy()
                headers.update({
                    "CRVAL1": crval1,
                    "CRPIX1": crpix1,
                    "CDELT1": cdelt1
                })
                
                for key, value in headers.iteritems():
                    if len(key) > 8:
                        # To deal with ESO compatibility
                        hdu.header.update("HIERARCH {0}".format(key), value)
                    
                    else:
                        try:
                            hdu.header.update(key, value)

                        except ValueError:
                            logger.warn("Could not save header: %s = %s".format(
                                key, value))
                
                hdu.writeto(filename, clobber=clobber)


def _cross_correlate(dispersion, observed_flux, template_flux):
    """
    Return a redshift by cross correlation of a template and observed spectra.

    :param dispersion:
        The dispersion points of the observed and template fluxes.

    :type dispersion:
        :class:`numpy.array`

    :param observed_flux:
        The observed fluxes for each dispersion point.

    :type observed_flux:
        :class:`numpy.array`

    :param template_flux:
        The template fluxes for each dispersion point.

    :type template_flux:
        :class:`numpy.array`

    :returns:
        A redshift.

    :rtype:
        float
    """

    # Be forgiving, although we shouldn't have to be.
    N = np.min(map(len, [dispersion, observed_flux, template_flux]))

    # Ensure an even number of points
    if N % 2 > 0:
        N -= 1

    dispersion = dispersion[:N]
    observed_flux = observed_flux[:N]
    template_flux = template_flux[:N]

    assert len(dispersion) == len(observed_flux)
    assert len(observed_flux) == len(template_flux)
    
    # Set up z array
    m = len(dispersion) / 2
    z_array = dispersion/dispersion[N/2] - 1.0
    
    # Apodize edges
    edge_buffer = 0.1 * (dispersion[-1] - dispersion[0])
    low_w_indices = np.nonzero(dispersion < dispersion[0] + edge_buffer)[0]
    high_w_indices = np.nonzero(dispersion > dispersion[-1] - edge_buffer)[0]

    apod_curve = np.ones(N, dtype='d')
    apod_curve[low_w_indices] = (1.0 + np.cos(np.pi*(1.0 - \
        (dispersion[low_w_indices] - dispersion[0])/edge_buffer)))/2.
    apod_curve[high_w_indices] = (1.0 + np.cos(np.pi*(1.0 - \
        (dispersion[-1] - dispersion[high_w_indices])/edge_buffer)))/2.

    apod_observed_flux = observed_flux * apod_curve
    apod_template_flux = template_flux * apod_curve

    fft_observed_flux = np.fft.fft(apod_observed_flux)
    fft_template_flux = np.fft.fft(apod_template_flux)
    template_flux_corr = (fft_observed_flux * fft_template_flux.conjugate())   \
        / np.sqrt(np.inner(apod_observed_flux, apod_observed_flux)             \
        * np.inner(apod_template_flux, apod_template_flux))

    correlation = np.fft.ifft(template_flux_corr).real

    # Reflect about zero
    ccf = np.zeros(N)
    ccf[:N/2] = correlation[N/2:]
    ccf[N/2:] = correlation[:N/2]
    
    # Get height and redshift of best peak
    h = ccf.max()
    
    # Scale the CCF
    ccf -= ccf.min()
    ccf *= (h/ccf.max())

    z_best = z_array[ccf.argmax()]    
    z_err = (np.ptp(z_array[np.where(ccf >= 0.5*h)])/2.35482)**2

    return z_best