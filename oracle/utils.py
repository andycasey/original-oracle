# coding: utf-8

""" General utilities """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["atomic_number", "gaussian", "voigt", "CachedProperty"]

import numpy as np
from scipy.special import wofz


class CachedProperty(object):    
    """ Computes attribute value and caches it in the instance. """

    def __init__(self, method, name=None):
        self.method = method
        self.name = name or method.__name__
        self.__doc__ = method.__doc__
    
    def __get__(self, inst, cls): 
        if inst is None:
            return self
        elif self.name in inst.__dict__:
            return inst.__dict__[self.name]
        else:
            result = self.method(inst)
            inst.__dict__[self.name]=result
            return result    
    
    def __set__(self, inst, value):
        raise AttributeError("This property is read-only")
    
    def __delete__(self,inst):
        del inst.__dict__[self.name]


def gaussian(mu, sigma, x):
    return np.exp(-(x - mu)**2 / (2*sigma**2))


def voigt(mu, fwhm, amplitude, shape, x):
    """Evaluates a Voigt absorption profile across the `x`
    values with a given local continuum.

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    :param mu:
        The centroid of the line.

    :param fwhm:
        The full-width half-maximum of the Gaussian profile.

    :param amplitude:
        The amplitude of the line.

    :param shape:
        The shape parameter.

    :param x:
        An array of x-points.
    """
    n = len(x) if not isinstance(x, float) else 1

    profile = 1 / wofz(np.zeros((n)) + 1j * np.sqrt(np.log(2.0)) * shape).real
    profile *= amplitude * wofz(2*np.sqrt(np.log(2.0)) * (x - mu)/fwhm \
        + 1j * np.sqrt(np.log(2.0))*shape).real
    return profile
    

def atomic_number(element):
    """
    Return the atomic number of a given element.

    :param element:
        The short-hand notation for the element (e.g., Fe).

    :type element: str
    """
    
    if not isinstance(element, (unicode, str)):
        raise TypeError("element must be represented by a string-type")

    periodic_table = """H                                                  He
                        Li Be                               B  C  N  O  F  Ne
                        Na Mg                               Al Si P  S  Cl Ar
                        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                        Fr Ra Lr Rf Db Sg Bh Hs Mt Ds Rg Cn UUt"""
    
    lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
    actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
    
    periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
        .replace(" Ra ", " Ra " + actinoids + " ").split()
    del actinoids, lanthanoids
    
    if element not in periodic_table:
        return ValueError("element '{0}' is not known".format(element))

    return periodic_table.index(element) + 1