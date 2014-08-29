# coding: utf-8

""" General utilities """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["atomic_number", "reflect_about"]

import numpy as np
from scipy.special import wofz


def reflect_about(a, limits):
    """
    Similar to :func:`numpy.clip`, except it just reflects about some limiting axes.

    :param a:
        The array of values to reflect.

    :type a:
        :class:`numpy.array`

    :param limits:
        The upper and lower limits to reflect about. Use ``None`` for no limit.

    :type limits:
        A two length tuple or list-type.

    :returns:
        The reflected array.

    :rtype:
        :class:`numpy.array`
    """

    lower, upper = limits
    if lower is not None:
        a[a < lower] = lower + (lower - a[a < lower])
    if upper is not None:
        a[a > upper] = upper - (a[a > upper] - upper)
    return a


def atomic_number(element):
    """
    Return the atomic number of a given element.

    :param element:
        The short-hand notation for the element (e.g., Fe).

    :type element:
        str

    :returns:
        The atomic number for a given element.

    :rtype:
        int
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