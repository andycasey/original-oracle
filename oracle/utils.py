# coding: utf-8

""" General utilities """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["atomic_number", "element", "reflect_about", "latexify"]

import logging

import numpy as np
from scipy.special import wofz

logger = logging.getLogger(__name__)

def stellar_jacobian(stellar_parameters, *args):
    """ Approximate the Jacobian of the stellar parameters and
    minimisation parameters, based on calculations from the Sun """

    logger.info("Updated approximation of the Jacobian")

    teff, vt, logg, feh = stellar_parameters[:4]

    # This is the black magic.
    full_jacobian = np.array([
        [ 5.4393e-08*teff - 4.8623e-04, -7.2560e-02*vt + 1.2853e-01,  1.6258e-02*logg - 8.2654e-02,  1.0897e-02*feh - 2.3837e-02],
        [ 4.2613e-08*teff - 4.2039e-04, -4.3985e-01*vt + 8.0592e-02, -5.7948e-02*logg - 1.2402e-01, -1.1533e-01*feh - 9.2341e-02],
        [-3.2710e-08*teff + 2.8178e-04,  3.8185e-03*vt - 1.6601e-02, -1.2006e-02*logg - 3.5816e-03, -2.8592e-05*feh + 1.4257e-03],
        [-1.7822e-08*teff + 1.8250e-04,  3.5564e-02*vt - 1.1024e-01, -1.2114e-02*logg + 4.1779e-02, -1.8847e-02*feh - 1.0949e-01]
    ])
    return full_jacobian.T


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


def latexify(labels, default_latex_labels=None):
    """
    Return a LaTeX-ified label.

    Args:
        labels (str or list-type of str objects): The label(s) to latexify. 
        default_latex_labels (dict): Dictionary of common labels to use.

    Returns:
        LaTeX-ified label.
    """

    common_labels = {
        "teff": "$T_{\\rm eff}$ (K)",
        "feh": "[Fe/H]",
        "logg": "$\log{g}$",
        "alpha": "[$\\alpha$/Fe]",
        "xi": "$\\xi$ (km s$^{-1}$)"
    }

    if default_latex_labels is not None:
        common_labels.update(default_latex_labels)
    
    listify = True
    if isinstance(labels, str):
        listify = False
        labels = [labels]

    latex_labels = []
    for label in labels:

        if label.startswith("doppler_sigma_"):
            color = ["blue", "green", "red", "ir"][int(label.split("_")[-1])]
            latex_labels.append("$\sigma_{\\rm doppler," + color + "}$ ($\\AA{}$)")

        else:
            latex_labels.append(common_labels.get(label, label))

    if not listify:
        return latex_labels[0]

    return latex_labels


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


def element(atomic_number):
    """
    Return the element of a given atomic number.

    :param atomic_number:
        The atomic number for the element in question (e.g., 26).

    :type atomic_number:
        int-like

    :returns:
        The short-hand element for a given atomic number.

    :rtype:
        str
    """

    atomic_number = int(atomic_number)
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
    return periodic_table[atomic_number - 1]