# coding: utf-8

from __future__ import absolute_import, print_function

""" Test ThereminModel Solver """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

from glob import glob

import theremin


config = {
    "model": {
      "continuum": False,
      "redshift": False,
      "doppler_broadening": True
    },
    "ThereminModel": {
      "clean_line_list_filename": "../si/linedata/HERMES_clean.dat",
      "blend_line_list_filename": "../si/linedata/HERMES_blending.dat",
    }
}
filenames = glob("../tests/data/benchmarks/18Sco/18Sco_narval*noresample*")
data = map(theremin.specutils.Spectrum1D.load, filenames)
data.sort(key=lambda s: s.disp[0])

model = theremin.SpectrumModel(config, data)
