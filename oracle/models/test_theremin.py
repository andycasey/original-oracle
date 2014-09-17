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
    "mask": [
      [4720.35, 4725]
    ],
    "ThereminModel": {
      # line_list_filename can be binary or ascii -- we will figure out which
      "line_list_filename": "~/codes/oracle/oracle/tests/line_list.data",
      "atomic_lines": [
        # This contains our atomic lines of interest. These must be located in
        # the line_list_filename as well. 
        [4720.149, 26.1],
        [4731.453, 26.1],
        [4758.118, 22.0],
        [4759.270, 22.0],
        [4764.524, 22.1],
        [4778.255, 22.0],
        [4781.711, 22.0],
        [4788.757, 26.0],
        [4793.962, 26.0],
        [4794.360, 26.0],
        [4797.975, 22.0],
        [4798.532, 22.1],
        [4802.880, 26.0],
        [4808.148, 26.0],
        [4820.411, 22.0],
        [4833.197, 26.1],
        [4874.010, 22.1],
        [4890.759, 26.0],
        [4891.492, 26.0],
        [5651.469, 26.0],
        [5652.318, 26.0],
        [5661.346, 26.0],
        [5679.023, 26.0],
        [5680.240, 26.0],
        [5689.460, 22.0],
        [5696.090, 26.0],
        [5705.465, 26.0],
        [5716.450, 22.0],
        [5720.436, 22.0],
        [5731.762, 26.0],
        [5732.296, 26.0],
        [5741.848, 26.0],
        [5775.081, 26.0],
        [5778.453, 26.0],
        [5849.684, 26.0],
        [5853.148, 26.0],
        [5855.077, 26.0],
        [5858.778, 26.0],
        [5866.451, 22.0],
        [6481.870, 26.0],
        [6494.980, 26.0],
        [6498.939, 26.0],
        [6516.080, 26.1],
        [6518.367, 26.0],
        [6546.239, 26.0],
        [6592.914, 26.0],
        [6593.870, 26.0],
        [6597.561, 26.0],
        [6599.105, 22.0],
        [6609.110, 26.0],
        [6627.545, 26.0],
        [6648.080, 26.0],
        [6677.987, 26.0],
        [6699.142, 26.0],
        [6703.567, 26.0],
        [6713.745, 26.0],
        [6716.666, 22.0],
        [6733.151, 26.0],
        [6739.521, 26.0],
        [7852.677, 22.0],
        [7710.364, 26.0],
        [7748.269, 26.0],
        [7711.723, 26.1],
        [7723.210, 26.0]
      ]
    }
}
filenames = glob("../tests/data/benchmarks/18Sco/18Sco_narval*noresample*")
data = map(theremin.specutils.Spectrum1D.load, filenames)

model = theremin.SpectrumModel(config, data)
