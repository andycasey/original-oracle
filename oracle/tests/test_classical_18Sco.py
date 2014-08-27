# coding: utf-8

""" Infer the properties of the Sun using the ClassicalModel model. """


import os
import unittest

from glob import glob
from textwrap import dedent

import specutils
import models

# Load spectra
# Create a StellarSpectrum Model

# Create a ClassicalModel model.
# Get initial guess for all parameters.
# Optimse stellar parameters
# Infer stellar parameters

class Infer_18Sco_ClassicalModel(unittest.TestCase):

    def setUp(self):
        """ Load the spectra and establish the StellarSpectrum model. """

        # Load the data.
        path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "data/benchmarks")
        self.data = map(specutils.Spectrum1D.load, 
            glob(os.path.join(path, "18Sco_?_?.txt")))

        assert len(self.data) > 0

        configuration = dedent("""
            model:
              redshift: yes
              outliers: yes
              continuum: no 
            
            settings:
              max_sampler_threads: 32

            balance:
              atomic_lines:
               - [4742.789, 22.00]
               - [4758.118, 22.00]
               - [4759.270, 22.00]
               - [4778.255, 22.00]
               - [4781.711, 22.00]
               - [4797.975, 22.00]
               - [4805.414, 22.00]
               - [4820.411, 22.00]
               - [4870.125, 22.00]
               - [4885.079, 22.00]
               - [5689.460, 22.00]
               - [5716.450, 22.00]
               - [5720.436, 22.00]
               - [5739.469, 22.00]
               - [5766.345, 22.00]
               - [5823.685, 22.00]
               - [5866.451, 22.00]
               - [6599.105, 22.00]
               - [6716.666, 22.00]
               - [7852.677, 22.00]
               - [4719.515, 22.10]
               - [4764.524, 22.10]
               - [4779.985, 22.10]
               - [4798.532, 22.10]
               - [4849.168, 22.10]
               - [4874.010, 22.10]
               - [6606.949, 22.10]
               - [4788.757, 26.00]
               - [4793.962, 26.00]
               - [4794.360, 26.00]
               - [4802.880, 26.00]
               - [4808.148, 26.00]
               - [4890.755, 26.00]
               - [4891.492, 26.00]
               - [5646.684, 26.00]
               - [5651.469, 26.00]
               - [5652.318, 26.00]
               - [5661.346, 26.00]
               - [5679.023, 26.00]
               - [5680.240, 26.00]
               - [5696.090, 26.00]
               - [5704.733, 26.00]
               - [5705.465, 26.00]
               - [5720.898, 26.00]
               - [5724.455, 26.00]
               - [5731.762, 26.00]
               - [5732.296, 26.00]
               - [5741.848, 26.00]
               - [5752.032, 26.00]
               - [5775.081, 26.00]
               - [5778.453, 26.00]
               - [5806.724, 26.00]
               - [5809.218, 26.00]
               - [5821.888, 26.00]
               - [5848.066, 26.00]
               - [5848.127, 26.00]
               - [5849.684, 26.00]
               - [5853.148, 26.00]
               - [5855.077, 26.00]
               - [5858.778, 26.00]
               - [5859.586, 26.00]
               - [5861.110, 26.00]
               - [5862.357, 26.00]
               - [6481.870, 26.00]
               - [6498.939, 26.00]
               - [6518.367, 26.00]
               - [6546.239, 26.00]
               - [6592.914, 26.00]
               - [6593.870, 26.00]
               - [6597.561, 26.00]
               - [6609.110, 26.00]
               - [6627.545, 26.00]
               - [6646.932, 26.00]
               - [6648.081, 26.00]
               - [6653.853, 26.00]
               - [6677.987, 26.00]
               - [6699.142, 26.00]
               - [6703.567, 26.00]
               - [6705.101, 26.00]
               - [6710.319, 26.00]
               - [6713.745, 26.00]
               - [6725.357, 26.00]
               - [6726.667, 26.00]
               - [6730.291, 26.00]
               - [6733.151, 26.00]
               - [7710.364, 26.00]
               - [7748.269, 26.00]
               - [7751.109, 26.00]
               - [7780.556, 26.00]
               - [7802.473, 26.00]
               - [7807.909, 26.00]
               - [7844.559, 26.00]
               - [7879.758, 26.00]
               - [4720.149, 26.10]
               - [4731.453, 26.10]
               - [4833.197, 26.10]
               - [6516.080, 26.10]
               - [7711.723, 26.10]
              wavelength_tolerance: 0.1
              wavelength_contribution: 2
              upper_smoothing_sigma: 0.5
            """)
        self.model = models.StellarSpectrum(configuration, self.data)


    def runTest(self):

        """
        foobar = models.AbsorptionProfile(4754.04, method="gaussian", outliers=True)
        xopt = foobar.optimise(self.data[0], np.ones(len(self.data[0].disp)))

        xopt_dict = dict(zip(foobar.parameters, xopt))
        foobar_spec = foobar(dispersion=self.data[0].disp, **xopt_dict)

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(self.data[0].disp, self.data[0].flux, 'k')
        ax.plot(foobar_spec[:,0], foobar_spec[:,1], 'b')

        raise a

        opt_theta = self.model.channels[0].optimise()
        result = self.model.channels[0].infer(opt_theta)


        raise a
        """

        optimised_parameters = self.model.optimise(self.data)

        raise a
        posterior, sampler, info = self.model.infer(self.data,
            optimised_parameters)


if __name__ == "__main__":
    foobar = Infer_18Sco_ClassicalModel()
    foobar.setUp()
    print("model foobar is set up")
    foobar.runTest()