import unittest

import numpy as np

import nanoepitools.math as nem

tol = 10e-6


class TestPlotArchive(unittest.TestCase):

    def test_llr_to_p(self):
        llrs = np.arange(-20, 20, 0.1)
        p = nem.llr_to_p(llrs)
        llr_rec = np.log((p / (1 - p)))
        max_deviation = np.abs(llr_rec - llrs).max()
        self.assertLess(max_deviation, tol)

    def test_llr_to_p_prior(self):
        llrs = np.arange(-20, 20, 0.1)
        p = nem.llr_to_p(llrs, prior=0.2)
        llr_rec = np.log((p / (1 - p))) - np.log(0.2/0.8)
        max_deviation = np.abs(llr_rec - llrs).max()
        self.assertLess(max_deviation, tol)

    def test_llr_to_p(self):
        p = np.arange(0.001, 0.999, 0.1)
        llr = nem.p_to_llr(p)
        p_rec = 1 / (1 + np.exp(-llr))
        max_deviation = np.abs(p - p_rec).max()
        self.assertLess(max_deviation, tol)

    def test_llr_to_p_with_prior(self):
        p = np.arange(0.001, 0.999, 0.1)
        llr = nem.p_to_llr(p, prior=0.2)
        p_rec = 1 / (1 + np.exp(-(llr + np.log(0.2/0.8))))
        max_deviation = np.abs(p - p_rec).max()
        self.assertLess(max_deviation, tol)

    def test_llr_to_p_and_back(self):
        llrs = np.arange(-20, 20, 0.1)
        p = nem.llr_to_p(llrs, prior=0.33)
        llr_rec = nem.p_to_llr(p, prior=0.33)
        max_deviation = np.abs(llr_rec - llrs).max()
        self.assertLess(max_deviation, tol)