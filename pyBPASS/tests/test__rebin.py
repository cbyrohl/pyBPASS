"""
Tests for _rebin Fortran extension module.
"""
from unittest import TestCase
from .. import _rebin
import numpy as np


class Testrebin_sed(TestCase):

    def test_std(self):
        lam = np.linspace(1, 10)
        L_lam = np.ones_like(lam)
        bin_edges = np.array([4, 5])
        lam_new = np.empty((1))
        L_lam_new = np.empty((1))
        _rebin.rebin.rebin_sed(lam, L_lam, bin_edges, lam_new, L_lam_new)
        self.assertAlmostEqual(
            4.5, lam_new[0]
        )
        self.assertAlmostEqual(
            1.0, L_lam_new[0]
        )
        return

    def test_edge_cases(self):
        lam = np.linspace(1, 10)
        L_lam = np.ones_like(lam)

        # bin edges at ends of lam
        bin_edges = np.array([1, 5, 10])
        lam_new = np.empty((2))
        L_lam_new = np.empty((2))
        _rebin.rebin.rebin_sed(lam, L_lam, bin_edges, lam_new, L_lam_new)
        self.assertAlmostEqual(
            3, lam_new[0]
        )
        self.assertAlmostEqual(
            7.5, lam_new[1]
        )
        self.assertAlmostEqual(
            1.0, L_lam_new[0]
        )
        self.assertAlmostEqual(
            1.0, L_lam_new[1]
        )

        # zero width bin
        bin_edges = np.array([1, 5, 5, 10])
        lam_new = np.empty((3))
        L_lam_new = np.empty((3))
        _rebin.rebin.rebin_sed(lam, L_lam, bin_edges, lam_new, L_lam_new)
        for j, x in enumerate([3, 5, 7.5]):
            self.assertAlmostEqual(
                x, lam_new[j]
            )
        for j, x in enumerate([1, 0, 1]):
            self.assertAlmostEqual(
                1.0, L_lam_new[0]
            )
        return

    def test_L_conservation(self):
        lam = np.linspace(1, 5, num=50)
        L_lam = np.random.random(len(lam))
        lam_new = np.empty((1))
        L_lam_new = np.empty_like(lam_new)
        bin_edges = np.array([lam[0], lam[-1]])
        _rebin.rebin.rebin_sed(lam, L_lam, bin_edges, lam_new, L_lam_new)
        self.assertAlmostEqual(
            np.trapz(L_lam, x=lam), L_lam_new[0]*np.diff(bin_edges)[0]
        )
        return
