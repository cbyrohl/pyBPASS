"""
Tests for _sed_tools Fortran extension module.
"""
from unittest import TestCase
from .. import _sed_tools
import numpy as np


class TestBin_sed(TestCase):

    def test_std(self):
        lam = np.linspace(1, 10)
        L_lam = np.ones_like(lam)
        bin_edges = np.array([4, 5])
        L_lam_new = np.empty((1))
        _sed_tools.sed_tools.bin_sed(lam, L_lam, bin_edges, L_lam_new)
        self.assertAlmostEqual(
            1.0, L_lam_new[0]
        )
        return

    def test_edge_cases(self):
        lam = np.linspace(1, 10)
        L_lam = np.ones_like(lam)

        # bin edges at ends of lam
        bin_edges = np.array([1, 5, 10])
        L_lam_new = np.empty((2))
        _sed_tools.sed_tools.bin_sed(lam, L_lam, bin_edges, L_lam_new)
        self.assertAlmostEqual(
            1.0, L_lam_new[0]
        )
        self.assertAlmostEqual(
            1.0, L_lam_new[1]
        )

        # zero width bin
        bin_edges = np.array([1, 5, 5, 10])
        L_lam_new = np.empty((3))
        _sed_tools.sed_tools.bin_sed(lam, L_lam, bin_edges, L_lam_new)
        for j, x in enumerate([1, 0, 1]):
            self.assertAlmostEqual(
                1.0, L_lam_new[0]
            )
        return

    def test_L_conservation(self):
        lam = np.linspace(1, 5, num=50)
        L_lam = np.random.random(len(lam))
        L_lam_new = np.empty((1))
        bin_edges = np.array([lam[0], lam[-1]])
        _sed_tools.sed_tools.bin_sed(lam, L_lam, bin_edges, L_lam_new)
        self.assertAlmostEqual(
            np.trapz(L_lam, x=lam), L_lam_new[0]*np.diff(bin_edges)[0]
        )
        return

    def test_L_conservation_2(self):
        lam = np.linspace(1, 5, num=50)
        L_lam = np.random.random(len(lam))
        bin_edges = np.array([lam[0], 3, lam[-1]])
        L_lam_new = np.empty((2))
        _sed_tools.sed_tools.bin_sed(lam, L_lam, bin_edges, L_lam_new)
        self.assertAlmostEqual(
            np.trapz(L_lam, x=lam),
            np.sum(L_lam_new*np.diff(bin_edges)),
            places=4
        )
        return
