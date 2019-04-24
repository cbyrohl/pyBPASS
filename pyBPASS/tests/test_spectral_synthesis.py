"""
Tests for the spectral_synthesis module. Requires local copies of the BPASS
data releases it is supposed to work with.
"""

from unittest import TestCase
from .. import spectral_synthesis
from . import config
import numpy as np
import os


class TestBPASSsedDatabase(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path = config.BPASSpath
        cls.version = "2.2.1"
        cls.bin_chab300_z040 = np.loadtxt(
            os.path.join(
                cls.path,
                "bpass_v2.2.1_imf_chab300/spectra-bin-imf_chab300.z040.dat.gz"
            )
        )
        cls.bin_chab300_z030 = np.loadtxt(
            os.path.join(
                cls.path,
                "bpass_v2.2.1_imf_chab300/spectra-bin-imf_chab300.z030.dat.gz"
            )
        )
        cls.db_chab300_bin = spectral_synthesis.BPASSsedDatabase(
            cls.path,
            cls.version,
            "chab300",
            "bin"
        )
        return

    @classmethod
    def tearDownClass(cls):
        return

    def test_interpolate_clip(self):
        db = self.__class__.db_chab300_bin

        # clip Z
        z = 1
        age = 1e6
        with self.assertWarns(UserWarning):
            lam, sed = db.interpolate(z, age)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.any(sed > 0.0)
            )
        self.assertTrue(
            np.allclose(sed, self.__class__.bin_chab300_z040[:, 1], atol=0.0)
        )

        # clip age
        z = 0.040
        age = 1e5
        with self.assertWarns(UserWarning):
            lam, sed = db.interpolate(z, age)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.allclose(sed, self.__class__.bin_chab300_z040[:, 1], atol=0.0)
        )
        return

    def test_interpolate(self):
        db = self.__class__.db_chab300_bin

        # at given grid point
        z = 0.040
        age = 1e6
        lam, sed = db.interpolate(z, age)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.allclose(sed, self.__class__.bin_chab300_z040[:, 1], atol=0.0)
        )

        # between two points in Z
        z = 0.035
        age = 1e6
        lam, sed = db.interpolate(z, age)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.allclose(
                sed,
                0.5*(self.__class__.bin_chab300_z030[:, 1] +
                     self.__class__.bin_chab300_z040[:, 1]),
                atol=0.0
            )
        )

        # between two points in age
        z = 0.030
        age = 10**6.05
        lam, sed = db.interpolate(z, age)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z030[:, 0])
        )
        self.assertTrue(
            np.allclose(
                sed,
                0.5*(self.__class__.bin_chab300_z030[:, 1] +
                     self.__class__.bin_chab300_z030[:, 2]),
                atol=0.0
            )
        )

        # between four points
        z = 0.035
        age = 10**6.05
        lam, sed = db.interpolate(z, age)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.allclose(
                sed,
                0.5*(
                    self.__class__.bin_chab300_z040[:, 1] +
                    self.__class__.bin_chab300_z030[:, 2]
                ),
                atol=0.0
            )
        )
        return

    def test_interpolate_array(self):
        db = self.__class__.db_chab300_bin

        z = np.array([0.040, 0.035])
        age = np.array([1e6, 1e6])
        lam, sed = db.interpolate(z, age)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.allclose(sed[0, :],
                        self.__class__.bin_chab300_z040[:, 1],
                        atol=0.0)
        )
        self.assertTrue(
            np.allclose(
                sed[1, :],
                0.5*(self.__class__.bin_chab300_z030[:, 1] +
                     self.__class__.bin_chab300_z040[:, 1]),
                atol=0.0
            )
        )
        return

    def test_interpolate_with_mass(self):
        db = self.__class__.db_chab300_bin

        z = 0.040
        age = 1e6
        mass = 1e2/1e6
        with self.assertWarns(UserWarning):
            lam, sed = db.interpolate(z, age, masses=mass)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.allclose(sed,
                        mass*self.__class__.bin_chab300_z040[:, 1], atol=0.0)
        )
        return

    def test_interpolate_with_mass_array(self):
        db = self.__class__.db_chab300_bin

        z = np.array([0.040, 0.035])
        age = np.array([1e6, 1e6])
        mass = np.array([1, 1e3])
        lam, sed = db.interpolate(z, age, masses=mass)
        self.assertTrue(
            np.all(lam == self.__class__.bin_chab300_z040[:, 0])
        )
        self.assertTrue(
            np.allclose(sed[0, :],
                        mass[0]*self.__class__.bin_chab300_z040[:, 1],
                        atol=0.0)
        )
        self.assertTrue(
            np.allclose(
                sed[1, :],
                mass[1]*0.5*(self.__class__.bin_chab300_z030[:, 1] +
                             self.__class__.bin_chab300_z040[:, 1]),
                atol=0.0
            )
        )
        return
