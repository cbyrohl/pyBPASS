"""
Tests for the spectral_synthesis module. Requires local copies of the BPASS
data releases it is supposed to work with.
"""

from unittest import TestCase
from .. import spectral_synthesis
from . import config
from ddt import ddt, data, unpack
import numpy as np
import os
import astropy.units as u


@ddt
class TestBPASSsedDatabase(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path = config.BPASSpath
        cls.version = "2.2.1"

        bpath = os.path.join(cls.path, "bpass_v2.2.1_imf_chab300")

        path = os.path.join(bpath, "spectra-bin-imf_chab300.z040.dat.gz")
        if not os.path.isfile(path):
            path = path.rstrip(".gz")
        cls.bin_chab300_z040 = np.loadtxt(path)

        path = os.path.join(bpath, "spectra-bin-imf_chab300.z030.dat.gz")
        if not os.path.isfile(path):
            path = path.rstrip(".gz")
        cls.bin_chab300_z030 = np.loadtxt(path)

        cls.db_chab300_bin = spectral_synthesis.BPASSsedDatabase(
            cls.path,
            cls.version,
            "chab300",
            "bin"
        )

        cls.lam_min = (
            (200*u.eV).to(u.angstrom, equivalencies=u.spectral())
        ).value
        cls.lam_max = (
            (13.6*u.eV).to(u.angstrom, equivalencies=u.spectral())
        ).value
        cls.db_chab300_bin_lam_cut = spectral_synthesis.BPASSsedDatabase(
            cls.path,
            cls.version,
            "chab300",
            "bin",
            lam_min=cls.lam_min,
            lam_max=cls.lam_max
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
        # Delaunay triangulation is not unique for regular grid, can in
        # principle get two answers
        opt1 = 0.5*(
            self.__class__.bin_chab300_z040[:, 1] +
            self.__class__.bin_chab300_z030[:, 2]
        )
        opt2 = 0.5*(
            self.__class__.bin_chab300_z040[:, 2] +
            self.__class__.bin_chab300_z030[:, 1]
        )
        self.assertTrue(
            np.allclose(sed, opt1, atol=0.0)
            or
            np.allclose(sed, opt2, atol=0.0)
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

    @data(
        1.0645117e-09,
        8.359301e-10,
        8.000402e-10,
    )
    def test_nan_spectra(self, Z):
        """
        Some of the values for Z used here produced NaN spectra in early
        versions when not much attention was paid to data types.
        """
        db = self.__class__.db_chab300_bin

        age = 1e7
        M = 1e6
        Z = np.array(Z, dtype=np.float32)

        with self.assertWarns(UserWarning) as cm:
            lam, sed = db.interpolate(Z, age, masses=M)
        self.assertEqual(
            len(cm.warnings), 1
        )
        self.assertTrue(
            'Input metallicities for interpolation outside of available'
            in str(cm.warnings[0].message)
        )
        self.assertFalse(
            np.any(np.isnan(sed))
        )

        return

    @unpack
    @data(
        (2e-5, 1e6),
        (0.04, 1e11),
        (1e-3, 1e8)
    )
    def test_float32(self, Z, age):
        db = self.__class__.db_chab300_bin

        age = np.array(age, dtype=np.float32)
        z = np.array(Z, dtype=np.float32)

        lam, sed = db.interpolate(z, age)
        self.assertFalse(
            np.any(np.isnan(sed))
        )
        return

    def test_lam_cut(self):
        db = self.__class__.db_chab300_bin
        db_cut = self.__class__.db_chab300_bin_lam_cut

        # at given grid point
        z = 0.040
        age = 1e6
        lam, sed = db.interpolate(z, age)
        lam_cut, sed_cut = db_cut.interpolate(z, age)
        idx = np.where(
            (lam >= self.__class__.lam_min) &
            (lam <= self.__class__.lam_max)
        )
        self.assertTrue(
            np.all(lam[idx] == lam_cut)
        )
        self.assertTrue(
            np.allclose(sed[idx], sed_cut, atol=0.0)
        )
        return

    def test_em_decline(self):
        """
        Expect ionizing emissivity to do down with age.
        """
        db_cut = self.__class__.db_chab300_bin_lam_cut

        Z = 0.040

        age = 1e6
        _, sed_young = db_cut.interpolate(Z, age)

        age = 1e11
        _, sed_old = db_cut.interpolate(Z, age)

        self.assertTrue(
            np.all(sed_young > sed_old)
        )
        return


class TestBPASSemRatesDatabase(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = config.BPASSpath
        cls.version = "2.2.1"
        cls.db_chab300_bin = spectral_synthesis.BPASSemRatesDatabase(
            cls.path,
            cls.version,
            "chab300",
            "bin"
        )

        bpath = os.path.join(cls.path, "bpass_v2.2.1_imf_chab300")

        path = os.path.join(bpath, "ionizing-bin-imf_chab300.z040.dat.gz")
        if not os.path.isfile(path):
            path = path.rstrip(".gz")
        cls.bin_chab300_z040 = np.loadtxt(path)

        path = os.path.join(bpath, "ionizing-bin-imf_chab300.z030.dat.gz")
        if not os.path.isfile(path):
            path = path.rstrip(".gz")
        cls.bin_chab300_z030 = np.loadtxt(path)

        path = os.path.join(bpath, "ionizing-bin-imf_chab300.z020.dat.gz")
        if not os.path.isfile(path):
            path = path.rstrip(".gz")
        cls.bin_chab300_z020 = np.loadtxt(path)

    def test_interpolate_array(self):
        db = self.__class__.db_chab300_bin

        z = np.array([0.040, 0.035])
        age = np.array([1e6, 1e6])
        res = db.interpolate(z, age)
        self.assertEqual(
            (2, 4), res.shape
        )
        return

    def test_interpolate(self):
        db = self.__class__.db_chab300_bin

        # at given grid point
        z = 0.040
        age = 1e6
        res = db.interpolate(z, age)
        self.assertTrue(
            np.allclose(
                res,
                10**(self.__class__.bin_chab300_z040[0, 1:]), atol=0.0)
        )

        # between two points in Z
        z = 0.035
        age = 1e6
        res = db.interpolate(z, age)
        self.assertTrue(
            np.allclose(
                res,
                10**(
                    0.5*(self.__class__.bin_chab300_z030[0, 1:] +
                         self.__class__.bin_chab300_z040[0, 1:])
                ),
                atol=0.0
            )
        )

        # between two points in age
        z = 0.030
        age = 10**6.05
        res = db.interpolate(z, age)
        self.assertTrue(
            np.allclose(
                res,
                10**(
                    0.5*(self.__class__.bin_chab300_z030[0, 1:] +
                         self.__class__.bin_chab300_z030[1, 1:])
                ),
                atol=0.0
            )
        )

        # between four points
        z = 0.025
        age = 10**6.15
        res = db.interpolate(z, age)
        # Delaunay triangulation is not unique for regular grid, can in
        # principle get two answers
        opt1 = 10**(
            0.5*(
                self.__class__.bin_chab300_z020[1, 1:] +
                self.__class__.bin_chab300_z030[2, 1:]
            )
        )
        opt2 = 10**(
            0.5*(
                self.__class__.bin_chab300_z020[2, 1:] +
                self.__class__.bin_chab300_z030[1, 1:]
            )
        )

        self.assertTrue(
            np.allclose(res, opt1, atol=0.0)
            or
            np.allclose(res, opt2, atol=0.0)
        )

        # between four points
        z = 0.035
        age = 10**6.05
        res = db.interpolate(z, age)
        # Delaunay triangulation is not unique for regular grid, can in
        # principle get two answers
        opt1 =             10**(
            0.5*(
                self.__class__.bin_chab300_z030[1, 1:] +
                self.__class__.bin_chab300_z040[0, 1:]
            )
        )
        opt2 = 10**(
            0.5*(
                self.__class__.bin_chab300_z030[0, 1:] +
                self.__class__.bin_chab300_z040[1, 1:]
            )
        )

        self.assertTrue(
            np.allclose(res, opt1, atol=0.0)
            or
            np.allclose(res, opt2, atol=0.0)
        )

        return

    def test_Q_0_decline(self):
        """
        Expect Q_0 to go down as stars get older.
        """
        db = self.__class__.db_chab300_bin

        Z = 1.0
        with self.assertWarns(UserWarning):
            young = db.interpolate(Z, 0)
            old = db.interpolate(Z, 1e9)
        self.assertTrue(
            np.all(young > old)
        )

        Z = 0.0
        with self.assertWarns(UserWarning):
            young = db.interpolate(Z, 0)
            old = db.interpolate(Z, 1e9)
        self.assertTrue(
            np.all(young > old)
        )
        return


class TestBin_spectra(TestCase):

    def test_std(self):
        wave = np.linspace(1, 1000, num=2000)
        SEDs = np.random.random((10, len(wave)))
        bins = np.linspace(100, 500, num=5)
        edges = False
        wave_new, SEDs_new = spectral_synthesis.bin_spectra(
            wave, SEDs, bins, edges=edges
        )

        self.assertTrue(
            np.allclose(bins, wave_new)
        )
        self.assertEqual(
            (SEDs.shape[0], len(bins)), SEDs_new.shape
        )
        self.assertTrue(
            np.all(SEDs_new >= 0)
        )
        return

    def test_L_conservation(self):
        wave = np.linspace(1, 1000, num=20000, endpoint=True)
        SEDs = np.random.random((20, len(wave)))
        bins = np.linspace(1, 1000, num=10, endpoint=True)
        edges = True
        wave_new, SEDs_new = spectral_synthesis.bin_spectra(
            wave, SEDs, bins, edges=edges
        )
        self.assertEqual(
            (SEDs.shape[0], len(bins)-1), SEDs_new.shape
        )
        self.assertTrue(
            np.all(SEDs_new >= 0)
        )
        self.assertTrue(
            np.allclose(
                np.trapz(SEDs, x=wave, axis=1),
                np.sum(SEDs_new*(bins[1:]-bins[:-1]), axis=1)
            )
        )
        return
