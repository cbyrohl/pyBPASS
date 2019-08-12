import numpy as _np
from pyBPASS.database import BPASSdatabase as _BPASSdatabase

_constants = {
    '2.2.1': {
        'L_sun': 3.848e26
    }
}


class BPASSsedDatabase(_BPASSdatabase):
    """
    Representation of and interface for single BPASS SED grid defined by an IMF
    and a population type.

    Parameters
    ----------
    dbdtype : dtype, optional
        Data type to use for arrays read from disk. Defaults to numpy's
        float64.
    lam_min : float, optional
        Minimum wavelength [angstrom] of the interpolated spectra. Defaults to
        None, in which case the minimum available wavelength in the database
        will be taken.
    lam_max : float, optional
        Maximum wavelength [angstrom] of the interpolated spectra. Defaults to
        None, in which case the maximum available wavelength in the database
        will be taken.

    Attributes
    ----------
    log_ages : list
        Log of the population ages [yr] at which BPASS provides SEDs for this
        population.
    wavelengths : array
        The wavelengths [angstrom] at which the SEDs are sampled.
    SEDgrid : array
        Array of fluxes [L_sun/angstrom] provided by BPASS as a function of
        metallicity, age and wavelength.
    Lsun : float
        The solar luminosity value [J/s] assumed in the BPASS version an
        instance of this class corresponds to.
    """
    def __init__(self, path, version, imf, popType, dbdtype=_np.float64,
                 lam_min=None, lam_max=None):
        super().__init__(path, version, imf, popType, "spectra", dbdtype)
        self.Lsun = _constants[self.version]['L_sun']
        self._constructGrid(dbdtype, lam_min, lam_max)
        self._constructInterpolator()
        return

    def _constructGrid(self, dbdtype, lam_min, lam_max):
        """
        Load all available SEDs into memory.

        Builds a 3D array of flux [L_sun/angstrom] as function of
        (Z,age,lambda) for a simple stellar population of 1e6 M_sun.
        """

        # load first file to extract some info
        fArr = _np.loadtxt(self.flist[0], dtype=dbdtype)

        wavelengths = fArr[:, 0]
        if lam_min is not None:
            idx_min = _np.searchsorted(wavelengths, lam_min, side='left')
        else:
            idx_min = None
        if lam_max is not None:
            idx_max = _np.searchsorted(wavelengths, lam_max, side='right')
        else:
            idx_max = None
        self.wavelengths = fArr[:, 0][idx_min:idx_max]
        self.log_ages = _np.array([
            (6+0.1*(n-2)) for n in range(2, fArr.shape[1]+1)
        ], dtype=dbdtype)
        self._aMin = self.log_ages[0]
        self._aMax = self.log_ages[-1]

        # can now allocate the actual grid
        self.grid = _np.empty(
            (len(self.metallicities),
             len(self.log_ages),
             len(self.wavelengths)),
            dtype=dbdtype
        )

        # already have the first file loaded
        self.grid[0, :, :] = fArr[idx_min:idx_max, 1:].T
        # load the rest
        for j, f in enumerate(self.flist[1:]):
            self.grid[j+1, :, :] = \
                _np.loadtxt(f, dtype=dbdtype)[idx_min:idx_max, 1:].T
        return

    def interpolate(self, metallicities, ages, masses=1):
        """
        Interpolate spectra for stellar populations.

        The units of the returned SEDs are Solar luminosities per
        Angstrom. Computes SEDs for stellar populations formed in single
        instantaneous bursts.

        Interpolation is done in metallicity-log(age) space by triangulation
        and subsequent barycentric interpolation.

        Parameters
        ----------
        metallicities : array, shape `(N)` or float
            The metallicities of the stellar populations given as the mass
            fractions in metals.
        ages : array, shape `(N)` or float
            The ages of the stellar populations [yr].
        masses : optional, array, shape `(N)` or float
            The masses of the stellar populations [1e6 M_sun].

        Returns
        -------
         : array, shape `(N_lam)`
            Array of wavelengths [angstrom] provided by the database.
         : array, shape `(N, N_lam)`
            Interpolated SED [L_sun/angstrom] normalized to 1e6 M_sun for each
            input stellar population at all wavelengths provided by the
            database.
        """
        return self.wavelengths, self._interpolate(metallicities, ages, masses)


class BPASSionRatesDatabase(_BPASSdatabase):
    def __init__(self, path, version, imf, popType, dbdtype=_np.float64,
                 lam_min=None, lam_max=None):
        super().__init__(path, version, imf, popType, "ionizing", dbdtype)
        self._constructGrid(dbdtype)
        self._constructInterpolator()
        return

    def _constructGrid(self, dbdtype):
        # load first file to extract some info
        fArr = _np.loadtxt(self.flist[0], dtype=dbdtype)

        self.log_ages = fArr[:, 0][:]
        self._aMin = self.log_ages[0]
        self._aMax = self.log_ages[-1]

        # can now allocate the actual grid
        self.grid = _np.empty(
            (len(self.metallicities),
             len(self.log_ages),
             fArr.shape[1]-1),
            dtype=dbdtype
        )

        # already have the first file loaded
        self.grid[0, :, :] = fArr[:, 1:]
        # load the rest
        for j, f in enumerate(self.flist[1:]):
            self.grid[j+1, :, :] = \
                _np.loadtxt(f, dtype=dbdtype)[:, 1:]
        return

    def interpolate(self, metallicities, ages, masses=1):
        """
        Interpolate ionizing photon rates for stellar populations.

        Computes:
            log(Nion):
                ionizing photon production rate
            log(L_alpha in ergs/s):
                Balmer H line luminosity, assuming =log(Nion/s)-11.87
            log (L_FUV in ergs/s/A):
                luminosity in the FUV band (mean flux from 1556 to 1576A)
            log (L_NUV in ergs/s/A):
                luminosity in the NUV band (mean flux from 2257 to 2277A)

        Interpolation is done in metallicity-log(age) space by triangulation
        and subsequent barycentric interpolation.

        Parameters
        ----------
        metallicities : array, shape `(N)` or float
            The metallicities of the stellar populations given as the mass
            fractions in metals.
        ages : array, shape `(N)` or float
            The ages of the stellar populations [yr].
        masses : optional, array, shape `(N)` or float
            The masses of the stellar populations [1e6 M_sun].

        Returns
        -------
         : array, shape `(N, 4)`
            Interpolated values as specified above and normalized to 1e6 M_sun
            for each input stellar population.
        """
        return self._interpolate(metallicities, ages, masses)
