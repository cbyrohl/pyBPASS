import warnings as _warnings
import numpy as _np
from pyBPASS.database import BPASSdatabase as _BPASSdatabase
from pyBPASS import _sed_tools

_constants = {
    '2.2.1': {
        'L_sun': 3.848e26
    }
}


class BPASSsedDatabase(_BPASSdatabase):
    """
    This class models a database for the SEDs as a function of metallicity and
    age of stellar populations for fixed IMF and population type.

    Attributes
    ----------
    wavelengths : array
        The wavelengths [angstrom] at which the SEDs are sampled.
    SEDgrid : array
        Array of fluxes [L_sun/angstrom] provided by BPASS as a function of
        metallicity, age and wavelength.
    Lsun : float
        The solar luminosity value [J/s] assumed in the BPASS version an
        instance of this class corresponds to.
    """
    def __init__(self, basepath, version, imf, popType, dbdtype=_np.float64,
                 lam_min=None, lam_max=None):
        """
        Parameters
        ----------
        lam_min : float, optional
            Minimum wavelength [angstrom] of the interpolated spectra. Defaults
            to None, in which case the minimum available wavelength in the
            database will be taken.
        lam_max : float, optional
            Maximum wavelength [angstrom] of the interpolated spectra. Defaults
            to None, in which case the maximum available wavelength in the
            database will be taken.
        """
        if lam_min is not None and lam_max is not None:
            if lam_min > lam_max:
                raise ValueError("lam_min is larger than lam_max!")
        super().__init__(basepath, version, imf, popType, "spectra", dbdtype)
        self.Lsun = _constants[self.version]['L_sun']
        self._constructGrid(dbdtype, lam_min, lam_max)
        self._constructInterpolator()
        return

    def _constructAges(self):
        first_row = _np.loadtxt(self.flist[0], max_rows=1)

        # compute ages that rows correspond to according to BPASS v2.2.1 manual
        self.log_ages = _np.array([
            (6+0.1*(n-2)) for n in range(2, first_row.shape[0]+1)
        ], dtype=self.dbdtype)
        self._aMin = self.log_ages[0]
        self._aMax = self.log_ages[-1]
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
            Interpolated SED [L_sun/angstrom] for each input stellar population
            at all wavelengths provided by the database.
        """
        masses = _check_masses(masses)
        return self.wavelengths, self._interpolate(metallicities, ages)*masses


class BPASSemRatesDatabase(_BPASSdatabase):
    """
    This class models a database for different emissivities as a function of
    metallicity and age of stellar populations for fixed IMF and population
    type.

    Attributes
    ----------
    grid : array
        Array of [Nion, L_alpha, L_FUV, L_NUV] as provided by BPASS as a
        function of metallicity and age.
    """
    def __init__(self, path, version, imf, popType, dbdtype=_np.float64):
        super().__init__(path, version, imf, popType, "ionizing", dbdtype)
        self._constructGrid(dbdtype)
        self._constructInterpolator()
        return

    def _constructAges(self):
        # load first file to extract some info
        fArr = _np.loadtxt(self.flist[0], dtype=self.dbdtype)

        self.log_ages = fArr[:, 0][:]
        self._aMin = self.log_ages[0]
        self._aMax = self.log_ages[-1]
        return

    def _constructGrid(self, dbdtype):
        # load first file to extract some info
        fArr = _np.loadtxt(self.flist[0], dtype=dbdtype)

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
        Interpolate emission rates for stellar populations.

        Computes:
            Nion in 1/s:
                ionizing photon production rate
            L_alpha in ergs/s:
                Balmer H line luminosity, assuming =log(Nion/s)-11.87
            L_FUV in ergs/s/A:
                luminosity in the FUV band (mean flux from 1556 to 1576A)
            L_NUV in ergs/s/A:
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
            Interpolated values as specified above for each input stellar
            population.
        """
        masses = _check_masses(masses)
        return masses*10**(self._interpolate(metallicities, ages))


def _check_masses(masses):
    if _np.amin(masses) < 1e-3:
        _warnings.warn(
            "Input masses below 1000 M_sun! For such small populations,"
            " single stars can contribute a significant fraction of the"
            " population mass and re-scaling BPASS values averaged over"
            " more massive populations likely yields incorrect results."
        )
    try:
        # to be able to multiply array of interpolated values by it
        masses = masses[:, None]
    except TypeError as e:
        if 'not subscriptable' not in str(e):
            raise
    return masses


def bin_spectra(wave, spectra, bins, edges=False):
    """
    Bin spectra conserving luminosity.

    Given SEDs sampled at certain wavelengths/frequencies will compute their
    values in given wavelength/frequency bins. The new SED values are bin
    averages computed using trapezoidal integration. This ensures that the
    luminosity per bin is conserved. Of course, only downsampling really makes
    sense here, i.e. the input SEDs should be well sampled compared to the
    desired output bins.

    Effectively converts input spectra to step functions of
    wavelength/frequency. Note, in particular, that this means that only
    rectangle rule integration can sensibly be performed on the output
    spectra. Higher order integration methods are not meaningful.

    Parameters
    ----------
    wave : array, shape `(N_wave)`
        Wavelengths or frequencies at which spectra are known. Assumed to be
        sorted in ascending order.
    spectra : array, shape `(N, N_wave)`
        The SEDs to resample given as L_lambda [Energy/Time/Wavelength] or L_nu
        [Energy/Time/Frequency] in accordance with `wave`.
    bins : array, shape `(N_bins)`
        The bins to which to resample spectra. Either values in the bins or
        their edges. See `edges`. Assumed to be sorted in ascending
        order. Required to lie within the range provided by `wave`.
    edges : bool, optional
        Whether the values given in `bins` are bin edges or values in the
        bins. If `True`, `N_wave_new=N_bins-1`. If `False`, `N_wave_new=N_bins`
        and in this case bin edges are constructed such that they always lie
        between neighbouring points. The first/last bin is assumed to be
        symmetric around the first/last value in bins.

    Returns
    -------
    wave_new : array, shape `(N_wave_new)`
        The wavelength/frequency values to which spectra were binned. If edges
        is `False`, this will be identical to `bins`. Otherwise it will be the
        bin centers.
    spectra_new : array, shape `(N, N_wave_new)`
        The binned spectra.
    """
    if edges:
        bin_edges = bins
        wave_new = (bins[1:] + bins[:-1])/2
    else:
        wave_new = bins
        bin_edges = _np.empty((len(bins) + 1))
        bin_edges[1:-1] = (bins[1:] + bins[:-1])/2
        bin_edges[0] = bins[0] - (bin_edges[1]-bins[0])
        bin_edges[-1] = bins[-1] + (bins[-1]-bin_edges[-2])
    if not _np.alltrue((bin_edges >= wave[0]) & (bin_edges <= wave[-1])):
        raise ValueError("bin_edges outside of valid range!")
    spectra_new = _np.empty((spectra.shape[0], len(bin_edges) - 1))
    for j, spec in enumerate(spectra):
        _sed_tools.sed_tools.bin_sed(
            wave, spec, bin_edges, spectra_new[j, :]
        )
    return wave_new, spectra_new
