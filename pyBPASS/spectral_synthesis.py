import warnings as _warnings
import glob as _glob
import os as _os
import numpy as _np
from scipy import interpolate as _interpolate
from pyBPASS.database import BPASSdatabase as _BPASSdatabase


class BPASSsedDatabase(_BPASSdatabase):
    """
    Representation and interface for single BPASS SED grid defined by an IMF
    and a population type.

    Attributes
    ----------
    imf : string
        The IMF of this SED grid.
    popType : string
        The population type of this SED grid. Can be single ('sin') or binary
        ('bin').
    metallicities : list
        The metallicities at which BPASS provides SEDs for this population. In
        units of Z_sun.
    log_ages : list
        Log of the population ages [yr] at which BPASS provides SEDs for this
        population.
    wavelengths : array
        The wavelengths [angstrom] at which the SEDs are sampled.
    """

    def __init__(self, path, version, imf, popType):
        super().__init__(path, version)

        # IMF specific subfolder of BPASS data release
        self.folder = \
            'bpass_v'+self.version+'_imf'+imf if not imf[0].isalpha() \
            else 'bpass_v'+self.version+'_imf_'+imf
        if not _os.path.isdir(_os.path.join(self.path, self.folder)):
            raise ValueError(_os.path.join(self.path, self.folder) +
                             "is not a directory.")

        self.imf = imf
        self.popType = popType
        self._constructGrid()
        self._constructInterpolator()
        return

    def _constructGrid(self):
        """
        Load all available SEDs into memory.

        Builds a 3D array of flux in L_sun/angstrom as function of
        (Z,age,lambda) for a simple stellar population of 1e6 M_sun.
        """
        reg = \
            "spectra-"+self.popType+"-imf"+self.imf+"*" \
            if not self.imf[0].isalpha() \
            else "spectra-"+self.popType+"-imf_"+self.imf+"*"
        flist = _glob.glob(_os.path.join(self.path, self.folder, reg))
        flist.sort(
            key=lambda x: self._zFromFilename(_os.path.basename(x))
        )
        self.metallicities = [
            self._zFromFilename(_os.path.basename(x)) for x in flist
        ]
        if not sorted(self.metallicities):
            raise RuntimeError(
                "Sorting available tables by metallicity failed! "
                "Found the following values: "+str(self.metallicities))
        self._zMin = self.metallicities[0]
        self._zMax = self.metallicities[-1]

        fArr = _np.loadtxt(flist[0])
        self.wavelengths = fArr[:, 0]
        self.log_ages = [(6+0.1*(n-2)) for n in range(2, fArr.shape[1]+1)]
        self._aMin = self.log_ages[0]
        self._aMax = self.log_ages[-1]

        self.SEDgrid = _np.zeros(
            (len(self.metallicities),
             len(self.log_ages),
             len(self.wavelengths))
        )
        self.SEDgrid[0, :, :] = fArr[:, 1:].T
        for j, f in enumerate(flist[1:]):
            self.SEDgrid[j+1, :, :] = _np.loadtxt(f)[:, 1:].T
        return

    def _constructInterpolator(self):
        """
        Construction an interpolator on metallicity-age grid. Use
        `LinearNDInterpolator` because it can interpolate vector-valued
        quantities.
        """
        zz, aa = _np.meshgrid(self.metallicities, self.log_ages, indexing='ij')
        self._points = _np.stack((zz, aa), -1).reshape((-1, 2))
        self._interpolator = _interpolate.LinearNDInterpolator(
            self._points, self.SEDgrid.reshape((-1, len(self.wavelengths)))
        )
        return

    def interpolate(self, metallicities, ages):
        """
        Interpolate spectra for stellar populations.

        The units of the returned SEDs are Solar Luminosities per
        Angstrom. Computes SEDs for stellar populations formed in single
        instantaneous bursts. Normalised to a population mass of 1e6 M_sun.

        Interpolation is done in metallicity-log(age) space by triangulation
        and subsequent barycentric interpolation.

        Parameters
        ----------
        metallicities : array, shape `(N)`
            The metallicities of the stellar populations in units of solar
            metallicity.
        ages : array, shape `(N)`
            The ages of the stellar populations [yr].

        Returns
        -------
         : array, shape `(N_lam)`
            Array of wavelengths [angstrom] provided by the database.
         : array, shape `(N, N_lam)`
            Interpolated SED [L_sun/angstrom] normalized to 1e6 M_sun for each
            input stellar population at all wavelengths provided by the
            database.
        """
        # clipping
        if _np.amax(metallicities) > self._zMax or \
           _np.amin(metallicities) < self._zMin:
            _warnings.warn(
                "Input metallicities for spectral synthesis outside of "
                "available range "+str(self._zMin)+", "+str(self._zMax) +
                " provided. They will be clipped."
            )
            metallicities = _np.clip(metallicities, self._zMin, self._zMax)
        if _np.amax(ages) > 10**(self._aMax) or \
           _np.amin(ages) < 10**(self._aMin):
            _warnings.warn(
                "Input ages for spectral synthesis outside of available "
                "range "+str(10**self._aMin)+", "+str(10**self._aMax) +
                " provided. They will be clipped."
            )
            ages = _np.clip(ages, 10**self._aMin, 10**self._aMax)

        return self.wavelengths, \
            self._interpolator(metallicities, _np.log10(ages))


if __name__ == "__main__":
    c = BPASSsedDatabase(
        '/freya/u/mglatzle/data/BPASS/BPASSv2.2.1_release-07-18-Tuatara/',
        '2.2.1',
        'chab300',
        'bin')
