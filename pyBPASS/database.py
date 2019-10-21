import warnings as _warnings
import re as _re
import os as _os
import glob as _glob
import numpy as _np
from scipy import interpolate as _interpolate


class BPASSdatabase(object):
    """
    Parameters
    ----------
    path : string
        Path to BPASS data release.
    version : string
        Version of the BPASS data release.
    imf : string
        Which stellar initial mass function to use.
    popType : string
        The population type to use.

    Attributes
    ----------
    imf : string
        The IMF of this SED grid.
    popType : string
        The population type of this SED grid. Can be single ('sin') or binary
        ('bin').
    metallicities : list
        The metallicities at which BPASS provides data for this IMF and
        population type combination. Given as the mass fraction in metals.
    flist : list
        List of files that belong to this database.
    """

    def __init__(self, path, version, imf, popType, key, dbdtype):
        if not _os.path.isdir(path):
            raise ValueError(path + " is not a directory.")
        self.path = _os.path.abspath(path)
        self.version = version
        self.imf = imf
        self.popType = popType

        # IMF specific subfolder of BPASS data release
        self.folder = \
            'bpass_v'+self.version+'_imf'+imf if not imf[0].isalpha() \
            else 'bpass_v'+self.version+'_imf_'+imf
        if not _os.path.isdir(_os.path.join(self.path, self.folder)):
            raise ValueError(_os.path.join(self.path, self.folder) +
                             " is not a directory.")

        # wildcarded filename to glob for
        reg = \
            key+"-"+self.popType+"-imf"+self.imf+"*" \
            if not self.imf[0].isalpha() \
            else key+"-"+self.popType+"-imf_"+self.imf+"*"

        # list of files that comprise our database
        flist = _glob.glob(_os.path.join(self.path, self.folder, reg))
        flist.sort(
            key=lambda x: self._zFromFilename(_os.path.basename(x))
        )
        self.metallicities = _np.array([
            self._zFromFilename(_os.path.basename(x)) for x in flist
        ], dtype=dbdtype)
        if _np.any(_np.diff(self.metallicities) <= 0):
            raise RuntimeError(
                "Sorting available tables by metallicity failed! "
                "Found the following values: "+str(self.metallicities))
        self.flist = flist
        self._zMin = self.metallicities[0]
        self._zMax = self.metallicities[-1]
        return

    def _zFromFilename(self, fname):
        zStr = _re.search("\\.z(.{3})\\.", fname).group(1)
        try:
            z = float("0."+zStr)
        except ValueError:
            z = float("1e-"+zStr[-1])
        return z

    def _constructInterpolator(self):
        """
        Construct an interpolator on metallicity-age grid. Use
        `LinearNDInterpolator` because it can interpolate vector-valued
        quantities.
        """
        zz, aa = _np.meshgrid(self.metallicities, self.log_ages, indexing='ij')
        self._points = _np.stack((zz, aa), -1).reshape((-1, 2))
        self._interpolator = _interpolate.LinearNDInterpolator(
            self._points, self.grid.reshape((-1, self.grid.shape[2]))
        )
        return

    def _interpolate(self, metallicities, ages, masses=1):
        # argument checking
        try:
            # to be able to multiply array of interpolated values by it
            masses = masses[:, None]
        except TypeError as e:
            if 'not subscriptable' not in str(e):
                raise
        # clipping
        if _np.amax(metallicities) > self._zMax or \
           _np.amin(metallicities) < self._zMin:
            _warnings.warn(
                "Input metallicities for spectral synthesis outside of "
                "available range "+str(self._zMin)+", "+str(self._zMax) +
                " provided. They will be clipped."
            )
            # need to do a copy and cast here
            metallicities = _np.array(
                metallicities, dtype=self.metallicities.dtype
            )
            _np.clip(metallicities, self._zMin, self._zMax, out=metallicities)
        if _np.amax(ages) > 10**(self._aMax) or \
           _np.amin(ages) < 10**(self._aMin):
            _warnings.warn(
                "Input ages for spectral synthesis outside of available "
                "range "+str(10**self._aMin)+", "+str(10**self._aMax) +
                " [yr] provided. They will be clipped."
            )
            # need to do a copy and cast here
            ages = _np.array(
                ages, dtype=self.log_ages.dtype
            )
            _np.clip(ages, 10**self._aMin, 10**self._aMax, out=ages)
        if _np.amin(masses) < 1e-3:
            _warnings.warn(
                "Input masses below 1000 M_sun! For such small populations,"
                " single stars can contribute a significant fraction of the"
                " population mass and re-scaling BPASS values averaged over"
                " more massive populations likely yields incorrect results."
            )
        return self._interpolator(metallicities, _np.log10(ages))*masses
