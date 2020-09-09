import abc as _abc
import warnings as _warnings
import re as _re
import os as _os
import glob as _glob
import numpy as _np
from scipy import interpolate as _interpolate


class BPASSdatabase(_abc.ABC):
    """
    This class models a database that contains a certain property of stellar
    populations as a function of age and metallicity. A database instance for a
    fixed BPASS version is uniquely defined by a property (e.g. 'spectra'), an
    IMF and a population type (binary or not). This class is abstract and
    intended to be a base class. Specialized subclasses should be used for each
    property.

    Attributes
    ----------
    version : string
        BPASS version this database uses.
    imf : string
        The initial mass function (IMF) this database uses.
    popType : string
        The population type uses.
    metallicities : array
        The metallicities at which BPASS provides data for this IMF and
        population type combination. Given as the mass fraction in metals.
    log_ages : array
        log10 of the ages [yr] at which BPASS provides data for this IMF and
        population type combination.
    flist : list of strings
        List of files that belong to this database.
    dbdtype : dtype
        Data type used by this database.
    """

    @_abc.abstractmethod
    def __init__(self, basepath, version, imf, popType, key,
                 dbdtype=_np.float64):
        """
        Parameters
        ----------
        basepath : string
            Path to BPASS data release.
        version : string
            Version of the BPASS data release.
        imf : string
            Which stellar initial mass function to use.
        popType : string
            The population type to use.
        key : string
            The stellar population property this database contains.
        dbdtype : dtype, optional
            Data type to use for arrays read from disk. Defaults to numpy's
            float64.
        """
        if not _os.path.isdir(basepath):
            raise ValueError(basepath + " is not a directory.")
        self.version = version
        self.imf = imf
        self.popType = popType
        self.dbdtype = dbdtype

        self._constructFlist(basepath, key)
        self._constructMetallicities()
        self._constructAges()

        return

    def _constructFlist(self, basepath, key):
        # folder containing files for given parameter configuration
        folder = _os.path.join(
            basepath,
            'bpass_v'+self.version+'_imf'+self.imf if not self.imf[0].isalpha()
            else 'bpass_v'+self.version+'_imf_'+self.imf
        )
        if not _os.path.isdir(folder):
            raise ValueError(f"{folder} is not a directory.")
        # wildcarded filename to glob for
        reg = \
            key+"-"+self.popType+"-imf"+self.imf+"*" \
            if not self.imf[0].isalpha() \
            else key+"-"+self.popType+"-imf_"+self.imf+"*"

        # list of files that comprise our database
        flist = _glob.glob(_os.path.join(folder, reg))
        flist.sort(
            key=lambda x: self._zFromFilename(_os.path.basename(x))
        )
        self.flist = flist
        return

    def _constructMetallicities(self):
        self.metallicities = _np.array([
            self._zFromFilename(_os.path.basename(x)) for x in self.flist
        ], dtype=self.dbdtype)
        if _np.any(_np.diff(self.metallicities) <= 0):
            raise RuntimeError(
                "Sorting available files by metallicity failed! "
                "Found the following values: "+str(self.metallicities) +
                " and the following files: "+str(self.flist)
            )
        self._zMin = self.metallicities[0]
        self._zMax = self.metallicities[-1]
        return

    @_abc.abstractmethod
    def _constructAges(self):
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

    def _interpolate(self, metallicities, ages):
        # clipping
        if _np.amax(metallicities) > self._zMax or \
           _np.amin(metallicities) < self._zMin:
            _warnings.warn(
                "Input metallicities for interpolation outside of "
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
                "Input ages for interpolation outside of available "
                "range "+str(10**self._aMin)+", "+str(10**self._aMax) +
                " [yr] provided. They will be clipped."
            )
            # need to do a copy and cast here
            ages = _np.array(
                ages, dtype=self.log_ages.dtype
            )
            _np.clip(ages, 10**self._aMin, 10**self._aMax, out=ages)
        return self._interpolator(metallicities, _np.log10(ages))
