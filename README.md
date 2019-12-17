Tools for interaction with [BPASS](https://bpass.auckland.ac.nz/9.html) data
releases. Implemented on [DR
2.2.1](https://drive.google.com/drive/folders/1BS2w9hpdaJeul6-YtZum--F4gxWIPYXl)
so they might have to be slightly modified for other versions.

# Spectral synthesis
To compute SEDs for stellar populations the mass, metallicity and age of which
are known:

```python
from pyBPASS.spectral_synthesis import BPASSsedDatabase

path = 'path-to-BPASS-data-release'
version = '2.2.1'
imf = 'chab300'
popType = 'bin'

# takes some time since it loads complete SED grid from disk
db = BPASSsedDatabase(path, version, imf, popType)


wavelengths, seds = db.interpolate(metallicities, ages, masses)
```

Similarly, to compute the emission rates of photons in different bands for the above
stellar populations:
```python
from pyBPASS.spectral_synthesis import BPASSionRatesDatabase

db = BPASSionRatesDatabase(path, version, imf, popType)

Nion = db.interpolate(metallicities, ages, masses)
```

To bin spectra (e.g. for use in a radiative transfer code):
```python
from pyBPASS.spectral_synthesis import bin_spectra

wave_new, seds_new = bin_spectra(wavelengths, seds, bins)
```

For more information, see docstrings. They should be pretty clear.

# Installation
`numpy` is required for the installation process itself, not only as
dependency. So please make sure it is installed (`pip install numpy`) before
proceeding. As far as I know, there are only hacky ways of bootstrapping
`numpy` during setup so I decided against implementing any of them.

To install:
```sh
git clone git@gitlab.mpcdf.mpg.de:mglatzle/pybpass.git
cd pyBPASS
pip install .
```
If you want to install in `develop` mode instead:
```sh
pip install -e .
```

# Tests
Tests can be run by `python setup.py test`. Most of them require a local
database copy. The path can be set in `pyBPASS.tests.config`.
