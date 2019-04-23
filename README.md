Tools for interaction with BPASS data releases. Implemented on DR 2.2.1 so they
might not work for other versions.

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


wavlengths, seds = db.interpolate(metallicities, ages, masses)
```

For more info see docstrings.
