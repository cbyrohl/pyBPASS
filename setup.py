try:
    import numpy.distutils.core
except ImportError as e:
    raise ImportError(
        "Import of numpy.distutils failed, can not proceed with setup."
        " Please make sure you have a current and functioning version"
        " of numpy installed.\n"
        " Caught the following exception: " + str(e)
    )
import numpy.distutils.extension

rebin_ext = numpy.distutils.extension.Extension(
    name="pyBPASS._rebin",
    sources=["fsrc/rebin.f90"]
)

with open('README.md') as fh:
    long_desc = fh.read()

numpy.distutils.core.setup(
    name='pyBPASS',
    version='1.0.1',
    author_email='mglatzle@mpa-garching.mpg.de',
    description='Python tools for BPASS data releases.',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    license='GPL-3.0',
    packages=['pyBPASS'],
    install_requires=[
        'numpy',
        'scipy',
    ],
    tests_require=[
        'pytest'
    ],
    ext_modules=[rebin_ext],
)
