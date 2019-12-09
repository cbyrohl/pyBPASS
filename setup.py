# setuptools import here is required to be able to run tests with `python
# setup.py test`. If I interpret stackoverflow.com/a/41896134/ correctly, this
# also means that setuptools will be used for the actual installation, even
# though we use
# numpy.distutils.setup. [[github.com/wafo-project/pywafo/issues/14][This]]
# issue on github might provide more insights. I consider how it currently
# works good enough, though, so will not continue looking into it.
import setuptools

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

sed_ext = numpy.distutils.extension.Extension(
    name="pyBPASS._sed_tools",
    sources=["fsrc/sed_tools.f90"]
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
    ext_modules=[sed_ext],
)
