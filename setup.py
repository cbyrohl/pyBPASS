import setuptools

with open('README.md') as fh:
    long_desc = fh.read()

setuptools.setup(
    name='pyBPASS',
    version='1.0.0',
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
)
