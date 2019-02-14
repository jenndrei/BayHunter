#!/usr/bin/env python
try:
    from numpy.distutils.core import Extension
    from numpy.distutils.core import setup
except ImportError:
    raise ImportError('Numpy needs to be installed or updated.')

setup(
    name="BayHunter",
    version="1.0",
    author="Jennifer Dreiling",
    author_email="jennifer.dreiling@gfz-potsdam.de",
    description=("Transdimensional Bayesian Inversion of RF and/or SWD."),
    install_requires=[],
    url="",
    packages=['BayHunter'],
    package_dir={
        'BayHunter': 'src'},

    scripts=['src/scripts/baywatch'],

    package_data={
        'BayHunter': ['defaults/*'], },

    ext_modules=[
        Extension(
            name='BayHunter.surfdisp96_ext',
            sources=['src/surfdisp96.f'],
            extra_f77_compile_args='-O3 -ffixed-line-length-none -fbounds-check -m64'.split(), # noqa
            f2py_options=['only:', 'surfdisp96', ':'],
            language='f77')
        ]
)
