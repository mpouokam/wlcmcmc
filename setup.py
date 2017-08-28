#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import re

from setuptools import setup, find_packages

from distutils.core import Extension

import numpy

classifiers = [
    'Development Status :: 2 - Pre-Alpha',
    'Programming Language :: Python',
    'Natural Language :: English',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'License :: Other/Proprietary License',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: Implementation :: CPython',
    'Private :: Do Not Upload',  # This prevents accidental uploads to pypi
]

# When in doubt...  pip install -e .[all]
requirements = {
    'requirements': [
    ],

    'setup': [
        'pip',
    ],

    'tests': [
        'pytest==3.2.0'
    ],
}

# Developers should probably run:  pip install .[dev]
#   Overrides and expands dev above
requirements['dev'] = [
    r for k, reqs in requirements.items() for r in reqs
    if k not in ['requirements']
]

# All is for usability:  pip install .[all]
requirements['all'] = [
    r for k, reqs in requirements.items() for r in reqs
]

# Find package files
packages = find_packages()

extras = {k: v for k, v in requirements.items() if k != 'requirements'}

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

setup(
    name="wlcmcmc",
    version="0.1.0",
    description="Worm-Like Chain Generation (Markov Chain Monte Carlo)",
    long_description="Uses Markov Chain Monte Carlo to generate Worm-Like Chains",
    url="https://github.com/rhymeswithlion/wlcmcmc",
    license='',
    author="Brian Cruz",
    author_email="brian@math.berkeley.edu",

    # Package Properties
    packages=packages,
    include_package_data=True,
    install_requires=requirements['requirements'],
    extras_require=extras,

    zip_safe=False,
    keywords="wlcmcmc",
    classifiers=classifiers,
    test_suite='tests',
    setup_requires=requirements.get('setup') or [],
    tests_require=requirements.get('tests') or [],

    # Extension modules
    ext_modules=[
        Extension(b"wlcmcmc/_wlcmcmclib",
                  [b"wlcmcmc/wlcmcmclib.i",
                   b"wlcmcmc/wlcmcmclib.cpp",
                   b"wlcmcmc/mtrand/mtrand.cpp"
                   ],
                  include_dirs=[numpy_include],
                  swig_opts=['-c++'],
                  )
    ]
)
