#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from setuptools import setup, find_packages
from codecs import open
from os import path, system
from re import compile as re_compile

def read(filename):
    kwds = {"encoding": "utf-8"} if sys.version_info[0] >= 3 else {}
    with open(filename, **kwds) as fp:
        contents = fp.read()
    return contents

# Get the version information.
here = path.abspath(path.dirname(__file__))
version = "0.1"

## TODO require some other things
setup(
    name="caga",
    version=version,
    author="Alex Ji",
    description="Caterpillar-GAMMA modules",
    long_description=read(path.join(here, "README.md")),
    url="https://github.com/caterpillarproject/caga",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    keywords="astronomy",
    packages=find_packages(exclude=["documents", "tests"]),
    install_requires=[
        "numpy",
        ],
)
