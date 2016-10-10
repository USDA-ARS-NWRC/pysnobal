#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, glob

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils.extension import Extension
from Cython.Distutils import build_ext

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]


cmdclass = {}
ext_modules = []

# make sure we're using GCC
os.environ["CC"] = "gcc"

loc = 'pysnobal_c/libsnobal'
# sources=[os.path.join(loc, val) for val in ["_adj_layers.c"]]
sources = glob.glob(os.path.join(loc, '*.c'))
# sources = ['_adj_layers.c']
ext_modules += [
                Extension(
                    "pysnobal_c.libsnobal",
                    sources,
                    include_dirs=["pysnobal_c/h"],
                    extra_compile_args=['-fopenmp', '-O3'],
                    extra_link_args=['-fopenmp', '-O3'],
                    )
                ]

cmdclass.update({ 'build_ext': build_ext })

setup(
    name='pysnobal_c',
    version='0.1.0',
    description="Python wrapper of the Snobal point model",
    long_description=readme + '\n\n' + history,
    author="Scott Havens",
    author_email='scott.havens@ars.usda.gov',
    url='https://gitlab.com/ars-snow/pysnobal',
    packages=[
        'pysnobal_c', 'pysnobal_c.libsnobal'
    ],
#     package_dir={'pysnobal':
#                  'pysnobal'},
    include_package_data=True,
    install_requires=requirements,
    license="ISCL",
    zip_safe=False,
    keywords='pysnobal_c',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    cmdclass = cmdclass,
    ext_modules = ext_modules,
)