#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import sys
import numpy

from setuptools import setup, find_packages, Extension

# from distutils.extension import Extension
from Cython.Distutils import build_ext

with open('README.md') as readme_file:
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

# make sure we're using GCC
if "CC" not in os.environ:
    os.environ["CC"] = "gcc"

if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')

# ------------------------------------------------------------------------------
# Compiling the C code for the Snobal library

loc = 'pysnobal/c_snobal/libsnobal'
cwd = os.getcwd()
sources = glob.glob(os.path.join(loc, '*.c'))

loc = 'pysnobal/c_snobal'
extra_cc_args = ['-fopenmp', '-O3', '-L./pysnobal']
sources += [os.path.join(loc, val) for val in ["snobal.pyx"]]
extensions = [
    Extension(
        "pysnobal.c_snobal.snobal",
        sources,
        # libraries=["snobal"],
        include_dirs=[
            numpy.get_include(),
            'pysnobal/c_snobal',
            'pysnobal/c_snobal/h'
        ],
        # runtime_library_dirs=['{}'.format(os.path.join(cwd,'pysnobal'))],
        extra_compile_args=extra_cc_args,
        extra_link_args=extra_cc_args,
    )
]

cmdclass.update({'build_ext': build_ext})

setup(
    name='pysnobal',
    version='0.2.0',
    description="Python wrapper of the Snobal mass and "
                "energy balance snow model",
    long_description=readme + '\n\n' + history,
    author="Scott Havens",
    author_email='scott.havens@ars.usda.gov',
    url='https://github.com/USDA-ARS-NWRC/pysnobal',
    packages=find_packages(exclude=["tests"]),
    include_package_data=True,
    install_requires=requirements,
    license="CC0 1.0",
    zip_safe=False,
    keywords='pysnobal',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: CC0 1.0',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
    ],
    test_suite='tests',
    tests_require=test_requirements,
    cmdclass=cmdclass,
    ext_modules=extensions,
)
