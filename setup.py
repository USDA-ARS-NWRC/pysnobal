#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
import sys

import numpy
# from distutils.extension import Extension
from Cython.Distutils import build_ext
from setuptools import Extension, find_packages, setup

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
os.environ["CC"] = "gcc"

if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')

# ------------------------------------------------------------------------------
# Compiling the C code for the Snobal libary

# the C code
loc = 'pysnobal/c_snobal/libsnobal'
cwd = os.getcwd()
sources = glob.glob(os.path.join(loc, '*.c'))

# # compile c_snobal.pyx
# loc = 'pysnobal/c_snobal'
# c_snobal_source = [os.path.join(loc, val) for val in ["c_snobal.pyx"]]
# extensions = [
#     Extension(
#         "pysnobal.c_snobal.c_snobal",
#         sources + c_snobal_source,
#         # libraries=["snobal"],
#         include_dirs=[numpy.get_include(), 'pysnobal/c_snobal',
#                       'pysnobal/c_snobal/h'],
#         # runtime_library_dirs=['{}'.format(os.path.join(cwd,'pysnobal'))],
#         extra_compile_args=['-fopenmp', '-O3', '-L./pysnobal'],
#         extra_link_args=['-fopenmp', '-O3', '-L./pysnobal']
#     )
# ]

# compile cython_snobal.pyx
loc = 'pysnobal/c_snobal'
cython_snobal_source = [os.path.join(loc, val)
                        for val in ["cython_snobal.pyx"]]
extensions = [
    Extension(
        "pysnobal.c_snobal.cython_snobal",
        sources + cython_snobal_source,
        # libraries=["snobal"],
        include_dirs=[numpy.get_include(), 'pysnobal/c_snobal',
                      'pysnobal/c_snobal/h'],
        # runtime_library_dirs=['{}'.format(os.path.join(cwd,'pysnobal'))],
        extra_compile_args=['-fopenmp', '-O3', '-L./pysnobal'],
        extra_link_args=['-fopenmp', '-O3', '-L./pysnobal']
    )
]

cmdclass.update({'build_ext': build_ext})

setup(
    name='pysnobal',
    version='0.2.0',
    description="Python wrapper of the Snobal mass and energy balance snow model",
    long_description=readme + '\n\n' + history,
    author="Scott Havens",
    author_email='scott.havens@ars.usda.gov',
    url='https://github.com/USDA-ARS-NWRC/pysnobal',
    packages=find_packages(exclude=["tests"]),
    include_package_data=True,
    package_data={'pysnobal': ['./pysnobal_core_config.ini']},
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
