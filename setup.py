#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, glob, sys
import numpy

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils.extension import Extension
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
ext_modules = []

# make sure we're using GCC
os.environ["CC"] = "gcc"

if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')

#------------------------------------------------------------------------------
# Compiling the C code for the Snobal libary

loc = 'pysnobal/libsnobal'
cwd = os.getcwd()
sources = glob.glob(os.path.join(loc, '*.c'))
loc = 'pysnobal'
sources += [os.path.join(loc, val) for val in ["snobal.pyx"]]
ext_modules += [
                Extension(
                    "pysnobal.snobal",
                    sources,
                    #libraries=["snobal"],
                    include_dirs=[numpy.get_include(),'pysnobal', 'pysnobal/h'],
                    #runtime_library_dirs=['{}'.format(os.path.join(cwd,'pysnobal'))],
                    extra_compile_args=['-fopenmp', '-O3', '-L./pysnobal'],
                    extra_link_args=['-fopenmp', '-O3', '-L./pysnobal']
                    )
                ]

cmdclass.update({ 'build_ext': build_ext })

setup(
    name='pysnobal',
    version='0.2.0',
    description="Python wrapper of the Snobal point model",
    long_description=readme + '\n\n' + history,
    author="Scott Havens",
    author_email='scott.havens@ars.usda.gov',
    url='https://gitlab.com/ars-snow/pysnobal',
    # packages=[
    #     'pysnobal', 'pysnobal.libsnobal'
    # ],
    packages=[
        'pysnobal'
    ],
#     package_dir={'pysnobal':
#                  'pysnobal'},
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
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
    ],
    test_suite='tests',
    tests_require=test_requirements,
    cmdclass = cmdclass,
    ext_modules = ext_modules,
)
