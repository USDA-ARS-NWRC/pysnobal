#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read()

setup_requirements = ['setuptools_scm', ]

test_requirements = ['pytest>=3', 'pytest-runner']

setup(
    author="USDA ARS Northwest Watershed Research Center",
    author_email='snow@ars.usda.gov',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Take 2 for pysnobal in pure python",
    entry_points={
        'console_scripts': [
            'pysnobal=pysnobal.cli:main',
        ],
    },
    install_requires=requirements,
    license="CC0 1.0",
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords='pysnobal',
    name='pysnobal',
    packages=find_packages(include=['pysnobal', 'pysnobal.*']),
    package_data={
        'pysnobal': [
            './pysnobal_core_config.ini'
        ]
    },
    use_scm_version={
        'local_scheme': 'node-and-date',
    },
    setup_requires=setup_requirements,
    test_suite='pysnobal.tests',
    tests_require=test_requirements,
    url='https://github.com/scotthavens/pysnobal',
    zip_safe=False,
)
