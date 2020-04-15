#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open("requirements.txt") as requirement_file:
    requirements_packages = requirement_file.readlines()

requirements = requirements_packages
setup_requirements = [ ]

test_requirements = []

setup(
    author="Liting Chen",
    author_email='lchen4@umcutrecht.nl',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="python wrapper with cmd line tool to covert data matrix for use of siCloneFit and pandas",
    entry_points={
        'console_scripts': [
            'siclonefitio=siclonefitio.cli:main',
        ],
    },
    install_requires=requirements,
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='siclonefitio',
    name='siclonefitio',
    packages=find_packages(include=['siclonefitio', 'siclonefitio.*', 'visual', 'visual.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/zztin/siclonefitio',
    version='0.1.0',
    zip_safe=False,
)
