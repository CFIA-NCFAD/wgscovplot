#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'rich',
    'typer',
    'jinja2',
    'pandas',
    'biopython',
    'requests',
    'prettytable',
    'simplejson'
]

setup(
    name="sequencing_coverage_plot",
    version="1.0.0",
    description="A Python library to generate interactive sequencing coverage plots",
    author="Hai Nguyen",
    author_email='nhhaidee@gmail.com',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Software Development :: Libraries",
    ],
    install_requires=requirements,
    packages=find_packages(include=['sequencing_coverage_plot']),
    zip_safe=False,
)