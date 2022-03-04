#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'typer',
    'rich',
    'jinja2',
    'pandas',
    'numpy',
    'biopython',
    'requests',
    'pydantic',
    'markdown',
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="Generate a standalone HTML file with an interactive coverage plot",
    entry_points={
        'console_scripts': [
            'wgscovplot=wgscovplot.cli:app',
        ],
    },
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='wgscovplot',
    name='wgscovplot',
    packages=find_packages(include=['wgscovplot', 'wgscovplot.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/nhhaidee/sequencing_coverage_plot',
    version='0.1.0',
    zip_safe=False,
)
