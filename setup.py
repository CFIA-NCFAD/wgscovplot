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
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    developer="Hai Nguyen",
    developer_email='nhhiadee@gmail.com',
    development_lead="Peter Kruczkiewicz",
    development_lead_email="peter.kruczkiewicz@gmail.com",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Generate a standalone HTML file with an interactive coverage plot",
    entry_points={
        'console_scripts': [
            'shicp=coverage_plot.cli:app',
        ],
    },
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='shicp',
    name='shicp',
    packages=find_packages(
        where='sequencing_coverage_plot',
        include=['coverage_plot']
    ),
    package_dir={"": "sequencing_coverage_plot"},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/nhhaidee/sequencing_coverage_plot',
    version='1.0.0dev',
    zip_safe=False,
)
