# wgscovplot

[![](https://img.shields.io/pypi/v/wgscovplot.svg)](https://pypi.org/project/wgscovplot/)
[![](https://github.com/nhhaidee/wgscovplot/workflows/CI/badge.svg?branch=master)](https://github.com/nhhaidee/wgscovplot/actions)
[![](https://img.shields.io/badge/License-Apache%20v2.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

**wgscovplot** generates interactive comparative sequencing coverage plots in self-contained, offline-friendly HTML files with optional annotation of variant calling results, PCR amplicon coverage and genetic features.

- [Live example output](https://nhhaidee.github.io)

![](https://raw.githubusercontent.com/nhhaidee/nhhaidee.github.io/master/wgscovplot.png)

## Installation

### From PyPI

Install from PyPI with `pip`

```
pip install wgscovplot
```

If the installation was successful, you should be able to type `wgscovplot --help` and get a help message on how to use the tool.

### Install from source

Clone the `wgscovplot` repository.

```
git clone https://github.com/nhhaidee/wgscovplot.git
```

Then change directory to `wgscovplot` and install.

```
cd wgscovplot
python setup.py install
```

## Features

- Compare sequencing coverage across multiple samples
- Fully-interactive plots with informative tooltips highlighting variant calling results and coverage statistics across all samples being shown  
- Easy-to-use: Simply provide a [nf-core/viralrecon], [peterk87/nf-virontus] Nextflow workflow results directory as input (`wgscovplot --input-dir /path/to/viralrecon/results`) and `wgscovplot` will figure out what files it needs to generate its interactive sequencing coverage plots 
- Annotate coverage plots with variant calling results from multiple different variant callers ([iVar](https://github.com/andersen-lab/ivar), [Nanopolish](https://github.com/jts/nanopolish), [Longshot](https://github.com/pjedge/longshot), [Medaka](https://github.com/nanoporetech/medaka)) and variant effect results from [SnpEff]/[SnpSift]

## Usage

Basic usage will output a `wgscovplot.html` file in the current directory:

```bash
wgscovplot --input-dir /path/to/viralrecon/results
```

Specify an NCBI Accession ID

```bash
wgscovplot \
  --input-dir /path/to/viralrecon/results \
  --ncbi-accession-id MN908947.3
```

Show help info with `$ wgscovplot --help`:

```
Usage: wgscovplot [OPTIONS]

Options:
  --input-dir PATH                Nextflow workflow results directory
                                  [required]
  --output-html PATH              Output File of Interactive HTML Coverage
                                  Plot  [default: wgscovplot.html]
  --ref-seq PATH                  Path to reference sequences
  --genbank PATH                  Genbank file contains gene features
  --ncbi-accession-id TEXT        NCBI accession id to fetch gene features
                                  and/or reference sequences
  --low-coverage-threshold INTEGER
                                  Low Coverage Threshold  [default: 10]
  --amplicon / --no-amplicon      Plot Amplicon Coverage Depth  [default:
                                  amplicon]
  --gene-feature / --no-gene-feature
                                  Plot Gene Features  [default: gene-feature]
  --gene-misc-feature / --no-gene-misc-feature
                                  Plot Miscellaneous Features  [default: no-
                                  gene-misc-feature]
  --dev / --no-dev                Run tool with debug mode  [default: no-dev]
  --verbose / --no-verbose        Verbose logs  [default: no-verbose]
  --version / --no-version        Print wgscovplot version and exit
  --install-completion [bash|zsh|fish|powershell|pwsh]
                                  Install completion for the specified shell.
  --show-completion [bash|zsh|fish|powershell|pwsh]
                                  Show completion for the specified shell, to
                                  copy it or customize the installation.
  --help                          Show this message and exit.

```

## Dependencies

- [Python](https://www.python.org/) (>=3.8)
    - [BioPython](https://github.com/biopython/biopython/)
    - [rich](https://rich.readthedocs.io/)
    - [typer](https://github.com/tiangolo/typer)
    - [Pandas](https://pandas.pydata.org/)
    - [requests](https://docs.python-requests.org/)
    - [jinja2]
- Javascript
    - [npm](https://www.npmjs.com/) or [yarn](https://yarnpkg.com/) for building JS bundle
    - [Echarts] for generating interactive plots
    - [select2] for nice select boxes
    - [Bootstrap](https://getbootstrap.com/) for styling

## Authors

* Development Lead: [Peter Kruczkiewicz]
* Software Developer: [Hai Nguyen]

## License

Copyright 2021 Canadian Food Inspection Agency of Canada, Government of Canada.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

[Peter Kruczkiewicz]: https://github.com/peterk87/
[Hai Nguyen]: https://github.com/nhhaidee/
[Echarts]: https://echarts.apache.org/en/index.html
[select2]: https://select2.org/
[jinja2]: https://jinja.palletsprojects.com/en/3.0.x/
[SnpEff]: https://pcingola.github.io/SnpEff/se_introduction/
[SnpSift]: https://pcingola.github.io/SnpEff/ss_introduction/
[Mosdepth]: https://github.com/brentp/mosdepth
[nf-core/viralrecon]: https://github.com/nf-core/viralrecon
[peterk87/nf-virontus]: https://github.com/peterk87/nf-virontus/
[Canadian Food Inspection Agency of Canada]: https://inspection.canada.ca/science-and-research/our-laboratories/ncfad-winnipeg/eng/1549576575939/1549576643836