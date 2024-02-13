# wgscovplot

[![](https://img.shields.io/pypi/v/wgscovplot.svg)](https://pypi.org/project/wgscovplot/)
[![](https://github.com/nhhaidee/wgscovplot/workflows/CI/badge.svg?branch=master)](https://github.com/nhhaidee/wgscovplot/actions)
[![](https://img.shields.io/badge/License-Apache%20v2.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

**wgscovplot** generates interactive comparative sequencing coverage plots in self-contained, offline-friendly HTML files with optional annotation of variant calling results, PCR amplicon coverage and genetic features.

[//]: # (TODO: Add screenshots of wgscovplot output)
[//]: # (TODO: Add updated example output for both segmented and non-segmented viruses)
- [Live example output](https://nhhaidee.github.io)

## Installation

### PyPI

Install from PyPI with `pip`

```bash
pip install wgscovplot
```

If the installation was successful, you should be able to type `wgscovplot --help` and get a help message on how to use the tool.

### Source

Use `pip` to install from source:

```bash
# optionally, create a virtual environment
python -m venv venv
source venv/bin/activate
# install from GitHub repo
pip install git+https://github.com/nhhaidee/wgscovplot.git
# run wgscovplot
wgscovplot --help
wgscovplot /path/to/results_folder
```

## Features

- Easy-to-use: Simply provide a Nextflow output directory containing  and `wgscovplot` will figure out what files it needs to generate its interactive sequencing coverage plots 
  - Compatible workflows: 
    - [nf-core/viralrecon]
    - [CFIA-NCFAD/nf-virontus] 
    - [CFIA-NCFAD/nf-flu]
- Fully-interactive plots featuring:
  - Zoom, scroll, pan, select regions of interest 
  - Informative tooltips highlighting variant calling results and coverage statistics across all samples being shown
  - Change the y-axis scale to linear or log scale
  - Select which samples to show (and which Influenza gene segments to show)
  - Highlight regions of interest (e.g. genetic features, primer/probe binding sites, low coverage regions)
- Annotate coverage plots with variant calling results from multiple different variant callers and variant effect results from [SnpEff]/[SnpSift]
  - Supported variant callers: [iVar](https://github.com/andersen-lab/ivar), [Nanopolish](https://github.com/jts/nanopolish), [Longshot](https://github.com/pjedge/longshot), [Medaka](https://github.com/nanoporetech/medaka), [Clair3](https://github.com/HKU-BAL/Clair3) 
- Compare sequencing coverage across multiple samples

## Usage

Basic usage will output a `wgscovplot.html` file in the current directory:

```bash
wgscovplot /path/to/results_folder
```

Show help info with `$ wgscovplot --help`:

```
 Usage: wgscovplot [OPTIONS] INPUT_DIR

╭─ Arguments ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *    input_dir      PATH  Directory containing Mosdepth and variant calling results from sequence analysis. For example, the output directory from execution of the nf-core/viralrecon or CFIA-NCFAD/nf-flu Nextflow workflow                │
│                           [default: None]                                                                                                                                                                                                    │
│                           [required]                                                                                                                                                                                                         │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --output-html             -o                          PATH                             wgscovplot HTML output file [default: wgscovplot.html]                                                                                                │
│ --primers-fasta           -p                          PATH                             FASTA file containing real-time PCR primer/probe sequences. [default: None]                                                                           │
│ --low-coverage-threshold  -l                          INTEGER                          Low sequencing coverage threshold. [default: 10]                                                                                                      │
│ --edit-distance           -d                          INTEGER                          The maximum differences or 'edits' allowed between real-time PCR primer/probe sequences and the sample sequences. [default: 0]                        │
│ --compress-depths             --no-compress-depths                                     Compress coverage depth arrays? [default: compress-depths]                                                                                            │
│ --verbose                 -v                                                           Verbose logs                                                                                                                                          │
│ --force                   -f                                                           Force overwrite of existing output files                                                                                                              │
│ --version                     --no-version                                             Print wgscovplot version and exit [default: no-version]                                                                                               │
│ --install-completion                                  [bash|zsh|fish|powershell|pwsh]  Install completion for the specified shell. [default: None]                                                                                           │
│ --show-completion                                     [bash|zsh|fish|powershell|pwsh]  Show completion for the specified shell, to copy it or customize the installation. [default: None]                                                    │
│ --help                                                                                 Show this message and exit.                                                                                                                           │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

 wgscovplot version 0.3.0; Python 3.11.6
```

## Dependencies

- [Python](https://www.python.org/) (>=3.9)
    - [BioPython](https://github.com/biopython/biopython/) for all things bioinformatics
    - [rich](https://rich.readthedocs.io/) for pretty terminal output
    - [typer](https://github.com/tiangolo/typer) for CLI
    - [Pandas](https://pandas.pydata.org/) for data wrangling
    - [jinja2] for HTML templating
    - [Edlib](https://github.com/Martinsos/edlib) for fuzzy string matching
- Typescript/Javascript
  - [Echarts] for performant generating interactive plots
  - [SolidJS](https://www.solidjs.com/) for reactive UI components
  - [Vite](https://vitejs.dev/) for TS/JS dev and building bundle
  - [Tailwind CSS](https://tailwindcss.com/) for styling

## Development

This project has two main components: a Python "backend" (CLI that spits out a templated HTML with built JS embedded) and a Javascript frontend.

Python backend development is done in the `wgscovplot` directory. 

Web frontend development is done in the `web` directory. The frontend is built with [Vite](https://vitejs.dev/), [SolidJS](https://www.solidjs.com/) and [ECharts]. 

### Environment

Python development is recommended with PyCharm. Jetbrains IDEs have great support for Python development and virtual environments.

Jetbrains IDEs work great for Typescript/Javascript development as well, but any editor will do if you have Vite live reload enabled. 

### Setup

```bash
# clone repo
git clone https://github.com/nhhaidee/wgscovplot.git
cd wgscovplot
# optionally, create a virtual environment
python -m venv venv
source venv/bin/activate
# install dev dependencies
pip install hatch

# start shell with Hatch
hatch shell

# run linting with Hatch
hatch run lint:all

# run tests with Hatch
hatch run cov
```

### Frontend development

See [web/README.md](web/README.md) for more details.

## Authors

* Development Lead: [Peter Kruczkiewicz]
* Software Developer: [Hai Nguyen]

## License

Copyright 2024 Canadian Food Inspection Agency of Canada, Government of Canada.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

[Peter Kruczkiewicz]: https://github.com/peterk87/
[Hai Nguyen]: https://github.com/nhhaidee/
[Echarts]: https://echarts.apache.org/en/index.html
[select2]: https://select2.org/
[jinja2]: https://jinja.palletsprojects.com/
[SnpEff]: https://pcingola.github.io/SnpEff/se_introduction/
[SnpSift]: https://pcingola.github.io/SnpEff/ss_introduction/
[Mosdepth]: https://github.com/brentp/mosdepth
[nf-core/viralrecon]: https://github.com/nf-core/viralrecon
[peterk87/nf-virontus]: https://github.com/peterk87/nf-virontus/
[CFIA-NCFAD/nf-flu]: https://github.com/CFIA-NCFAD/nf-flu/
[Canadian Food Inspection Agency of Canada]: https://inspection.canada.ca/science-and-research/our-laboratories/ncfad-winnipeg/eng/1549576575939/1549576643836