# The Whole Genome Sequencing Coverage Plot

## Introduction

**The Whole Genome Sequencing Coverage Plot (wgscovplot)** is a tool to generate HTML Interactive Coverage Plot given coverage depth information, variants and DNA Gene features

## Dependencies

### Software Version

- The tool is still in the development phases with current version is ```1.0.0.dev0```

### Limitation

- In the development phase, the tool is designed for SARS-COV2 virus. The supports for other virus such as Avian Influenza, African Swine Fever Virus (ASFV), Foot and Mouth Disease (FMD) will be added soon.

### Programming Languages

- Python (>=3.7)
- HTML/CSS/Javascript

### Tools/libraries

- [Echarts]
- Python libraries: Biopython, rich, typer, pandas, requests
- [select2]
- [jinja2]

## Installation

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

### Install by pip

```
pip install wgscovplot
```

### Install by conda

```
conda install -c bioconda wgscovplot
```

If the installation was successful, you should be able to type `wgscovplot --help` and get a help message on how to use the tool.


## Usage

```
Usage: wgscovplot [OPTIONS]

Options:
  --input-dir PATH                Nextflow workflow results directory
                                  [required]
  --output-html PATH              Output File of Interactive HTML Coverage
                                  Plot  [default: wgscovplot.html]
  --ref-seq PATH                  Reference Sequences  [required]
  --genbank PATH                  Genbank file contains features of reference
                                  sequence
  --amplicon / --no-amplicon      Plot Amplicon Coverage Depth  [default: no-
                                  amplicon]
  --gene-feature / --no-gene-feature
                                  Plot Gene Feature  [default: no-gene-
                                  feature]
  --verbose / --no-verbose        Verbose logs  [default: no-verbose]
  --version / --no-version        Print wgscovplot version and exit
  --install-completion [bash|zsh|fish|powershell|pwsh]
                                  Install completion for the specified shell.
  --show-completion [bash|zsh|fish|powershell|pwsh]
                                  Show completion for the specified shell, to
                                  copy it or customize the installation.
  --help                          Show this message and exit.

```

### Other data:

- Reference Sequences (reference.fasta)
- Genbank file contains gene features (sequence_genbank.gb)

### Output:

The tool will generate the Coverage Plot for samples in HTML file

### Command

```
wgscovplot --input-dir /path/to/nextflow_results_folder --ref-seq reference.fasta --genbank sequence_genbank.gb --gene-feature
```
#### For amplicon coverage plot
```
wgscovplot --input-dir /path/to/nextflow_results_folder --ref-seq reference.fasta --genbank sequence_genbank.gb --gene-feature --amplicon
```

## Authors

* Development Lead: [Peter Kruczkiewicz]
* Software Developer: [Hai Nguyen]

## License

Copyright 2021 Canadian Food Inspection Agency of Canada, Government of Canada.

Distributed under the MIT license.

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->


[Peter Kruczkiewicz]: https://github.com/peterk87/
[Hai Nguyen]: https://github.com/nhhaidee/
[Echarts]: https://echarts.apache.org/en/index.html
[select2]: https://select2.org/
[jinja2]: https://jinja.palletsprojects.com/en/3.0.x/
[Canadian Food Inspection Agency of Canada]: https://inspection.canada.ca/science-and-research/our-laboratories/ncfad-winnipeg/eng/1549576575939/1549576643836