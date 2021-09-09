# The Whole Genome Sequencing Coverage Plot

## Introduction

**The Whole Genome Sequencing Coverage Plot (wgscovplot)** is a tool to generate HTML Interactive Coverage Plot given coverage depth information, variants and DNA Gene features

## Dependencies

### Software Version

- The tool is still in the development phases with current version is ```dev1.0.0```

### Programming Languages

- Python (>=3.6)
- HTML/CSS/Javascript

### Tools/libraries

- [Echarts]
- Python libraries: Biopython, rich, typer, pandas, requests
- [select2]
- [jinja2]

## Installation

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
 -s, --samples-data PATH         List of Sample Names, Coverage, VCF File
                                 [required]
 -o, --output-html PATH          Output Interactive HTML Coverage Plot
                                 [default: coverage_plot.html]
 -r, --ref-seq PATH              Reference genome sequences file  [required]
 -b, --bed PATH                  Bed file
 -g, --genbank PATH              Genbank file contains features of reference
                                 sequence  [required]
 --verbose / --no-verbose        Verbose logs  [default: no-verbose]
 --version / --no-version        Print wgscovplot version and exit
 --install-completion [bash|zsh|fish|powershell|pwsh]
                                 Install completion for the specified shell.
 --show-completion [bash|zsh|fish|powershell|pwsh]
                                 Show completion for the specified shell, to
                                 copy it or customize the installation.
 --help                          Show this message and exit.
```

### Prepare Samples Data

#### The tool currently requires users to prepare samples data file (tsv file) is as the format below. We are updating the tool so that user just simply provide directory and file's identifiers

| Sample Name | Depth file                              | Variant vcf  
|-------------|-----------------------------------------|-------------------------------|
| Sample-01	  | path/to/sample-01-coverage-depth        | path/to/sample-01-variant-vcf | 
| Sample-02	  | path/to/sample-02-coverage-depth        | path/to/sample-02-variant-vcf | 
| Sample-03	  | path/to/sample-03-coverage-depth        | path/to/sample-03-variant-vcf |

#### Note that no header is needed

### Other data:

- Reference Sequences (reference.fasta)
- Genbank file contains gene features (sequence_genbank.gb)

### Output:

The tool will generate the Coverage Plot for samples in HTML file

### Command

```
wgscovplot -s samples.tsv -r reference.fasta -o sequencing_coverage_plot.html -g sequence_genbank.gb
```

## Authors

* Development lead: [Peter Kruczkiewicz]
* Software Developer: [Hai Nguyen]

## License

Copyright 2021 [Canadian Food Inspection Agency of Canada], Government of Canada.

Distributed under the MIT license.

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->


[Peter Kruczkiewicz]: https://github.com/peterk87/
[Hai Nguyen]: https://github.com/nhhaidee/
[Echarts]: https://echarts.apache.org/en/index.html
[select2]: https://select2.org/
[jinja2]: https://jinja.palletsprojects.com/en/3.0.x/
[Canadian Food Inspection Agency of Canada]: https://inspection.canada.ca/science-and-research/our-laboratories/ncfad-winnipeg/eng/1549576575939/1549576643836
