=================================================
Standalone HTML Interactive Coverage Plot (SHICP)
=================================================


Generate a standalone HTML file with an interactive coverage plot using Echarts_


* Free software: Apache Software License 2.0
* Documentation: https://shicp.readthedocs.io.


Features
--------

* Interactively view Coverage Plot using the power of Echarts_.
* DNA gene features are included with Coverage Plot which allow user to view, zoom in/out coverage depth for each specific region or the whole genome.
* A collapsible control menu allows to select sample to be displayed, set data zoom (start, end) for viewing region range
* Dynamic coverage statistic is also displayed in tooltips, providing rich information.

Usage
-----

Help output

.. code-block::

    $shicp --help
    Usage: shicp [OPTIONS]

    Options:
      -s, --samples-data PATH         List of Sample Names, Coverage, VCF File
                                      [required]

      -o, --output-html PATH          Output Interactive HTML Coverage Plot
                                      [default: coverage_plot.html]

      -r, --ref-seq PATH              Reference genome sequences file  [required]
      -b, --bed PATH                  Bed file
      -g, --genbank PATH              Genbank file contains features of reference
                                      sequence  [required]

      --verbose / --no-verbose        Verbose logs  [default: False]
      --version / --no-version        Print "shicp version" and exit
      --install-completion [bash|zsh|fish|powershell|pwsh]
                                      Install completion for the specified shell.
      --show-completion [bash|zsh|fish|powershell|pwsh]
                                      Show completion for the specified shell, to
                                      copy it or customize the installation.

      --help                          Show this message and exit.



*See an example shicp output here:*

- `sequencing_coverage_plot.html`_

.. _`sequencing_coverage_plot.html`: docs/data/sequencing_coverage_plot.html
.. _Echarts: https://echarts.apache.org/en/index.html