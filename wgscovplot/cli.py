# -*- coding: utf-8 -*-

"""Console script."""
import logging
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional
from rich.logging import RichHandler
from sys import version_info
from Bio.SeqIO.FastaIO import SimpleFastaParser
from .wgscovplot import (write_html_coverage_plot,
                         get_variant_data, get_depth_data, get_gene_feature, get_region_amplicon,
                         get_depth_amplicon, get_coverage_stat)
from wgscovplot import __version__
import markdown

app = typer.Typer()


def version_callback(value: bool):
    if value:
        typer.echo(f"wgscovplot version {__version__}")
        raise typer.Exit()


@app.command(
    epilog=f"wgscovplot version {__version__}; Python {version_info.major}.{version_info.minor}.{version_info.micro}"
)
def main(
        samples_data: Path = typer.Option(..., "-s", "--samples-data", help="List of Sample Names, Coverage, VCF File"),
        output_html: Path = typer.Option("coverage_plot.html",
                                         "-o", "--output-html",
                                         help="Output Interactive HTML Coverage Plot"),
        ref_seq: Path = typer.Option(...,
                                     "-r",
                                     "--ref-seq",
                                     help="Reference genome sequences file"),
        genbank: Path = typer.Option(...,
                                     "-g",
                                     "--genbank",
                                     help="Genbank file contains features of reference sequence"),
        amplicon: bool = typer.Option(False, help="Amplicon Coverage Plot"),
        verbose: bool = typer.Option(False, help="Verbose logs"),
        version: Optional[bool] = typer.Option(None,
                                               callback=version_callback,
                                               help=f'Print {"wgscovplot version"} and exit')
):
    from rich.traceback import install

    install(show_locals=True)

    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.INFO if not verbose else logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )
    # First read reference sequence
    if ref_seq:
        with open(ref_seq) as fh:
            for name, seq in SimpleFastaParser(fh):
                ref_seq = seq
    # Get DNA Gene features
    gene_feature = get_gene_feature(genbank)

    # Main Script
    if amplicon:
        df_samples_amplicon = pd.read_table(samples_data,
                                            names=['amplicon_region_file', 'amplicon_perbase_file', 'vcf_file'],
                                            index_col=0, header=None)
        df_samples_amplicon = df_samples_amplicon.fillna(0)
        # Get list of samples name
        samples_name = df_samples_amplicon.index.to_list()
        # Get depths and variant data
        depth_data = get_depth_data(df_samples_amplicon, len(ref_seq), amplicon)
        variant_data = get_variant_data(df_samples_amplicon)
        amplicon_data = get_depth_amplicon(df_samples_amplicon)
        # Get Coverage Statistics
        coverage_stat = get_coverage_stat(depth_data, len(ref_seq), 10)
        # Get Amplicon feature to plot them with DNA gene feature
        bedfile = df_samples_amplicon['amplicon_region_file'][0].strip()
        amplicon_feature = get_region_amplicon(bedfile)
        gene_feature_len = len(gene_feature)
        for amplicon_name, region in amplicon_feature.items():
            index = int(amplicon_name.split('_')[-1])
            if index % 2:
                gene_feature.append(
                    dict(name=amplicon_name,
                         value={
                             "idx": gene_feature_len + index - 1,
                             "start": region[0],
                             "end": region[1],
                             "level": 95,
                             "strand": 1,
                             "type": "amplicon_feature"
                         },
                         itemStyle={"color": 'violet'})
                )
            else:
                gene_feature.append(
                    dict(name=amplicon_name,
                         value={
                             "idx": gene_feature_len + index - 1,
                             "start": region[0],
                             "end": region[1],
                             "level": 80,
                             "strand": 1,
                             "type": "amplicon_feature"
                         },
                         itemStyle={"color": 'skyblue'})
                )
    else:

        df_samples = pd.read_table(samples_data, names=['coverage_depth_file', 'vcf_file'], index_col=0, header=None)
        df_samples = df_samples.fillna(0)
        # Get list of samples name
        samples_name = df_samples.index.to_list()
        # Get depths and variant data
        depth_data = get_depth_data(df_samples, len(ref_seq), amplicon)
        variant_data = get_variant_data(df_samples)
        amplicon_data = {}
        # Get Coverage Statistics
        coverage_stat = get_coverage_stat(depth_data, len(ref_seq), 10)
        '''
        samples_name ={}
        depth_data={}
        variant_data={}
        coverage_stat=[]
        amplicon_data = {}
        '''

    # Parse README to HTML save to them About Tab
    dirpath = Path(__file__).parent
    readme = dirpath / 'readme/README.md'
    with open(readme, "r", encoding="utf-8") as input_file:
        text = input_file.read()
    about_html = markdown.markdown(text, extensions=['tables', 'nl2br', 'extra'])

    # Save Coverage plot to HTML file
    write_html_coverage_plot(samples_name=samples_name,
                             depth_data=depth_data,
                             variant_data=variant_data,
                             ref_seq=ref_seq,
                             coverage_stat=coverage_stat,
                             gene_feature=gene_feature,
                             about_html=about_html,
                             output_html=output_html,
                             amplicon_data=amplicon_data,
                             amplicon=amplicon)
    logging.info(f'Wrote HTML Coverage Plot to "{output_html}"')


if __name__ == "__main__":
    app()