# -*- coding: utf-8 -*-

"""Console script."""
import logging
import typer
import pandas as pd
from pathlib import Path
from typing import Optional
from rich.logging import RichHandler
from sys import version_info
from Bio.SeqIO.FastaIO import SimpleFastaParser
from .wgscovplot import (write_html_coverage_plot, get_samples_name,
                         get_variant_data, get_depth_data, get_gene_feature, get_amplicon_regions,
                         get_amplicon_depths, get_amplicon_feature)
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
    if ref_seq:
        with open(ref_seq) as fh:
            for name, seq in SimpleFastaParser(fh):
                ref_seq = seq

    gene_feature = get_gene_feature(genbank)

    dirpath = Path(__file__).parent
    readme = dirpath / 'readme/README.md'
    with open(readme, "r", encoding="utf-8") as input_file:
        text = input_file.read()
    about_html = markdown.markdown(text, extensions=['tables', 'nl2br', 'extra'])

    if amplicon:
        df_samples_amplicon = pd.read_table(samples_data, names=['amplicon_region_file', 'amplicon_perbase_file', 'vcf_file'],
                                            index_col=0, header=None)
        df_samples_amplicon = df_samples_amplicon.fillna(0)
        samples_name = get_samples_name(df_samples_amplicon)
        depth_data, coverage_stat = get_amplicon_depths(len(ref_seq), df_samples_amplicon)
        variant_data = get_variant_data(df_samples_amplicon)
        # amplicon_feature = get_amplicon_feature (df_samples_amplicon)
        
        bedfile = df_samples_amplicon['amplicon_region_file'][0].strip()
        amplicon_feature = get_amplicon_regions(bedfile)
        gene_feature_len = len(gene_feature)

        for amplicon_name, region in amplicon_feature.items():
            index = int(amplicon_name.split('_')[-1])
            if index % 2:
                gene_feature.append(
                    dict(name=amplicon_name,
                         value=[gene_feature_len + index - 1, region[0], region[1], 90, 1],
                         itemStyle={"color": 'violet'})
                )
            else:
                gene_feature.append(
                    dict(name=amplicon_name,
                         value=[gene_feature_len + index - 1, region[0], region[1], 80, 1],
                         itemStyle={"color": 'skyblue'})
                )
    else:
        df_samples = pd.read_table(samples_data, names=['coverage_depth_file', 'vcf_file'], index_col=0, header=None)
        df_samples = df_samples.fillna(0)
        samples_name = get_samples_name(df_samples)
        depth_data, coverage_stat = get_depth_data(df_samples, 10)
        variant_data = get_variant_data(df_samples)
    write_html_coverage_plot(amplicon=amplicon, samples_name=samples_name, depth_data=depth_data, variant_data=variant_data,
                             ref_seq=ref_seq, coverage_stat=coverage_stat, gene_feature=gene_feature,
                             about_html=about_html, output_html=output_html)
    logging.info(f'Wrote HTML Coverage Plot to "{output_html}"')


if __name__ == "__main__":
    app()
