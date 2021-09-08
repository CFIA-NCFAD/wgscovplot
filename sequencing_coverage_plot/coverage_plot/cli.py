# -*- coding: utf-8 -*-

"""Console script."""
import logging
from pathlib import Path
from typing import Optional
import typer
from rich.logging import RichHandler
from sys import version_info
from Bio.SeqIO.FastaIO import SimpleFastaParser
from .coverage_plot import (write_html_coverage_plot, prepare_data, get_gene_feature)
from coverage_plot import __version__

app = typer.Typer()


def version_callback(value: bool):
    if value:
        typer.echo(f"shicp version {__version__}")
        raise typer.Exit()


@app.command(
    epilog=f"shicp version {__version__}; Python {version_info.major}.{version_info.minor}.{version_info.micro}"
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
        bed: Optional[Path] = typer.Option(None,
                                           "-b",
                                           "--bed",
                                           help="Bed file"),
        genbank: Path = typer.Option(...,
                                     "-g",
                                     "--genbank",
                                     help="Genbank file contains features of reference sequence"),
        verbose: bool = typer.Option(False, help="Verbose logs"),
        version: Optional[bool] = typer.Option(
            None,
            callback=version_callback,
            help=f'Print "shicp version" and exit',
        )
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
    samples_name, depth_data, variant_data, coverage_stat = prepare_data(samples_data)
    write_html_coverage_plot(samples_name=samples_name, depth_data=depth_data, variant_data=variant_data,
                             ref_seq=ref_seq, coverage_stat=coverage_stat, gene_feature=gene_feature,
                             output_html=output_html)
    logging.info(f'Wrote HTML Coverage Plot to "{output_html}"')


if __name__ == "__main__":
    app()
