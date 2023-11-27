# -*- coding: utf-8 -*-

"""Console script."""
import logging
import typer
from pathlib import Path
from typing import Optional
from rich.logging import RichHandler
from sys import version_info
from wgscovplot.wgscovplot import run
from wgscovplot import __version__
from textual.app import App, ComposeResult
from textual.widgets import Button,Input, Label, Static

app = typer.Typer()
logger = logging.getLogger(__name__)
VERSION = f"wgscovplot version {__version__}; Python {version_info.major}.{version_info.minor}.{version_info.micro}"

def version_callback(value: bool):
    if value:
        typer.echo(f"wgscovplot version {__version__}")
        raise typer.Exit()


def check_dir_exists_callback(path: Path) -> Path:
    if not (path.exists() or path.is_dir()):
        raise typer.BadParameter(f'An existing Nextflow workflow results directory must be specified! '
                                 f'"{path}" does not exist or is not a directory!')
    return path


@app.command(
    epilog=VERSION
)
def main(
        input_dir: Path = typer.Argument(...,
                                         callback=check_dir_exists_callback,
                                         help="Directory containing Mosdepth and variant calling results from "
                                              "sequence analysis. For example, the output directory from execution of "
                                              "the nf-core/viralrecon or CFIA-NCFAD/nf-flu Nextflow workflow"),
        output_html: Path = typer.Option("wgscovplot.html", help="wgscovplot HTML output file"),
        primers_fasta: Path = typer.Option(None, help="FASTA file containing real-time PCR primer/probe sequences."),
        low_coverage_threshold: int = typer.Option(default=10, help="Low sequencing coverage threshold."),
        show_amplicons: bool = typer.Option(default=True,
                                            help="Show amplicon positions and coverage along with sequencing coverage "
                                                 "plots."),
        show_gene_features: bool = typer.Option(default=True,
                                                help="Show genetic features such as CDS, genes, 5'/3'UTR along with "
                                                     "sequencing coverage plots."),
        is_segmented: bool = typer.Option(default=False,
                                          help="Output coverage plots for segmented viruses like Influenza A virus."),
        dev: bool = typer.Option(default=False, help="Run wgscovplot in development output mode."),
        interactive: bool = typer.Option(default=False, help="Run wgscovplot in interactive mode with textual."),
        max_primer_mismatches: int = typer.Option(default=0,
                                                  help="The maximum differences or 'edits' allowed between real-time "
                                                       "PCR primer/probe sequences and the sample sequences."),
        verbose: bool = typer.Option(default=False, help="Verbose logs"),
        version: Optional[bool] = typer.Option(None,
                                               callback=version_callback,
                                               help=f'Print {"wgscovplot version"} and exit')
):
    init_logging(verbose)
    logger.info(VERSION)
    logger.info(f'{input_dir=}')
    logger.info(f'{output_html=}')
    run(input_dir, low_coverage_threshold, is_segmented, primers_fasta,
        max_primer_mismatches, dev,
        output_html)

def init_logging(verbose):
    from rich.traceback import install
    install(show_locals=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.INFO if not verbose else logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


if __name__ == "__main__":
    app()
