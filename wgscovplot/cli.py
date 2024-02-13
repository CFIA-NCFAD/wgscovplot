"""Console script."""

import logging
from pathlib import Path
from sys import version_info
from typing import Optional

import typer
from rich.logging import RichHandler

from wgscovplot.__about__ import __version__
from wgscovplot.wgscovplot import run

app = typer.Typer()
logger = logging.getLogger(__name__)
VERSION = f"wgscovplot version {__version__}; Python {version_info.major}.{version_info.minor}.{version_info.micro}"


def version_callback(value: bool):
    if value:
        typer.echo(f"wgscovplot version {__version__}")
        raise typer.Exit()


def check_dir_exists_callback(path: Path) -> Path:
    if not (path.exists() or path.is_dir()):
        msg = (
            f"An existing Nextflow workflow results directory must be specified! '{path}' does not exist or is not a"
            " directory!"
        )
        raise typer.BadParameter(msg)
    return path


@app.command(epilog=VERSION)
def main(
    input_dir: Path = typer.Argument(
        ...,
        callback=check_dir_exists_callback,
        help="Directory containing Mosdepth and variant calling results from "
        "sequence analysis. For example, the output directory from execution of "
        "the nf-core/viralrecon or CFIA-NCFAD/nf-flu Nextflow workflow",
    ),
    output_html: Path = typer.Option("wgscovplot.html", "-o", "--output-html", help="wgscovplot HTML output file"),
    primers_fasta: Path = typer.Option(
        None, "-p", "--primers-fasta", help="FASTA file containing real-time PCR primer/probe sequences."
    ),
    low_coverage_threshold: int = typer.Option(
        10, "-l", "--low-coverage-threshold", help="Low sequencing coverage threshold."
    ),
    edit_distance: int = typer.Option(
        0,
        "-d",
        "--edit-distance",
        help="The maximum differences or 'edits' allowed between real-time "
        "PCR primer/probe sequences and the sample sequences.",
    ),
    compress_depths: bool = typer.Option(default=True, is_flag=True, help="Compress coverage depth arrays?"),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Verbose logs"),
    force: bool = typer.Option(
        False, "-f", "--force", is_flag=True, show_default=False, help="Force overwrite of existing output files"
    ),
    version: Optional[bool] = typer.Option(  # noqa: ARG001
        None, callback=version_callback, help=f'Print {"wgscovplot version"} and exit'
    ),
):
    init_logging(verbose)
    logger.info(VERSION)
    if output_html.exists() and not force:
        msg = f"Output file '{output_html}' exists! Use --force to overwrite."
        logger.error(msg)
        raise typer.Exit(1)
    logger.info(f"{input_dir=}")
    logger.info(f"{output_html=}")
    run(input_dir, low_coverage_threshold, primers_fasta, edit_distance, output_html, compress_depths)


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
