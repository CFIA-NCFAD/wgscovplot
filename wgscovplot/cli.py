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

app = typer.Typer()


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
    epilog=f"wgscovplot version {__version__}; Python {version_info.major}.{version_info.minor}.{version_info.micro}"
)
def main(
        input_dir: Path = typer.Option(..., callback=check_dir_exists_callback, help="Nextflow workflow results "
                                                                                     "directory"),
        output_html: Path = typer.Option("wgscovplot.html", help="Output File of Interactive HTML Coverage Plot"),
        ref_seq: Path = typer.Option(None, help="Path to reference sequences"),
        genbank: Path = typer.Option(None, help="Genbank file contains gene features"),
        ncbi_accession_id: str = typer.Option(default="", help="NCBI accession id to fetch gene features "
                                                               "and/or reference sequences"),
        amplicon: bool = typer.Option(default=False, help="Plot Amplicon Coverage Depth"),
        gene_feature: bool = typer.Option(default=False, help="Plot Gene Features"),
        gene_misc_feature: bool = typer.Option(default=False, help="Plot Miscellaneous Features"),
        verbose: bool = typer.Option(default=False, help="Verbose logs"),
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
    run(input_dir, ref_seq, genbank, ncbi_accession_id, amplicon, gene_feature, gene_misc_feature, output_html)


if __name__ == "__main__":
    app()
