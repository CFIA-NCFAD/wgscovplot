#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
from os.path import exists
from typer.testing import CliRunner

from sequencing_coverage_plot.coverage_plot.cli import app

runner = CliRunner()

dirpath = Path(__file__).parent
input_ref = dirpath/'data/nCoV-2019.reference.fasta'
input_genbank = dirpath/'data/sequence_sars_cov2.gb'


def test_cli():
    assert input_ref.exists()
    assert input_genbank.exists()
    result = runner.invoke(app)
    result = runner.invoke(app)
    assert result.exit_code != 0
    assert 'Error: Missing option' in result.output
    help_result = runner.invoke(app, ['--help'])
    assert help_result.exit_code == 0
    assert 'Show this message and exit.' in help_result.output
