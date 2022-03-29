#!/usr/bin/env python
# -*- coding: utf-8 -*-
from os.path import exists
from pathlib import Path
from typer.testing import CliRunner
from wgscovplot.cli import app

runner = CliRunner()

dirpath = Path(__file__).parent
input_ref = dirpath / 'data/nCoV-2019.reference.fasta'
input_genbank = dirpath / 'data/sequence_sars_cov2.gb'
input_dir = dirpath / 'data/tools'


def test_cli():
    assert input_ref.exists()
    assert input_genbank.exists()
    result = runner.invoke(app)
    assert result.exit_code != 0
    assert 'Error: Missing option' in result.output
    help_result = runner.invoke(app, ['--help'])
    assert help_result.exit_code == 0
    assert 'Show this message and exit.' in help_result.output

    with runner.isolated_filesystem():
        out_html = 'wgscovplot_test1.html'
        test_result = runner.invoke(app, ['--input-dir', str(input_dir.resolve().absolute()),
                                          '--ref-seq', str(input_ref.absolute()),
                                          '--output-html', out_html,
                                          '--genbank', str(input_genbank.absolute()),
                                          '--gene-feature',
                                          '--no-amplicon'])
        assert test_result.exit_code == 0
        assert exists(out_html)

    with runner.isolated_filesystem():
        out_html = 'wgscovplot_test2.html'
        test_result = runner.invoke(app, ['--input-dir', str(input_dir.resolve().absolute()),
                                          '--ref-seq', str(input_ref.absolute()),
                                          '--output-html', out_html,
                                          '--genbank', str(input_genbank.absolute()),
                                          '--gene-feature'])
        assert test_result.exit_code == 0
        assert exists(out_html)

    with runner.isolated_filesystem():
        out_html = 'wgscovplot_test3.html'
        test_result = runner.invoke(app, ['--input-dir', str(input_dir.resolve().absolute()),
                                          '--ref-seq', str(input_ref.absolute()),
                                          '--output-html', out_html,
                                          '--gene-feature',
                                          '--no-amplicon'])
        assert test_result.exit_code == 1

    with runner.isolated_filesystem():
        out_html = 'wgscovplot_test4.html'
        ncbi_accession_id = 'MN908947.3'
        test_result = runner.invoke(app, ['--input-dir', str(input_dir.resolve().absolute()),
                                          '--ncbi-accession-id', str(ncbi_accession_id),
                                          '--output-html', out_html,
                                          '--gene-feature',
                                          '--no-amplicon'])
        assert test_result.exit_code == 0
        assert exists(out_html)

    with runner.isolated_filesystem():
        test_result = runner.invoke(app, ['--input-dir', str(input_dir.resolve().absolute())])
        assert test_result.exit_code == 0
        assert exists('wgscovplot.html')
