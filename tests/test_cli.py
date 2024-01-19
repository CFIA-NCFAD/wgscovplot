#!/usr/bin/env python
from os.path import exists
from pathlib import Path

from typer.testing import CliRunner

from wgscovplot.cli import app

runner = CliRunner()

dirpath = Path(__file__).parent
non_segmented_virus_data_dir = dirpath / "data/sars-cov-2"
segmented_virus_data_dir = dirpath / "data/iav"


def test_basic_cli():
    result = runner.invoke(app)
    assert result.exit_code != 0
    assert "Error" in result.output
    assert "Missing argument" in result.output
    help_result = runner.invoke(app, ["--help"])
    assert help_result.exit_code == 0
    assert "wgscovplot version" in help_result.output


def test_cli_non_segment_virus():
    with runner.isolated_filesystem():
        out_html = "wgscovplot_test1.html"
        test_result = runner.invoke(
            app,
            [
                str(non_segmented_virus_data_dir.resolve().absolute()),
                "--output-html",
                out_html,
            ],
        )
        assert test_result.exit_code == 0
        assert exists(out_html)


def test_cli_segment_virus():
    with runner.isolated_filesystem():
        out_html = "wgscovplot_test_segment_virus.html"
        test_result = runner.invoke(
            app, [str(segmented_virus_data_dir.resolve().absolute()), "--output-html", out_html]
        )
        assert test_result.exit_code == 0
        assert exists(out_html)
