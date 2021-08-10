from typing import Any
from pathlib import Path
import logging
from jinja2 import Environment, FileSystemLoader
from ..globals_parameters import Resources
from ..ultilities import utils


render_env = Environment(
    keep_trailing_newline=True,
    trim_blocks=True,
    lstrip_blocks=True,
    loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "templates")),
)


def write_utf8_html_file(file_name: Path, html_content: str):
    logging.info(f'Rendering to output file: "{file_name}"')
    with open(file_name, "w+", encoding="utf-8") as html_file:
        html_file.write(html_content)


def get_js_link(chart: Any) -> Any:
    chart.dependencies = [Resources.ECHARTS.value]
    return chart


def render_chart_to_file(template_file: str, chart: Any, path: Path, **kwargs):
    """
    Render a chart or page to local html files.

    :param chart: A Chart or Page object
    :param path: The destination file which the html code write to
    :param template_file: The name of template file.
    """
    tpl = render_env.get_template(template_file)
    html = tpl.render(chart=get_js_link(chart), **kwargs)
    write_utf8_html_file(path, html)


def render_html(
        chart,
        path: Path,
        template_file: Path,
        **kwargs
) -> Path:
    render_chart_to_file(template_file=template_file, chart=chart, path=path, **kwargs)
    return path
