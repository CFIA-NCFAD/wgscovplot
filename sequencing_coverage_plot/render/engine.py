from jinja2 import Environment
from typing import Optional, Any
from pathlib import Path
from ..globals_parameters import Resources, HTMLPageConfig
from ..ultilities import utils
import logging


def write_utf8_html_file(file_name: Path, html_content: str):
    logging.info(f'Rendering to output file: "{file_name}"')
    with open(file_name, "w+", encoding="utf-8") as html_file:
        html_file.write(html_content)


class RenderEngine:

    def __init__(self, env: Optional[Environment] = None):
        self.env = env or HTMLPageConfig.GLOBAL_ENV

    @staticmethod
    def get_js_link(chart: Any) -> Any:
        chart.dependencies = [Resources.ECHARTS]
        return chart

    def render_chart_to_file(self, template_file: str, chart: Any, path: Path, **kwargs):
        """
        Render a chart or page to local html files.

        :param chart: A Chart or Page object
        :param path: The destination file which the html code write to
        :param template_file: The name of template file.
        """
        tpl = self.env.get_template(template_file)
        html = utils.replace_placeholder(
            tpl.render(chart=self.get_js_link(chart), **kwargs)
        )
        write_utf8_html_file(path, html)


def render_html(
        chart,
        path: Path,
        template_file: Path,
        env: Optional[Environment] = None,
        **kwargs
) -> Path:
    RenderEngine(env).render_chart_to_file(template_file=template_file, chart=chart, path=path, **kwargs)
    return path
