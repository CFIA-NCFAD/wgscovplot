from typing import Any, Optional, Union
import uuid
import simplejson as json
from pathlib import Path
from ..ultilities import utils
from ..render import render
from ..globals_parameters import RenderType, ThemeType
import logging
from rich.logging import RichHandler
from ..options_dict import echarts_options


class InitOpts:
    def __init__(
            self,
            width: str = "900px",
            height: str = "500px",
            chart_id: Optional[str] = None,
            renderer: str = RenderType.CANVAS.value,
            page_title: str = "Sequencing Coverage Plot",
            theme: str = ThemeType.WHITE.value,
    ):
        self.opts: dict = {
            "width": width,
            "height": height,
            "chart_id": chart_id,
            "renderer": renderer,
            "page_title": page_title,
            "theme": theme,
        }


class BasePlot:

    def __init__(self, init_opts=InitOpts()):
        from rich.traceback import install
        install(show_locals=True)
        logging.basicConfig(
            format="%(message)s",
            datefmt="[%Y-%m-%d %X]",
            level=logging.INFO,
            handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
        )
        self.width = init_opts.opts.get("width", "900px")
        self.height = init_opts.opts.get("height", "500px")
        self.renderer = init_opts.opts.get("renderer", RenderType.CANVAS.value)
        self.page_title = init_opts.opts.get("page_title", "Interactive Coverage Plot")
        self.theme = init_opts.opts.get("theme", ThemeType.WHITE.value)
        self.chart_id = init_opts.opts.get("chart_id") or uuid.uuid4().hex
        self.options = echarts_options

    def get_options(self) -> dict:
        return utils.remove_key_with_none_value(self.options)

    def dump_options(self) -> str:
        return json.dumps(self.get_options(), indent=4, ignore_nan=True)

    def render_html(
        self,
        path: Path = "render.html",
        template_file: Path = "coverage_bar_chart.html",
        **kwargs,
    ) -> Path:
        logging.info(f'Jinja2 template file: "{template_file}"')
        self._prepare_render()
        return render.render_html(self, path, template_file, **kwargs)

    def _prepare_render(self):
        self.json_contents = self.dump_options()
