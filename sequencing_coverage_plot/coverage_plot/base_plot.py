from typing import Any, Optional, Union
import uuid
import simplejson as json
from pathlib import Path
from jinja2 import Environment
from ..ultilities import utils
from ..render import engine
from ..globals_parameters import HTMLPageConfig, RenderType, ThemeType
import logging
from rich.logging import RichHandler


class BasicOpts:
    __slots__ = ("opts",)

    def update(self, **kwargs):
        self.opts.update(kwargs)

    def get(self, key: str) -> Any:
        return self.opts.get(key)


class InitOpts(BasicOpts):
    def __init__(
        self,
        width: str = "900px",
        height: str = "500px",
        chart_id: Optional[str] = None,
        renderer: str = RenderType.CANVAS,
        page_title: str = HTMLPageConfig.PAGE_TITLE,
        theme: str = ThemeType.WHITE,
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

    def __init__(self, init_opts: Union[InitOpts, dict] = InitOpts()):
        from rich.traceback import install
        install(show_locals=True)
        logging.basicConfig(
            format="%(message)s",
            datefmt="[%Y-%m-%d %X]",
            level=logging.INFO,
            handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
        )
        _opts = init_opts
        if isinstance(init_opts, InitOpts):
            _opts = init_opts.opts
        self.width = _opts.get("width", "900px")
        self.height = _opts.get("height", "500px")
        self.renderer = _opts.get("renderer", RenderType.CANVAS)
        self.page_title = _opts.get("page_title", HTMLPageConfig.PAGE_TITLE)
        self.theme = _opts.get("theme", ThemeType.WHITE)
        self.chart_id = _opts.get("chart_id") or uuid.uuid4().hex
        self.options: dict = {}  # store all options and data for Echarts

    def get_options(self) -> dict:
        return utils.remove_key_with_none_value(self.options)

    def dump_options(self) -> str:
        return utils.replace_placeholder(
            json.dumps(self.get_options(), indent=4, ignore_nan=True)
        )

    def render_html(
        self,
        path: Path = "render.html",
        template_file: Path = "coverage_bar_chart.html",
        env: Optional[Environment] = None,
        **kwargs,
    ) -> Path:
        logging.info(f'Jinja2 template file: "{template_file}"')
        self._prepare_render()
        return engine.render_html(self, path, template_file, env, **kwargs)

    def _prepare_render(self):
        self.json_contents = self.dump_options()
