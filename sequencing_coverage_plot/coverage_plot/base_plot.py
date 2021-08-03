from typing import Any, Optional, Union
import uuid
import simplejson as json
from pathlib import Path
from jinja2 import Environment
from ..ultilities import utils
from ..render import engine
from ..globals_parameters import HTMLPageConfig, RenderType, ThemeType


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
        _opts = init_opts
        if isinstance(init_opts, InitOpts):
            _opts = init_opts.opts

        self.width = _opts.get("width", "900px")
        self.height = _opts.get("height", "500px")
        self.renderer = _opts.get("renderer", RenderType.CANVAS)
        self.page_title = _opts.get("page_title", HTMLPageConfig.PAGE_TITLE)
        self.theme = _opts.get("theme", ThemeType.WHITE)
        self.chart_id = _opts.get("chart_id") or uuid.uuid4().hex
        self.options: dict = {}

    def get_options(self) -> dict:
        return utils.remove_key_with_none_value(self.options)

    def render_html(
            self,
            output_html: Path = "coverage_plot.html",
            template_file: str = "coverage_bar_chart.html",
            env: Optional[Environment] = None,
            **kwargs,
    ) -> Path:
        return self
