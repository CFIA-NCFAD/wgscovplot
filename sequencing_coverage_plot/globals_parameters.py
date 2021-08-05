from pathlib import Path
from jinja2 import Environment, FileSystemLoader


class _Resources:
    ECHARTS: str = "https://cdn.jsdelivr.net/npm/echarts@5.1.2/dist/echarts.min.js"


class _RenderType:
    CANVAS: str = "canvas"
    SVG: str = "svg"


class _FileType:
    SVG: str = "svg"
    PNG: str = "png"
    PDF: str = "pdf"
    JPEG: str = "jpeg"
    HTML: str = "html"


class _ChartType:
    BAR: str = "bar"


class _ThemeType:
    LIGHT = "light"
    DARK = "dark"
    WHITE = "white"


class _HTMLPageConfig:
    PAGE_TITLE = "Interactive Coverage Plot"
    GLOBAL_ENV = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "render", "templates")),
    )


RenderType = _RenderType()
FileType = _FileType()
ThemeType = _ThemeType()
HTMLPageConfig = _HTMLPageConfig()
Resources = _Resources()
