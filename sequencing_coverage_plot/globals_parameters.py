from pathlib import Path
from jinja2 import Environment, FileSystemLoader


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
    HTML_ENV = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader()
    )


RenderType = _RenderType()
FileType = _FileType()
ThemeType = _ThemeType()
HTMLPageConfig = _HTMLPageConfig()
