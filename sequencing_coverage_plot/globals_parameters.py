from enum import Enum


class Resources(Enum):
    ECHARTS: str = "https://cdn.jsdelivr.net/npm/echarts@5.1.2/dist/echarts.min.js"


class RenderType(Enum):
    CANVAS: str = "canvas"
    SVG: str = "svg"


class FileType(Enum):
    SVG: str = "svg"
    PNG: str = "png"
    PDF: str = "pdf"
    JPEG: str = "jpeg"
    HTML: str = "html"


class ChartType(Enum):
    BAR: str = "bar"


class ThemeType(Enum):
    LIGHT = "light"
    DARK = "dark"
    WHITE = "white"
