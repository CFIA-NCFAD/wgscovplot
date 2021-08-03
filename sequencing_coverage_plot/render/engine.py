from jinja2 import Environment
from typing import Optional
from ..globals_parameters import Resources


def write_utf8_html_file(file_name: str, html_content: str):
    with open(file_name, "w+", encoding="utf-8") as html_file:
        html_file.write(html_content)


class RenderEngine:
    def __init__(self, env: Optional[Environment] = None):
        self.env = env
