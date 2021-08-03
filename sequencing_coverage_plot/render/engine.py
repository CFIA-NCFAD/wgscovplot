from jinja2 import Environment
from typing import Optional


class RenderEngine:
    def __init__(self, env: Optional[Environment] = None):
        self.env = env
