import logging
from pathlib import Path
from subprocess import run
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class CustomBuildHook(BuildHookInterface):
    def initialize(self, version: str, build_data: dict[str, Any]) -> None:  # noqa: ARG002
        logger.warning("CustomBuildHook: running 'npm run build' to compile wgscovplot.js")
        workdir = Path(self.root, "web")
        run(
            "npm run build",  # noqa: S607
            cwd=str(workdir),
            check=True,
            shell=True,  # noqa: S602
        )
        build_data["artifacts"].append("web/build/wgscovplot.js")
        logger.warning("Done! 'wgscovplot.js' built!")
