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
        has_npm = has_software("npm")
        has_bun = has_software("bun")
        js_installer = "npm" if has_npm else "bun" if has_bun else None
        if js_installer is None:
            msg = (
                "No JS installer found. Skipping wgscovplot.js compile and build. "
                "Please install 'npm' or 'bun' and try again."
            )
            raise RuntimeError(msg)
        run(
            f"{js_installer} install",
            cwd=str(workdir),
            check=True,
            shell=True,  # noqa: S602
        )
        run(
            f"{js_installer} run build",
            cwd=str(workdir),
            check=True,
            shell=True,  # noqa: S602
        )
        build_data["artifacts"].append("web/build/wgscovplot.js")
        logger.warning("Done! 'wgscovplot.js' built!")


def has_software(name: str) -> bool:
    return run(["which", name], capture_output=True, check=False).returncode == 0  # noqa: S603,S607
