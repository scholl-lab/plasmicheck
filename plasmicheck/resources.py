"""Resource path resolution for the plasmicheck package."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path


def get_resource_path(resource: str) -> Path:
    """Resolve a resource path within the plasmicheck package.

    Uses importlib.resources for reliable resolution whether running
    from source or from an installed package.

    Parameters
    ----------
    resource : str
        Path relative to the plasmicheck package root
        (e.g., "templates/report_template.html").
    """
    return Path(str(files("plasmicheck").joinpath(resource)))
