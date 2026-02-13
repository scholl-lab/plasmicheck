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

    Returns
    -------
    pathlib.Path
        A real filesystem path to the requested resource.

    Raises
    ------
    RuntimeError
        If the resource is not available as a real filesystem path
        (for example, when the package is loaded from a zip/pyz
        archive or other non-filesystem importer).
    """
    resource_ref = files("plasmicheck").joinpath(resource)

    # Only return a concrete filesystem path when the resource is
    # actually backed by the local filesystem.  For non-filesystem
    # importers (e.g., zipimport/pyz), callers should use
    # importlib.resources.as_file() directly.
    if isinstance(resource_ref, Path):
        return resource_ref

    raise RuntimeError(
        f"Resource {resource!r} is not available as a real filesystem "
        "path. Use importlib.resources.as_file() on the resource "
        "Traversable if you require a temporary filesystem path."
    )
