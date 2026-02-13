"""Unit tests for plasmicheck.resources."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from plasmicheck.resources import get_resource_path


class TestGetResourcePath:
    @pytest.mark.unit
    def test_resolves_template(self) -> None:
        path = get_resource_path("templates/report_template.html")
        assert path.exists()
        assert path.name == "report_template.html"

    @pytest.mark.unit
    def test_resolves_logo(self) -> None:
        path = get_resource_path("static/img/plasmicheck_logo_200px.png")
        assert path.exists()

    @pytest.mark.unit
    def test_resolves_config(self) -> None:
        path = get_resource_path("config.json")
        assert path.exists()

    @pytest.mark.unit
    def test_nonexistent_returns_path(self) -> None:
        path = get_resource_path("does_not_exist.txt")
        assert isinstance(path, Path)
        assert not path.exists()

    @pytest.mark.unit
    def test_raises_for_non_filesystem_importer(self) -> None:
        """Raises RuntimeError when the resource is not a real filesystem path."""
        fake_traversable = MagicMock()
        fake_traversable.__class__ = type("ZipPath", (), {})
        with patch("plasmicheck.resources.files") as mock_files:
            mock_files.return_value.joinpath.return_value = fake_traversable
            with pytest.raises(RuntimeError, match="not available as a real filesystem path"):
                get_resource_path("some_resource.txt")
