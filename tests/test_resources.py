"""Unit tests for plasmicheck.resources."""

from __future__ import annotations

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
        assert not path.exists()
