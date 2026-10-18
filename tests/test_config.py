"""Unit tests for plasmicheck.config."""

from __future__ import annotations

import pytest

from plasmicheck.config import _reset, get_config


class TestGetConfig:
    @pytest.mark.unit
    def test_returns_dict(self) -> None:
        cfg = get_config()
        assert isinstance(cfg, dict)

    @pytest.mark.unit
    def test_singleton_identity(self) -> None:
        _reset()
        a = get_config()
        b = get_config()
        assert a is b

    @pytest.mark.unit
    def test_contains_expected_keys(self) -> None:
        cfg = get_config()
        assert "default_threshold" in cfg
        assert "scoring" in cfg
        assert "alignment" in cfg

    @pytest.mark.unit
    def test_reset_clears_cache(self) -> None:
        a = get_config()
        _reset()
        b = get_config()
        assert a is not b
        assert a == b  # same content

    @pytest.mark.unit
    def test_removed_dead_keys_and_retained_figsize(self) -> None:
        cfg = get_config()
        # plot_sample_report.figsize must be retained for matplotlib_backend
        assert "figsize" in cfg["plot_sample_report"]
        assert cfg["plot_sample_report"]["figsize"] == {"width": 1000, "height": 400}
        # Dead keys in plot_summary and plot_sample_report must be removed
        assert "colors" not in cfg["plot_sample_report"]
        assert "title_box_plot" not in cfg["plot_sample_report"]
        assert "figsize" not in cfg["plot_summary"]
        assert "colors" not in cfg["plot_summary"]
        assert "plot_dimensions" not in cfg["plot_summary"]
