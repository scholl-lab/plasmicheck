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
