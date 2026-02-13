"""Centralized configuration for plasmicheck.

Loads config.json once at first access. All modules import from here
instead of loading config.json independently.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_config: dict[str, Any] | None = None
_CONFIG_PATH = Path(__file__).resolve().parent / "config.json"


def get_config() -> dict[str, Any]:
    """Return the loaded configuration, loading from disk on first call."""
    global _config
    if _config is None:
        with open(_CONFIG_PATH) as f:
            _config = json.load(f)
    return _config


def _reset() -> None:
    """Reset cached config (for testing only)."""
    global _config
    _config = None
