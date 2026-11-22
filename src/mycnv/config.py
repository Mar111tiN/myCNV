"""Parameter handling: packaged defaults (the somVarWES study config) deep-merged with an optional YAML file."""

from __future__ import annotations

import copy
from importlib import resources

import yaml


def _deep_update(base: dict, new: dict) -> dict:
    for key, value in new.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value
    return base


def _section(cfg: dict) -> dict:
    return cfg.get("CNV", cfg) if isinstance(cfg, dict) else {}


def default_config() -> dict:
    text = resources.files("mycnv").joinpath("data", "default_config.yaml").read_text()
    return _section(yaml.safe_load(text))


def load_config(path: str | None = None, overrides: dict | None = None) -> dict:
    """Defaults, then the YAML file at `path` (top-level `CNV:` optional), then `overrides`."""
    cfg = default_config()
    if path:
        with open(path) as fh:
            _deep_update(cfg, _section(yaml.safe_load(fh) or {}))
    if overrides:
        _deep_update(cfg, copy.deepcopy(overrides))
    return normalize_keys(cfg)


def normalize_keys(cfg: dict) -> dict:
    """Accept the legacy key `diffexp` as an alias of `diff_exp` in the rolling sections."""
    for section in cfg.get("rolling", {}).values():
        if isinstance(section, dict) and "diffexp" in section:
            section.setdefault("diff_exp", section.pop("diffexp"))
    return cfg


def dump_config(cfg: dict) -> str:
    return yaml.safe_dump({"CNV": cfg}, sort_keys=False)
