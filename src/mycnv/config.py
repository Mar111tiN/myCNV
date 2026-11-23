"""Parameter handling: packaged defaults (the somVarWES study config) deep-merged with an optional YAML file."""

from __future__ import annotations

import copy
from importlib import resources

import yaml

from .tools import log


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


# subtrees whose keys are free (column prefixes / aggregations of the rolling step)
_FREE = {("rolling", "cov", "data"), ("rolling", "snp", "data")}
_ALIASES = {"diffexp"}  # legacy spelling of rolling.*.diff_exp


def unknown_keys(new: dict, defaults: dict, path: tuple = ()) -> list[str]:
    """Dotted paths of keys in `new` that the defaults (= the keys mycnv reads) do not have."""
    out = []
    for key, value in new.items():
        sub = path + (key,)
        if key not in defaults:
            if key not in _ALIASES:
                out.append(".".join(sub))
        elif isinstance(value, dict) and isinstance(defaults[key], dict) and sub not in _FREE:
            out += unknown_keys(value, defaults[key], sub)
    return out


def _known(new: dict, defaults: dict, path: tuple = ()) -> dict:
    """`new` without the keys reported by unknown_keys."""
    out = {}
    for key, value in new.items():
        sub = path + (key,)
        if key not in defaults:
            if key in _ALIASES:
                out[key] = value
        elif isinstance(value, dict) and isinstance(defaults[key], dict) and sub not in _FREE:
            out[key] = _known(value, defaults[key], sub)
        else:
            out[key] = value
    return out


def load_config(path: str | None = None, overrides: dict | None = None) -> dict:
    """Defaults, then the YAML file at `path` (top-level `CNV:` optional), then `overrides`.

    Keys that mycnv does not read (e.g. the legacy `pileup`, `ascat`, `threads` entries, or a misspelt key) are
    reported with a warning and ignored. A bad `sex:` section raises ValueError here, before any command does work.
    """
    from .sex import thresholds

    cfg = default_config()
    if path:
        with open(path) as fh:
            user = _section(yaml.safe_load(fh) or {})
        extra = unknown_keys(user, cfg)
        if extra:
            log.warning("Config %s: ignoring keys mycnv does not use: %s", path, ", ".join(extra))
        _deep_update(cfg, _known(user, cfg))
    if overrides:
        _deep_update(cfg, copy.deepcopy(overrides))
    thresholds(cfg.get("sex"))
    return normalize_keys(cfg)


def normalize_keys(cfg: dict) -> dict:
    """Accept the legacy key `diffexp` as an alias of `diff_exp` in the rolling sections."""
    for section in cfg.get("rolling", {}).values():
        if isinstance(section, dict) and "diffexp" in section:
            section.setdefault("diff_exp", section.pop("diffexp"))
    return cfg


def dump_config(cfg: dict) -> str:
    return yaml.safe_dump({"CNV": cfg}, sort_keys=False)
