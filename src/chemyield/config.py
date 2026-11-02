"""YAML config loading with dotted-key command-line overrides."""

from __future__ import annotations

import copy
from pathlib import Path

import yaml


def deep_update(base: dict, upd: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in upd.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def parse_overrides(items) -> dict:
    """['split.kind=scaffold', 'seeds=[0,1]'] -> nested dict (values parsed as YAML)."""
    out: dict = {}
    for item in items or []:
        key, _, raw = item.partition("=")
        if not _:
            raise ValueError(f"Override must look like a.b=value, got '{item}'")
        cur = out
        parts = key.split(".")
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = yaml.safe_load(raw)
    return out


def load_config(path: str | Path, overrides=()) -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_config_path"] = str(path)
    return deep_update(cfg, parse_overrides(overrides))
