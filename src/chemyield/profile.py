"""Profile a candidate dataset BEFORE committing to it.

    python -m chemyield.profile --config configs/<candidate>.yaml [--sample 20000]
                                [--group-by scaffold role:ligand ...] [--out results/<name>_profile.md]

Runs the dataset's adapter + cleaning (same code as a real run), then reports the things that decide
whether the data is usable: how many reactions survive cleaning, what the yield distribution looks
like (0 % / near-100 % pile-ups), how many distinct components each role has, how many reactions
share the same product, and whether group / group-val splits can actually be built. --sample takes a
random subset of the RAW rows first (cleaning canonicalises every SMILES, so big files are slow).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .adapters import build_adapter
from .cleaning import clean_reactions
from .config import load_config
from .splits import _groups, make_split


def _group_row(df: pd.DataFrame, spec: str) -> str:
    groups = _groups(df, spec)
    sizes = np.array(sorted((len(v) for v in groups.values()), reverse=True))
    n = len(df)
    can_group = can_gval = "yes"
    try:
        make_split(df, "group", seed=0, group_by=spec)
    except ValueError as e:
        can_group = "no (" + str(e).split("(")[0].strip()[:40] + ")" if "too few" not in str(e) else "no"
    try:
        make_split(df, "group", seed=0, group_by=spec, val_by_group=True)
    except ValueError:
        can_gval = "no"
    return (f"| {spec} | {len(groups)} | {sizes[0] / n:.1%} | {int((sizes >= 20).sum())} | "
            f"{np.median(sizes):.0f} | {can_group} | {can_gval} |")


def profile(df: pd.DataFrame, report: dict | None = None, group_bys=("scaffold", "product"),
            name: str = "dataset") -> str:
    y = df["yield"].to_numpy(float)
    n = len(df)
    L = [f"# Profile: {name}", ""]
    if report:
        L += [f"Raw rows **{report['n_raw']}** -> clean rows **{report['n_clean']}**", ""]
        L += ["Cleaning steps (what was removed and why):", ""]
        L += [f"- {s['step']}: removed {s['removed']}, remaining {s['remaining']}"
              + (f" ({s['note']})" if s.get("note") else "") for s in report.get("steps", [])] + [""]
    L += ["## Yield", "",
          f"n = {n}; mean {y.mean():.1f}, std {y.std():.1f}, median {np.median(y):.1f}",
          f"- exactly 0: {np.mean(y <= 0):.1%}   |   >= 90: {np.mean(y >= 90):.1%}   |   "
          f">= 99: {np.mean(y >= 99):.1%}",
          "- deciles: " + ", ".join(f"{q:.0f}" for q in np.percentile(y, range(0, 101, 10))),
          "- (a piled-up distribution means a constant predictor is already hard to beat, and "
          "the benchmark's `mean baseline` interval width will be small)", ""]
    sizes = df["components"].map(len)
    L += ["## Components", "",
          f"components per reaction: min {sizes.min()}, median {int(sizes.median())}, max {sizes.max()}"]
    roles = report.get("unique_components_by_role") if report else None
    if roles:
        L += [f"- unique {r}: {k}" for r, k in roles.items()]
    L += [f"- unique products: {df['product'].nunique()} (reactions per product: median "
          f"{df.groupby('product').size().median():.0f})", ""]
    if report and report.get("role_signatures"):
        L += ["Role signatures (most common):", ""]
        L += [f"- {s['n_reactions']} reactions: {s['roles']}" for s in report["role_signatures"][:6]] + [""]
    L += ["## Group structure (can OOD splits be built?)", "",
          "| group_by | n groups | largest group | groups with >=20 rows | median group size | group split | group-val split |",
          "|---|---|---|---|---|---|---|"]
    for spec in group_bys:
        try:
            L.append(_group_row(df, spec))
        except Exception as e:  # unknown role etc.
            L.append(f"| {spec} | error: {str(e)[:60]} | | | | | |")
    L += ["", "Read: a group split needs a few dozen groups, none dominating (largest group "
          "well under 20 % of rows); `group-val` needs more groups still.", ""]
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--sample", type=int, help="profile a random subset of this many RAW rows")
    ap.add_argument("--group-by", nargs="*", default=["scaffold", "product"])
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    cfg = load_config(args.config, [])
    raw = build_adapter(cfg["dataset"]).load()
    if args.sample and args.sample < len(raw):
        raw = raw.sample(args.sample, random_state=0).reset_index(drop=True)
    df, report = clean_reactions(raw, cfg.get("cleaning"))
    text = profile(df, report, args.group_by, name=cfg["dataset"]["name"])
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"Wrote {args.out}")
    else:
        print(text)


if __name__ == "__main__":
    main()
