"""Dataset-agnostic reaction cleaning with an auditable report.

Steps (each is logged with how many rows it removed):
  1. drop rows with missing yield / components / product
  2. normalise yield scale (fraction -> percent), clip tiny overshoots, drop impossible values
  3. parse + canonicalise every SMILES (atom maps stripped); drop unparsable rows
  4. optional SMARTS role rules (e.g. split an unlabeled bag into aryl_halide / amine / ligand ...)
  5. element-balance sanity check (product elements must exist in the inputs)
  6. de-duplicate identical reactions (mean yield; flag disagreeing replicates)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from .chem import canonical_smiles, cached_mol, element_set

DEFAULTS = {
    "yield_scale": "auto",        # percent | fraction | auto
    "yield_range": [0.0, 100.0],
    "clip_tolerance": 5.0,        # values within this of the range are clipped, not dropped
    "drop_unparsable": True,
    "element_balance": "flag",    # drop | flag | off
    "dedupe": "mean",             # mean | median | first | off
    "conflict_std": 20.0,         # replicate std above this is a 'conflict'
    "drop_conflicts": False,
    "role_rules": [],             # [{role: ligand, smarts: "[P]", apply_to: [reactant, reagent]}]
}


def _apply_role_rules(role: str, smi: str, rules) -> str:
    mol = cached_mol(smi)
    for rule in rules:
        apply_to = rule.get("apply_to")
        if apply_to and role not in apply_to:
            continue
        patt = rule["_patt"]
        if patt is not None and mol is not None and mol.HasSubstructMatch(patt):
            return rule["role"]
    return role


def clean_reactions(df: pd.DataFrame, cfg: dict | None = None):
    cfg = {**DEFAULTS, **(cfg or {})}
    steps: list[dict] = []
    df = df.copy()

    def log(step, before, after, note=""):
        steps.append({"step": step, "removed": int(before - after), "remaining": int(after),
                      "note": note})

    n_raw = len(df)

    # 1. missing values -------------------------------------------------------
    df["yield"] = pd.to_numeric(df["yield"], errors="coerce")
    ok = (df["yield"].notna()
          & df["components"].map(lambda c: isinstance(c, list) and len(c) > 0)
          & df["product"].map(lambda p: isinstance(p, str) and p.strip() != ""))
    before = len(df)
    df = df[ok].reset_index(drop=True)
    log("missing yield/components/product", before, len(df))
    stats_before = df["yield"].describe().round(2).to_dict()

    # 2. yield scale ------------------------------------------------------------
    scale = cfg["yield_scale"]
    if scale == "auto":
        scale = "fraction" if df["yield"].quantile(0.99) <= 1.5 else "percent"
    if scale == "fraction":
        df["yield"] = df["yield"] * 100.0
    lo, hi = cfg["yield_range"]
    tol = cfg["clip_tolerance"]
    before = len(df)
    in_tol = df["yield"].between(lo - tol, hi + tol)
    n_clipped = int(((df["yield"] < lo) | (df["yield"] > hi))[in_tol].sum())
    df = df[in_tol].reset_index(drop=True)
    df["yield"] = df["yield"].clip(lo, hi)
    log("yield out of range", before, len(df),
        f"scale={scale}; {n_clipped} within ±{tol} clipped to [{lo},{hi}]")

    # 3. canonicalise -------------------------------------------------------------
    unsanitized: set[str] = set()
    new_comps, new_prod, keep = [], [], []
    for comps, prod in zip(df["components"], df["product"]):
        cc, good = [], True
        for role, smi in comps:
            can, san = canonical_smiles(smi)
            if can is None:
                good = False
                break
            if not san:
                unsanitized.add(can)
            cc.append((role, can))
        p_can, p_san = canonical_smiles(prod)
        if p_can is None:
            good = False
        elif not p_san:
            unsanitized.add(p_can)
        new_comps.append(cc)
        new_prod.append(p_can)
        keep.append(good)
    df["components"], df["product"] = new_comps, new_prod
    before = len(df)
    if cfg["drop_unparsable"]:
        df = df[np.array(keep)].reset_index(drop=True)
    log("unparsable SMILES", before, len(df),
        f"{len(unsanitized)} distinct molecules kept via unsanitized fallback (metal complexes etc.)")

    # 4. role rules -----------------------------------------------------------------
    rules = []
    for r in cfg["role_rules"]:
        r = dict(r)
        r["_patt"] = Chem.MolFromSmarts(r["smarts"])
        if r["_patt"] is None:
            raise ValueError(f"Bad SMARTS in role_rules: {r['smarts']}")
        rules.append(r)
    if rules:
        df["components"] = [[(_apply_role_rules(role, smi, rules), smi) for role, smi in comps]
                            for comps in df["components"]]

    # unique + sorted components => order-invariant, hashable reaction identity
    df["components"] = [sorted(set(map(tuple, c))) for c in df["components"]]

    # 5. element balance ---------------------------------------------------------------
    if cfg["element_balance"] != "off":
        bad = np.array([
            not element_set([p]) <= element_set([s for _, s in c])
            for c, p in zip(df["components"], df["product"])
        ])
        df["element_mismatch"] = bad
        before = len(df)
        if cfg["element_balance"] == "drop":
            df = df[~bad].reset_index(drop=True)
        log("element balance", before, len(df),
            f"{int(bad.sum())} reactions have product elements absent from inputs "
            f"({cfg['element_balance']})")

    # 6. de-duplicate ---------------------------------------------------------------------
    df["yield_std"] = 0.0
    df["n_replicates"] = 1
    if cfg["dedupe"] != "off":
        before = len(df)
        key = [(tuple(c), p) for c, p in zip(df["components"], df["product"])]
        df["_key"] = pd.factorize(pd.Series(key, dtype=object))[0]
        agg_fn = {"mean": "mean", "median": "median", "first": "first"}[cfg["dedupe"]]
        grouped = df.groupby("_key", sort=False)["yield"]
        agg = grouped.agg(agg_fn)
        std = grouped.std(ddof=0)
        cnt = grouped.size()
        first = df.drop_duplicates("_key").set_index("_key")
        first["yield"] = agg
        first["yield_std"] = std
        first["n_replicates"] = cnt
        n_conflict = int((std > cfg["conflict_std"]).sum())
        df = first.reset_index(drop=True)
        if cfg["drop_conflicts"]:
            df = df[df["yield_std"] <= cfg["conflict_std"]].reset_index(drop=True)
        log("duplicates", before, len(df),
            f"merged by {cfg['dedupe']}; {n_conflict} groups with replicate std > {cfg['conflict_std']}"
            + (" (dropped)" if cfg["drop_conflicts"] else " (kept, see yield_std)"))

    df = df.reset_index(drop=True)
    roles: dict[str, set[str]] = {}
    for comps in df["components"]:
        for role, smi in comps:
            roles.setdefault(role, set()).add(smi)
    from collections import Counter

    sig = Counter(tuple(sorted(Counter(r for r, _ in c).items())) for c in df["components"])
    report = {
        "role_signatures": [
            {"n_reactions": n, "roles": dict(k)} for k, n in sig.most_common(6)],
        "n_raw": n_raw,
        "n_clean": len(df),
        "steps": steps,
        "yield_before": stats_before,
        "yield_after": df["yield"].describe().round(2).to_dict(),
        "unique_components_by_role": {r: len(s) for r, s in sorted(roles.items())},
        "n_unique_products": int(df["product"].nunique()),
    }
    return df, report


# ------------------------------------------------------------------ persistence
def save_clean(df: pd.DataFrame, path: str | Path):
    out = df.copy()
    out["components"] = out["components"].map(json.dumps)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)


def load_clean(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["components"] = df["components"].map(lambda s: [tuple(x) for x in json.loads(s)])
    df["reaction_id"] = df["reaction_id"].astype(str)
    return df


def write_report(report: dict, path_md: str | Path):
    lines = ["# Cleaning report", "",
             f"Raw rows: **{report['n_raw']}** -> clean rows: **{report['n_clean']}**", "",
             "| step | removed | remaining | note |", "|---|---:|---:|---|"]
    for s in report["steps"]:
        lines.append(f"| {s['step']} | {s['removed']} | {s['remaining']} | {s['note']} |")
    lines += ["", "## Yield distribution (before -> after)", "",
              f"before: `{report['yield_before']}`", "", f"after: `{report['yield_after']}`", "",
              "## Unique components per role", ""]
    lines += [f"- {r}: {n}" for r, n in report["unique_components_by_role"].items()]
    lines += ["", f"Unique products: {report['n_unique_products']}", "",
              "## Role signatures (components per role, most common first)", ""]
    lines += [f"- {s['n_reactions']} reactions: {s['roles']}" for s in report["role_signatures"]]
    Path(path_md).parent.mkdir(parents=True, exist_ok=True)
    Path(path_md).write_text("\n".join(lines) + "\n")
    Path(str(path_md).replace(".md", ".json")).write_text(json.dumps(report, indent=2, default=str))
