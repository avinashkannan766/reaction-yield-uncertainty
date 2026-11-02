"""Run a grid of experiments (feature/model ablations x split types) and write a results table.

    python -m chemyield.benchmark --config configs/tdc_buchwald.yaml

The config's `experiments:` list overrides `features` / `model` per row; `benchmark.splits` lists
the split definitions to evaluate each row under. Cleaning runs once and is cached.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .config import deep_update, load_config
from .run import prepare_data, run_experiment


def _pm(v):
    return f"{v[0]:.3f} ± {v[1]:.3f}"


def run_benchmark(cfg: dict, root: str | Path = ".", only: list[str] | None = None,
                  only_splits: list[str] | None = None, suffix: str = "") -> pd.DataFrame:
    df = prepare_data(cfg, root=root)
    splits = cfg.get("benchmark", {}).get("splits") or [cfg.get("split", {"kind": "random"})]
    rows, folds = [], []
    for split in splits:
        split_name = split.get("name") or (split["kind"] + (f":{split['group_by']}" if split.get("group_by") else ""))
        if split.get("val_by_group"):
            split_name += "+gval"          # val/cal are whole unseen groups too
        if split.get("conformal", {}).get("method") == "group_kfold":
            split_name += "+gkf"           # calibration from group-K-fold out-of-group residuals
        if only_splits and not any(t in split_name for t in only_splits):
            continue
        for exp in cfg["experiments"]:
            if only and exp["name"] not in only:
                continue
            run_cfg = deep_update(cfg, {k: v for k, v in exp.items() if k != "name"})
            run_cfg["split"] = {k: v for k, v in split.items() if k != "conformal"}
            run_cfg["conformal"] = {**cfg.get("conformal", {}), **split.get("conformal", {})}
            run_cfg["seeds"] = split.get("seeds", cfg.get("seeds", [0]))
            run_cfg["experiment"] = {"name": f"{split_name.replace(':', '-')}__{exp['name']}"}
            try:
                s = run_experiment(run_cfg, df, root=root, quiet=True)
            except ValueError as e:
                if "too few usable groups" not in str(e):
                    raise
                print(f"[{split_name}] SKIPPED: {e}")
                break
            print(f"[{split_name}] {exp['name']:<28} R2 {_pm(s['r2'])}  RMSE {_pm(s['rmse'])}")
            if split.get("kind") == "loo":
                for m in s["per_seed"]:
                    folds.append({"split": split_name, "experiment": exp["name"],
                                  "held_out": m["held_out"], "n_test": m["n_test"],
                                  "r2": m["r2"], "rmse": m["rmse"], "mae": m["mae"],
                                  "coverage": m["coverage"]})
            rows.append({"split": split_name, "experiment": exp["name"],
                         **{k: s[k] for k in ("r2", "rmse", "mae", "spearman", "coverage", "mean_width")}})
    table = pd.DataFrame(rows)
    out = Path(root) / "results"
    out.mkdir(exist_ok=True)
    name = cfg["dataset"]["name"] + suffix
    flat = table.copy()
    for k in ("r2", "rmse", "mae", "spearman", "coverage", "mean_width"):
        flat[k + "_mean"], flat[k + "_std"] = flat[k].str[0], flat[k].str[1]
        flat = flat.drop(columns=k)
    flat.to_csv(out / f"{name}_benchmark.csv", index=False)
    lines = [f"# {name}: benchmark (mean ± std over seeds / leave-one-out folds)", ""]
    for split_name, g in table.groupby("split", sort=False):
        lines += [f"## split: {split_name}", "",
                  "| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |",
                  "|---|---|---|---|---|---|---|"]
        for _, r in g.iterrows():
            lines.append(f"| {r['experiment']} | {_pm(r['r2'])} | {_pm(r['rmse'])} | {_pm(r['mae'])} | "
                         f"{_pm(r['spearman'])} | {_pm(r['coverage'])} | {_pm(r['mean_width'])} |")
        lines.append("")
    if folds:
        fd = pd.DataFrame(folds)
        fd.to_csv(out / f"{name}_loo_folds.csv", index=False)
        lines += ["## leave-one-out folds (R² per held-out group)", "",
                  "| split | experiment | held out | n_test | R² | RMSE | coverage |",
                  "|---|---|---|---|---|---|---|"]
        for _, r in fd.iterrows():
            lines.append(f"| {r['split']} | {r['experiment']} | {str(r['held_out'])[:40]} | {r['n_test']} | "
                         f"{r['r2']:.3f} | {r['rmse']:.2f} | {r['coverage']:.2f} |")
        lines.append("")
    (out / f"{name}_benchmark.md").write_text("\n".join(lines))
    print(f"\nWrote {out / f'{name}_benchmark.md'}")
    return table


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--only", nargs="*", help="run only these experiment names")
    ap.add_argument("--only-splits", nargs="*", help="run only splits whose name contains one of these (e.g. +gkf)")
    ap.add_argument("--suffix", default="", help="appended to the results file names (keeps a partial run from overwriting the full one)")
    args = ap.parse_args(argv)
    run_benchmark(load_config(args.config, args.set), only=args.only,
                  only_splits=args.only_splits, suffix=args.suffix)


if __name__ == "__main__":
    main()
