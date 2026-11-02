"""Config-driven experiment runner.

    python -m chemyield.run --config configs/tdc_buchwald.yaml
    python -m chemyield.run --config configs/tdc_buchwald.yaml --set model.name=lightgbm split.kind=random

Pipeline: adapter -> clean (cached) -> split -> featurize (fit on train only) -> fit -> predict
-> conformal intervals -> metrics.  Swapping dataset = swapping the config's `dataset` block.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import features as _features  # noqa: F401  (registers featurizers)
from . import models as _models      # noqa: F401  (registers models)
from .adapters import build_adapter
from .cleaning import clean_reactions, load_clean, save_clean, write_report
from .config import load_config
from .features import FeatureSet, GraphFeaturizer
from .metrics import interval_metrics, regression_metrics
from .models import MODELS
from .splits import group_keys, held_out_label, make_split
from .uncertainty import SplitConformal


# ----------------------------------------------------------------------------- data
def prepare_data(cfg: dict, force: bool = False, root: str | Path = ".") -> pd.DataFrame:
    name = cfg["dataset"]["name"]
    out_dir = Path(root) / "data" / "processed" / name
    clean_path = out_dir / "clean.csv"
    if clean_path.exists() and not force:
        return load_clean(clean_path)
    raw = build_adapter(cfg["dataset"]).load()
    smp = cfg["dataset"].get("sample")       # {n: 50000, seed: 0}: random subset of the RAW rows
    if smp and smp["n"] < len(raw):
        raw = raw.sample(smp["n"], random_state=smp.get("seed", 0)).reset_index(drop=True)
        print(f"[sample] using {smp['n']} random raw rows (seed {smp.get('seed', 0)})")
    df, report = clean_reactions(raw, cfg.get("cleaning"))
    save_clean(df, clean_path)
    write_report(report, out_dir / "cleaning_report.md")
    print(f"[clean] {report['n_raw']} raw -> {report['n_clean']} clean rows "
          f"(report: {out_dir / 'cleaning_report.md'})")
    return df


# ------------------------------------------------------------------------ experiment
def _fmt(v):
    return "nan" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.3f}"


def _fit_predict(cfg: dict, seed: int, train: pd.DataFrame, val: pd.DataFrame,
                 targets: dict[str, pd.DataFrame]) -> dict:
    """Fit features + model on `train` (early stopping on `val`), predict every frame in `targets`.
    Featurizers are fit on `train` only, so each refit sees held-out chemistry as deployment would."""
    model_spec = dict(cfg["model"])
    model_spec.setdefault("seed", seed)
    model = MODELS.build(model_spec)
    fs = FeatureSet(cfg.get("features", [])).fit(train)
    if model.input_kind == "graph":
        gf = GraphFeaturizer(roles=model.roles, product=model.product).fit(train)

        def feats(d):
            extra = fs.transform(d) if fs.featurizers else None
            return gf.transform(d, extra=extra)
    else:
        def feats(d):
            return fs.transform(d)

    model.fit(feats(train), train["yield"].to_numpy(float), feats(val), val["yield"].to_numpy(float))
    return {k: model.predict(feats(d)) for k, d in targets.items()}


def _group_kfold_scores(df: pd.DataFrame, cfg: dict, seed: int, train: pd.DataFrame,
                        val: pd.DataFrame, group_by: str, k: int):
    """Out-of-group calibration data: split TRAIN's groups into K folds, refit on K-1 folds, predict
    the held-out fold. Every train row gets a prediction from a model that never saw its group."""
    keys = pd.Series([str(g) for g in group_keys(train, group_by)])
    uniq = sorted(keys.unique())
    k = min(k, len(uniq))
    if k < 2:
        raise ValueError(f"group_by='{group_by}' has too few usable groups ({len(uniq)} in train) "
                         "for group-K-fold calibration (needs >=2)")
    rng = np.random.default_rng(10_000 + seed)
    fold_of = {g: i % k for i, g in enumerate(rng.permutation(uniq))}
    fold = keys.map(fold_of).to_numpy()
    ys, mus, sigmas = [], [], []
    for f in range(k):
        tr, ho = train[fold != f].reset_index(drop=True), train[fold == f].reset_index(drop=True)
        mu_f, sg_f = _fit_predict(cfg, seed * 100 + f + 1, tr, val, {"ho": ho})["ho"]
        ys.append(ho["yield"].to_numpy(float))
        mus.append(mu_f)
        sigmas.append(sg_f)
    sigma = None if any(s is None for s in sigmas) else np.concatenate(sigmas)
    return np.concatenate(ys), np.concatenate(mus), sigma


def run_seed(df: pd.DataFrame, cfg: dict, seed: int) -> tuple[dict, pd.DataFrame]:
    sp = cfg.get("split", {})
    cf = cfg.get("conformal", {})
    method = cf.get("method", "split")
    if method not in ("split", "group_kfold"):
        raise ValueError(f"conformal.method must be 'split' or 'group_kfold', got '{method}'")
    idx = make_split(df, kind=sp.get("kind", "random"), seed=seed,
                     test_size=sp.get("test_size", 0.2), val_size=sp.get("val_size", 0.1),
                     cal_size=sp.get("cal_size", 0.1), group_by=sp.get("group_by"),
                     val_by_group=sp.get("val_by_group", False))
    if method == "group_kfold":
        if sp.get("kind") not in ("group", "loo"):
            raise ValueError("conformal.method=group_kfold needs split.kind 'group' or 'loo'")
        # calibration comes from out-of-group refits of the train data, so the cal rows go back
        # into training: no training groups are given up
        idx["train"] = np.sort(np.concatenate([idx["train"], idx["cal"]]))
    part = {k: df.iloc[v].reset_index(drop=True) for k, v in idx.items()}
    y = {k: v["yield"].to_numpy(float) for k, v in part.items()}

    t0 = time.time()
    pred = _fit_predict(cfg, seed, part["train"], part["val"], {"cal": part["cal"], "test": part["test"]})
    fit_seconds = time.time() - t0

    alpha = cf.get("alpha", 0.1)
    if method == "group_kfold":
        y_cal, mu_cal, sg_cal = _group_kfold_scores(df, cfg, seed, part["train"], part["val"],
                                                    sp["group_by"], int(cf.get("k", 5)))
        conf = SplitConformal(alpha=alpha).fit(y_cal, mu_cal, sg_cal)
    else:
        conf = SplitConformal(alpha=alpha).fit(y["cal"], *pred["cal"])
    mu, sigma = pred["test"]
    lo, hi = conf.interval(mu, sigma)

    held = held_out_label(df, sp["group_by"], seed) if sp.get("kind") == "loo" else None
    metrics = {"seed": seed, "held_out": held, "n_train": len(idx["train"]), "n_test": len(idx["test"]),
               "fit_seconds": round(fit_seconds, 1),
               **regression_metrics(y["test"], mu),
               **interval_metrics(y["test"], lo, hi, sigma, mu), "target_coverage": 1 - alpha}
    out = pd.DataFrame({"reaction_id": part["test"]["reaction_id"], "y": y["test"], "mu": mu,
                        "sigma": sigma if sigma is not None else np.nan, "lo": lo, "hi": hi})
    return metrics, out


def run_experiment(cfg: dict, df: pd.DataFrame | None = None, root: str | Path = ".",
                   save: bool = True, quiet: bool = False) -> dict:
    df = df if df is not None else prepare_data(cfg, root=root)
    exp = cfg.get("experiment", {}).get("name", "exp")
    seeds = cfg.get("seeds", [0])
    out_dir = Path(root) / "runs" / cfg["dataset"]["name"] / exp
    per_seed = []
    for s in seeds:
        m, preds = run_seed(df, cfg, s)
        per_seed.append(m)
        if save:
            (out_dir / f"seed_{s}").mkdir(parents=True, exist_ok=True)
            preds.to_csv(out_dir / f"seed_{s}" / "predictions.csv", index=False)
            (out_dir / f"seed_{s}" / "metrics.json").write_text(json.dumps(m, indent=2))
        if not quiet:
            print(f"[{exp}] seed {s}{' (held out ' + m['held_out'][:30] + ')' if m.get('held_out') else ''}: R2={_fmt(m['r2'])} RMSE={_fmt(m['rmse'])} "
                  f"MAE={_fmt(m['mae'])} cov={_fmt(m['coverage'])} ({m['fit_seconds']}s)")
    keys = ["r2", "rmse", "mae", "spearman", "coverage", "mean_width"]
    summary = {"experiment": exp, "n_seeds": len(seeds)}
    for k in keys:
        v = np.array([m.get(k, np.nan) for m in per_seed], float)
        summary[k] = ([float("nan"), float("nan")] if np.isnan(v).all()
                      else [float(np.nanmean(v)), float(np.nanstd(v))])
    if save:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "config.yaml").write_text(yaml.safe_dump(
            {k: v for k, v in cfg.items() if not k.startswith("_")}, sort_keys=False))
        (out_dir / "summary.json").write_text(json.dumps({**summary, "per_seed": per_seed}, indent=2))
    return {**summary, "per_seed": per_seed}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--set", nargs="*", default=[], metavar="a.b=value", help="override config keys")
    ap.add_argument("--force-clean", action="store_true", help="re-run the cleaning step")
    args = ap.parse_args(argv)
    cfg = load_config(args.config, args.set)
    df = prepare_data(cfg, force=args.force_clean)
    s = run_experiment(cfg, df)
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
