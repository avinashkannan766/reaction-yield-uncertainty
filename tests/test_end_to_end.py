"""Whole pipeline from a config, and the dataset-swap guarantee."""

import numpy as np
import pandas as pd

from chemyield.config import deep_update, parse_overrides
from chemyield.run import run_experiment, prepare_data
from chemyield.synthetic import make_dataframe


def _cfg(csv_path, name):
    return {
        "dataset": {"name": name, "adapter": "csv", "params": {
            "path": str(csv_path), "id_col": "id", "product_col": "product", "yield_col": "yield",
            "role_cols": {"aryl_halide": "aryl_halide", "amine": "amine", "ligand": "ligand",
                          "base": "base", "solvent": "solvent"}, "keep_cols": ["temperature"]}},
        "split": {"kind": "group", "group_by": "role:aryl_halide", "test_size": 0.2},
        "seeds": [0],
        "features": [{"name": "morgan", "n_bits": 64}],
        "model": {"name": "lightgbm", "n_estimators": 150},
    }


def test_full_run_writes_outputs_and_reports_sane_metrics(tmp_path):
    make_dataframe().sample(1500, random_state=1).to_csv(tmp_path / "a.csv", index=False)
    cfg = _cfg(tmp_path / "a.csv", "ds_a")
    s = run_experiment(cfg, root=tmp_path, quiet=True)
    assert s["r2"][0] > 0.2
    run_dir = tmp_path / "runs" / "ds_a" / "exp"
    assert (run_dir / "summary.json").exists() and (run_dir / "seed_0" / "predictions.csv").exists()
    preds = pd.read_csv(run_dir / "seed_0" / "predictions.csv")
    assert (preds["lo"] <= preds["mu"]).all() and (preds["mu"] <= preds["hi"]).all()
    assert (tmp_path / "data/processed/ds_a/cleaning_report.md").exists()


def test_swapping_dataset_needs_only_the_dataset_block(tmp_path):
    """A differently-shaped source (reaction-SMILES column, 'reactant'/'reagent' roles) runs
    through the SAME model/feature/split code with only the config's dataset block changed."""
    df = make_dataframe().sample(1200, random_state=2)
    rxn = [f"{h}.{a}>{l}.{b}.{s}>{p}" for h, a, l, b, s, p in
           zip(df.aryl_halide, df.amine, df.ligand, df.base, df.solvent, df["product"])]
    pd.DataFrame({"rid": range(len(df)), "rxn": rxn, "y": df["yield"].to_numpy() / 100.0}).to_csv(
        tmp_path / "b.csv", index=False)       # note: yields as fractions, different column names
    cfg = _cfg(tmp_path / "a.csv", "ds_b")
    cfg["dataset"] = {"name": "ds_b", "adapter": "csv", "params": {
        "path": str(tmp_path / "b.csv"), "id_col": "rid", "yield_col": "y", "reaction_smiles_col": "rxn"}}
    cfg["split"] = {"kind": "random"}
    out = prepare_data(cfg, root=tmp_path)
    assert out["yield"].max() > 1.5                         # fraction -> percent was auto-detected
    assert {r for c in out["components"] for r, _ in c} == {"reactant", "reagent"}
    s = run_experiment(cfg, out, root=tmp_path, quiet=True)
    assert s["r2"][0] > 0.2


def test_override_parsing():
    o = parse_overrides(["split.kind=scaffold", "seeds=[0,1]", "model.n_models=3"])
    assert o == {"split": {"kind": "scaffold"}, "seeds": [0, 1], "model": {"n_models": 3}}
    assert deep_update({"split": {"kind": "random", "test_size": 0.2}}, o)["split"] == {
        "kind": "scaffold", "test_size": 0.2}


def test_dataset_sample_option_limits_raw_rows(tmp_path):
    make_dataframe().sample(900, random_state=1).to_csv(tmp_path / "a.csv", index=False)
    cfg = _cfg(tmp_path / "a.csv", "ds_s")
    cfg["dataset"]["sample"] = {"n": 300, "seed": 0}
    out = prepare_data(cfg, root=tmp_path)
    assert len(out) <= 300
