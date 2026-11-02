"""conformal.method=group_kfold: calibrate on out-of-group residuals from K group-wise refits."""

import numpy as np
import pytest

from chemyield.run import run_seed, _group_kfold_scores, _fit_predict
from chemyield.splits import group_keys, make_split


def _cfg(method, kind="group", group_by="role:aryl_halide", **conf):
    return {"split": {"kind": kind, "group_by": group_by, "test_size": 0.2},
            "features": [{"name": "morgan", "n_bits": 64}],
            "model": {"name": "lightgbm", "n_estimators": 100},
            "conformal": {"alpha": 0.1, "method": method, **conf}}


def test_scores_cover_every_train_row_once_from_models_that_never_saw_the_group(clean_synth, monkeypatch):
    import chemyield.run as R
    idx = make_split(clean_synth, "group", 0, test_size=0.2, group_by="role:aryl_halide")
    train = clean_synth.iloc[idx["train"]].reset_index(drop=True)
    val = clean_synth.iloc[idx["val"]].reset_index(drop=True)
    seen = []
    real = R._fit_predict

    def spy(cfg, seed, tr, va, targets):
        ho = targets["ho"]
        g_tr = set(group_keys(tr, "role:aryl_halide")); g_ho = set(group_keys(ho, "role:aryl_halide"))
        assert not (g_tr & g_ho)                      # held-out fold shares no group with its fit
        seen.append(len(ho))
        return real(cfg, seed, tr, va, targets)

    monkeypatch.setattr(R, "_fit_predict", spy)
    y, mu, sg = _group_kfold_scores(clean_synth, _cfg("group_kfold"), 0, train, val, "role:aryl_halide", 4)
    assert sum(seen) == len(train) == len(y) == len(mu)


def test_runs_and_widens_intervals_vs_random_row_calibration(clean_synth):
    cov, width = {}, {}
    for method in ("split", "group_kfold"):
        res = [run_seed(clean_synth, _cfg(method, k=4), s)[0] for s in range(3)]
        cov[method] = np.mean([m["coverage"] for m in res])
        width[method] = np.mean([m["mean_width"] for m in res])
    assert cov["group_kfold"] > cov["split"] and width["group_kfold"] > width["split"], (cov, width)


def test_works_where_val_by_group_cannot(clean_synth):
    m, _ = run_seed(clean_synth, _cfg("group_kfold", kind="loo", group_by="role:ligand", k=3), 0)
    assert np.isfinite(m["coverage"])                 # 4 ligands: 3 left in train -> 3 folds


def test_needs_group_split_and_known_method(clean_synth):
    with pytest.raises(ValueError, match="group"):
        run_seed(clean_synth, {**_cfg("group_kfold"), "split": {"kind": "random"}}, 0)
    with pytest.raises(ValueError, match="method"):
        run_seed(clean_synth, _cfg("bogus"), 0)
