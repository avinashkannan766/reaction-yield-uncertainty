"""val_by_group: validation + calibration sets that are whole unseen groups."""

import numpy as np
import pytest

from chemyield.splits import group_keys, make_split


def _split(df, seed=0, **kw):
    return make_split(df, "group", seed, test_size=0.2, val_size=0.1, cal_size=0.1,
                      group_by="role:aryl_halide", **kw)


def test_val_cal_test_are_disjoint_whole_groups_unseen_in_train(clean_synth):
    keys = group_keys(clean_synth, "role:aryl_halide")
    s = _split(clean_synth, val_by_group=True)
    sets = {k: {keys[i] for i in s[k]} for k in ("train", "val", "cal", "test")}
    names = list(sets)
    for a in names:
        for b in names:
            if a < b:
                assert not (sets[a] & sets[b]), (a, b)          # no group shared between any two parts
    allidx = np.concatenate(list(s.values()))
    assert len(allidx) == len(set(allidx)) == len(clean_synth)  # nothing lost or duplicated
    n = len(clean_synth)
    assert 0.05 * n < len(s["val"]) < 0.25 * n and 0.05 * n < len(s["cal"]) < 0.25 * n
    assert len(s["train"]) > 0.3 * n


def test_default_is_unchanged_random_val_cal(clean_synth):
    keys = group_keys(clean_synth, "role:aryl_halide")
    s = _split(clean_synth)                                     # val_by_group off
    train_groups = {keys[i] for i in s["train"]}
    assert {keys[i] for i in s["val"]} <= train_groups          # val rows come from SEEN groups


def test_deterministic_and_seed_sensitive(clean_synth):
    a, b, c = (_split(clean_synth, seed=x, val_by_group=True) for x in (0, 0, 1))
    assert all(np.array_equal(a[k], b[k]) for k in a)
    assert not np.array_equal(a["val"], c["val"])


def test_too_few_groups_raises_the_skippable_error(clean_synth):
    with pytest.raises(ValueError, match="too few usable groups"):
        make_split(clean_synth, "loo", 0, group_by="role:ligand", val_by_group=True)   # 4 ligands
    with pytest.raises(ValueError, match="too few usable groups"):
        make_split(clean_synth, "group", 0, test_size=0.2, group_by="role:base", val_by_group=True)


def test_requires_group_or_loo_kind(clean_synth):
    with pytest.raises(ValueError, match="group"):
        make_split(clean_synth, "random", 0, val_by_group=True)


def test_group_calibration_raises_coverage_under_group_shift_and_widens_intervals(clean_synth):
    """The reason the option exists: random-row calibration under-covers a NEW group. Calibrating on
    whole unseen groups must give higher mean coverage (at the price of wider intervals)."""
    from chemyield.run import run_seed

    def run(gval):
        cov, width = [], []
        for seed in range(4):
            cfg = {"split": {"kind": "group", "group_by": "role:aryl_halide", "test_size": 0.2,
                             "val_by_group": gval},
                   "features": [{"name": "morgan", "n_bits": 64}],
                   "model": {"name": "lightgbm", "n_estimators": 150}, "conformal": {"alpha": 0.1}}
            m, _ = run_seed(clean_synth, cfg, seed)
            cov.append(m["coverage"])
            width.append(m["mean_width"])
        return np.mean(cov), np.mean(width)

    (cov_r, w_r), (cov_g, w_g) = run(False), run(True)
    assert cov_g > cov_r + 0.03, (cov_r, cov_g)
    assert w_g > w_r, (w_r, w_g)
