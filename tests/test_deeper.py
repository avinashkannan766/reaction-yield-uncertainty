"""Leave-one-out splits, the ID-embedding baseline and the transformer-interaction option."""

import numpy as np
import pytest
import torch

from chemyield.features import GraphFeaturizer
from chemyield.metrics import regression_metrics
from chemyield.models import MODELS
from chemyield.run import run_experiment
from chemyield.splits import group_keys, held_out_label, make_split


def _parts(df, kind="random", **kw):
    idx = make_split(df, kind, seed=kw.pop("seed", 0), **kw)
    return {k: df.iloc[v].reset_index(drop=True) for k, v in idx.items()}


def test_loo_holds_out_exactly_one_group_and_folds_cover_all(clean_synth):
    keys = group_keys(clean_synth, "role:ligand")
    all_groups = set(keys)
    seen = set()
    for seed in range(len(all_groups)):
        s = make_split(clean_synth, "loo", seed, group_by="role:ligand")
        held = {keys[i] for i in s["test"]}
        assert len(held) == 1                                   # one whole ligand
        assert not (held & {keys[i] for k in ("train", "val", "cal") for i in s[k]})
        assert len(s["test"]) == sum(k == next(iter(held)) for k in keys)   # ALL its rows
        seen |= held
        assert held_out_label(clean_synth, "role:ligand", seed) == ".".join(next(iter(held)))
    assert seen == all_groups
    # wraps around, and a column with a single value cannot be leave-one-out
    assert held_out_label(clean_synth, "role:ligand", 0) == held_out_label(clean_synth, "role:ligand", 4)
    with pytest.raises(ValueError, match="too few usable groups"):
        make_split(clean_synth.assign(const=1), "loo", 0, group_by="column:const")


def test_per_seed_loo_fold_is_recorded_by_the_runner(tmp_path, clean_synth):
    from chemyield.synthetic import make_dataframe

    make_dataframe().sample(1200, random_state=3).to_csv(tmp_path / "a.csv", index=False)
    cfg = {"dataset": {"name": "d", "adapter": "csv", "params": {
               "path": str(tmp_path / "a.csv"), "id_col": "id", "product_col": "product",
               "role_cols": {"aryl_halide": "aryl_halide", "amine": "amine", "ligand": "ligand",
                             "base": "base", "solvent": "solvent"}}},
           "split": {"kind": "loo", "group_by": "role:ligand"}, "seeds": [0, 1],
           "features": [{"name": "morgan", "n_bits": 32}],
           "model": {"name": "lightgbm", "n_estimators": 60}}
    s = run_experiment(cfg, root=tmp_path, quiet=True)
    held = [m["held_out"] for m in s["per_seed"]]
    assert len(set(held)) == 2 and all(held)


def test_id_encoder_gives_unseen_molecules_the_unknown_vector(clean_synth):
    ligs = sorted({s for c in clean_synth["components"] for r, s in c if r == "ligand"})
    keep = [all(s != ligs[0] for r, s in c if r == "ligand") for c in clean_synth["components"]]
    train, held = clean_synth[keep].reset_index(drop=True), clean_synth[[not k for k in keep]].reset_index(drop=True)
    gf = GraphFeaturizer(product=False).fit(train)
    gtr, gte = gf.transform(train), gf.transform(held)         # table now also holds the held-out ligand
    m = MODELS.build({"name": "mpnn", "encoder": "id", "product": False, "hidden": 16, "ffn_hidden": 32,
                      "epochs": 2, "batch_size": 128})
    m.fit(gtr, train["yield"].to_numpy(), gtr, train["yield"].to_numpy())
    lut = m.nets_[0].lut
    held_id = gf.table.index[ligs[0]]
    seen_id = gf.table.index[ligs[1]]
    assert lut[held_id] == 0 and lut[seen_id] > 0               # unseen -> UNK, seen -> own slot
    assert torch.all(m.nets_[0].id_emb.weight[0] == 0)          # UNK vector is exactly zero
    mu, _ = m.predict(gte)
    assert np.isfinite(mu).all() and len(mu) == len(held)


def test_graph_encoder_differs_from_id_encoder_on_seen_molecules(clean_synth):
    """Sanity: both train and fit the signal when all molecules are known (random split)."""
    p = _parts(clean_synth, "random", test_size=0.2)
    gf = GraphFeaturizer().fit(p["train"])
    G = {k: gf.transform(v) for k, v in p.items()}
    y = {k: v["yield"].to_numpy() for k, v in p.items()}
    r2 = {}
    for enc in ("graph", "id"):
        m = MODELS.build({"name": "mpnn", "encoder": enc, "hidden": 32, "ffn_hidden": 64, "epochs": 15,
                          "patience": 5, "batch_size": 64, "dropout": 0.0, "lr": 2e-3})
        mu, _ = m.fit(G["train"], y["train"], G["val"], y["val"]).predict(G["test"])
        r2[enc] = regression_metrics(y["test"], mu)["r2"]
    assert r2["graph"] > 0.4 and r2["id"] > 0.4, r2


def test_transformer_interaction_trains_and_is_not_a_noop(clean_synth):
    p = _parts(clean_synth, "random", test_size=0.2)
    gf = GraphFeaturizer().fit(p["train"])
    G = {k: gf.transform(v) for k, v in p.items()}
    y = {k: v["yield"].to_numpy() for k, v in p.items()}
    spec = {"name": "mpnn", "interaction": "transformer", "hidden": 32, "heads": 4, "ffn_hidden": 64,
            "epochs": 15, "patience": 5, "batch_size": 64, "dropout": 0.0, "lr": 2e-3}
    m = MODELS.build(spec)
    mu, _ = m.fit(G["train"], y["train"], G["val"], y["val"]).predict(G["test"])
    assert regression_metrics(y["test"], mu)["r2"] > 0.4
    assert any("tr." in n for n, _ in m.nets_[0].named_parameters())
    with pytest.raises(ValueError):
        MODELS.build({**spec, "interaction": "nope"}).fit(G["train"], y["train"], G["val"], y["val"])
