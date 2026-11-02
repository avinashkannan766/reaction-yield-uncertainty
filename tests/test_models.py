import numpy as np
import pytest

from chemyield.features import FeatureSet, GraphFeaturizer
from chemyield.metrics import regression_metrics
from chemyield.models import MODELS
from chemyield.splits import make_split
from chemyield.uncertainty import SplitConformal


def _split(df):
    idx = make_split(df, "random", seed=0, test_size=0.2)
    return {k: df.iloc[v].reset_index(drop=True) for k, v in idx.items()}


@pytest.mark.parametrize("spec", [{"name": "ridge"}, {"name": "rf", "n_estimators": 60},
                                  {"name": "lightgbm", "n_estimators": 200, "n_models": 2}])
def test_tabular_models_beat_mean_baseline(clean_synth, spec):
    p = _split(clean_synth)
    fs = FeatureSet([{"name": "morgan", "n_bits": 64}]).fit(p["train"])
    X = {k: fs.transform(v) for k, v in p.items()}
    y = {k: v["yield"].to_numpy() for k, v in p.items()}
    mu, sigma = MODELS.build(spec).fit(X["train"], y["train"], X["val"], y["val"]).predict(X["test"])
    base = MODELS.build({"name": "mean"}).fit(X["train"], y["train"]).predict(X["test"])[0]
    assert regression_metrics(y["test"], mu)["rmse"] < 0.8 * regression_metrics(y["test"], base)["rmse"]
    assert sigma is None or (sigma >= 0).all()


def test_mpnn_learns_the_hidden_signal(clean_synth):
    """The D-MPNN must clearly beat the mean baseline and give finite ensemble sigma."""
    p = _split(clean_synth)
    model = MODELS.build({"name": "mpnn", "hidden": 48, "ffn_hidden": 96, "epochs": 15,
                          "patience": 5, "n_models": 2, "batch_size": 64, "dropout": 0.0, "lr": 2e-3})
    gf = GraphFeaturizer().fit(p["train"])
    G = {k: gf.transform(v) for k, v in p.items()}
    y = {k: v["yield"].to_numpy() for k, v in p.items()}
    model.fit(G["train"], y["train"], G["val"], y["val"])
    mu, sigma = model.predict(G["test"])
    m = regression_metrics(y["test"], mu)
    assert m["r2"] > 0.5, m
    assert sigma is not None and np.isfinite(sigma).all()


def test_mpnn_extra_features_path(clean_synth):
    p = _split(clean_synth)
    fs = FeatureSet([{"name": "conditions", "numeric_cols": ["temperature"]}]).fit(p["train"])
    gf = GraphFeaturizer().fit(p["train"])
    G = {k: gf.transform(v, extra=fs.transform(v)) for k, v in p.items()}
    y = {k: v["yield"].to_numpy() for k, v in p.items()}
    model = MODELS.build({"name": "mpnn", "hidden": 16, "ffn_hidden": 32, "epochs": 2, "batch_size": 128})
    mu, _ = model.fit(G["train"], y["train"], G["val"], y["val"]).predict(G["test"])
    assert mu.shape == (len(p["test"]),) and np.isfinite(mu).all()


def test_conformal_coverage_on_exchangeable_data():
    rng = np.random.default_rng(0)
    y = rng.normal(50, 10, 4000)
    mu = y + rng.normal(0, 5, y.size)
    conf = SplitConformal(alpha=0.1, lo=-1e9, hi=1e9).fit(y[:2000], mu[:2000])
    lo, hi = conf.interval(mu[2000:])
    assert abs(((y[2000:] >= lo) & (y[2000:] <= hi)).mean() - 0.9) < 0.03
