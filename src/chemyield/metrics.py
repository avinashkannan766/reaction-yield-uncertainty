"""Point-prediction and interval metrics."""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def regression_metrics(y, mu) -> dict:
    y, mu = np.asarray(y, float), np.asarray(mu, float)
    err = mu - y
    ss_res, ss_tot = float((err ** 2).sum()), float(((y - y.mean()) ** 2).sum())
    return {
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "rmse": float(np.sqrt((err ** 2).mean())),
        "mae": float(np.abs(err).mean()),
        "spearman": float(spearmanr(y, mu).statistic) if len(y) > 2 and np.ptp(mu) > 0 else float("nan"),
    }


def interval_metrics(y, lo, hi, sigma=None, mu=None) -> dict:
    y, lo, hi = (np.asarray(a, float) for a in (y, lo, hi))
    out = {"coverage": float(((y >= lo) & (y <= hi)).mean()), "mean_width": float((hi - lo).mean())}
    if sigma is not None and mu is not None and np.ptp(sigma) > 0:
        out["err_sigma_spearman"] = float(
            spearmanr(np.abs(np.asarray(mu) - y), np.asarray(sigma)).statistic)
    return out
