"""Split-conformal prediction intervals on top of any regressor."""

from __future__ import annotations

import numpy as np


class SplitConformal:
    """Normalised if the model supplies sigma (e.g. ensemble std), else plain |residual|.

    Guarantee holds only if calibration and test are exchangeable; under group
    (out-of-distribution) splits expect under-coverage - report it, don't hide it.
    """

    def __init__(self, alpha: float = 0.1, lo: float = 0.0, hi: float = 100.0, eps: float = 1.0):
        self.alpha, self.lo, self.hi, self.eps = alpha, lo, hi, eps
        self.q_: float | None = None

    def fit(self, y_cal, mu_cal, sigma_cal=None):
        scores = np.abs(np.asarray(y_cal) - np.asarray(mu_cal))
        if sigma_cal is not None:
            scores = scores / (np.asarray(sigma_cal) + self.eps)
        n = len(scores)
        level = min(1.0, np.ceil((n + 1) * (1 - self.alpha)) / n)
        self.q_ = float(np.quantile(scores, level, method="higher"))
        return self

    def interval(self, mu, sigma=None):
        mu = np.asarray(mu, float)
        width = self.q_ * (np.asarray(sigma) + self.eps) if sigma is not None else self.q_
        return np.clip(mu - width, self.lo, self.hi), np.clip(mu + width, self.lo, self.hi)
