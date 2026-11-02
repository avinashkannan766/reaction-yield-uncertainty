"""Model interface. `input_kind` tells the runner what to feed: 'tabular' (ndarray) or 'graph'."""

from __future__ import annotations

from ..registry import Registry

MODELS = Registry("model")


class Model:
    input_kind = "tabular"

    def fit(self, X, y, X_val=None, y_val=None):  # pragma: no cover - interface
        raise NotImplementedError

    def predict(self, X):
        """Return (mu, sigma). sigma is None when the model has no native uncertainty."""
        raise NotImplementedError
