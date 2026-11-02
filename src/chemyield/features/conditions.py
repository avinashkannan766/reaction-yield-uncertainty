"""Condition features: numeric columns (temperature, time ...) and component identity."""

from __future__ import annotations

import numpy as np

from .base import FEATURIZERS, Featurizer


@FEATURIZERS.register("conditions")
class ConditionFeaturizer(Featurizer):
    """numeric_cols  : standardised with train statistics (constant columns become 0)
    onehot_roles  : multi-hot of which molecule fills each role (vocabulary from train only,
                    so a held-out ligand is all-zeros - exactly what the model will face)
    """

    def __init__(self, numeric_cols=(), onehot_roles=()):
        self.numeric_cols = list(numeric_cols)
        self.onehot_roles = list(onehot_roles)

    def fit(self, df):
        self.mean_, self.std_ = {}, {}
        for c in self.numeric_cols:
            v = df[c].astype(float)
            self.mean_[c] = float(v.mean())
            s = float(v.std(ddof=0))
            self.std_[c] = s if s > 0 else 1.0
        self.vocab_ = {}
        for r in self.onehot_roles:
            self.vocab_[r] = sorted({s for comps in df["components"] for role, s in comps if role == r})
        self.names_ = list(self.numeric_cols)
        self.names_ += [f"is_{r}_{i}" for r in self.onehot_roles for i in range(len(self.vocab_[r]))]
        return self

    def transform(self, df):
        cols = []
        for c in self.numeric_cols:
            v = (df[c].astype(float).to_numpy() - self.mean_[c]) / self.std_[c]
            cols.append(v[:, None])
        for r in self.onehot_roles:
            pos = {s: i for i, s in enumerate(self.vocab_[r])}
            block = np.zeros((len(df), len(pos)), np.float32)
            for row, comps in enumerate(df["components"]):
                for role, smi in comps:
                    if role == r and smi in pos:
                        block[row, pos[smi]] = 1.0
            cols.append(block)
        if not cols:
            return np.zeros((len(df), 0), np.float32)
        return np.concatenate(cols, axis=1).astype(np.float32)
