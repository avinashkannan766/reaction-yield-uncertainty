"""Featurizer interface + registry. Tabular featurizers return (n_rows, n_features) arrays."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..registry import Registry

FEATURIZERS = Registry("featurizer")


def roles_in(df: pd.DataFrame) -> list[str]:
    return sorted({role for comps in df["components"] for role, _ in comps})


def unique_smiles(df: pd.DataFrame) -> list[str]:
    seen = {smi for comps in df["components"] for _, smi in comps}
    seen |= {s for p in df["product"] for s in p.split(".")}
    return sorted(seen)


class Featurizer:
    names_: list[str]

    def fit(self, df: pd.DataFrame):
        return self

    def transform(self, df: pd.DataFrame) -> np.ndarray:  # pragma: no cover - interface
        raise NotImplementedError


class FeatureSet:
    """Concatenates several featurizers built from config specs."""

    def __init__(self, specs: list[dict]):
        self.featurizers = [FEATURIZERS.build(s) for s in specs]

    def fit(self, df):
        for f in self.featurizers:
            f.fit(df)
        return self

    def transform(self, df) -> np.ndarray:
        if not self.featurizers:
            return np.zeros((len(df), 0), dtype=np.float32)
        return np.concatenate([f.transform(df) for f in self.featurizers], axis=1).astype(np.float32)

    def fit_transform(self, df):
        return self.fit(df).transform(df)

    @property
    def names(self) -> list[str]:
        return [n for f in self.featurizers for n in f.names_]
