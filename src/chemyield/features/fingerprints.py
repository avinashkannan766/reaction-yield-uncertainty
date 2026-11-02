"""Morgan (ECFP-like) count fingerprints, pooled per component role."""

from __future__ import annotations

import numpy as np
from rdkit.Chem import rdFingerprintGenerator

from ..chem import cached_mol
from .base import FEATURIZERS, Featurizer, roles_in


@FEATURIZERS.register("morgan")
class MorganFeaturizer(Featurizer):
    """One block per role (sum of component fingerprints) + a block for the product.

    Role-wise blocks keep 'which molecule played which part' explicit, which a single
    fingerprint of the whole reaction mixture would lose.
    """

    def __init__(self, radius: int = 2, n_bits: int = 512, count: bool = True,
                 roles="auto", product: bool = True):
        self.radius, self.n_bits, self.count = radius, n_bits, count
        self.roles_cfg, self.product = roles, product
        self._gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        self._cache: dict[str, np.ndarray] = {}

    def fit(self, df):
        self.roles_ = roles_in(df) if self.roles_cfg == "auto" else list(self.roles_cfg)
        blocks = list(self.roles_) + (["product"] if self.product else [])
        self.names_ = [f"morgan_{b}_{i}" for b in blocks for i in range(self.n_bits)]
        return self

    def _fp(self, smi: str) -> np.ndarray:
        if smi not in self._cache:
            mol = cached_mol(smi)
            if mol is None:
                v = np.zeros(self.n_bits, np.float32)
            elif self.count:
                v = self._gen.GetCountFingerprintAsNumPy(mol).astype(np.float32)
            else:
                v = self._gen.GetFingerprintAsNumPy(mol).astype(np.float32)
            self._cache[smi] = v
        return self._cache[smi]

    def transform(self, df):
        idx = {r: i for i, r in enumerate(self.roles_)}
        n_blocks = len(self.roles_) + (1 if self.product else 0)
        X = np.zeros((len(df), n_blocks, self.n_bits), np.float32)
        for row, (comps, prod) in enumerate(zip(df["components"], df["product"])):
            for role, smi in comps:
                if role in idx:
                    X[row, idx[role]] += self._fp(smi)
            if self.product:
                for smi in prod.split("."):
                    X[row, -1] += self._fp(smi)
        return X.reshape(len(df), -1)
