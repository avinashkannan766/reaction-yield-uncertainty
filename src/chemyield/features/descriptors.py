"""Physchem descriptors (RDKit) and user-supplied lookup tables, pooled per role."""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit.Chem import Descriptors

from ..chem import canonical_smiles, cached_mol
from .base import FEATURIZERS, Featurizer, roles_in

DEFAULT_DESCRIPTORS = [
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds",
    "RingCount", "NumAromaticRings", "FractionCSP3", "HeavyAtomCount", "MolMR",
    "LabuteASA", "BertzCT", "NumHeteroatoms",
]
_DESC_FUNCS = dict(Descriptors._descList)


@FEATURIZERS.register("rdkit_desc")
class RDKitDescriptorFeaturizer(Featurizer):
    """Mean descriptor value per role (+ component count) and for the product. NaN if absent."""

    def __init__(self, descriptors=None, roles="auto", product: bool = True):
        self.descriptors = list(descriptors or DEFAULT_DESCRIPTORS)
        unknown = [d for d in self.descriptors if d not in _DESC_FUNCS]
        if unknown:
            raise ValueError(f"Unknown RDKit descriptors: {unknown}")
        self.roles_cfg, self.product = roles, product
        self._cache: dict[str, np.ndarray] = {}

    def fit(self, df):
        self.roles_ = roles_in(df) if self.roles_cfg == "auto" else list(self.roles_cfg)
        blocks = list(self.roles_) + (["product"] if self.product else [])
        self.names_ = [f"desc_{b}_{d}" for b in blocks for d in self.descriptors]
        self.names_ += [f"n_{r}" for r in self.roles_]
        return self

    def _desc(self, smi):
        if smi not in self._cache:
            mol = cached_mol(smi)
            vals = []
            for d in self.descriptors:
                try:
                    vals.append(float(_DESC_FUNCS[d](mol)))
                except Exception:
                    vals.append(np.nan)
            self._cache[smi] = np.nan_to_num(np.array(vals, np.float32), nan=np.nan,
                                             posinf=np.nan, neginf=np.nan)
        return self._cache[smi]

    def transform(self, df):
        idx = {r: i for i, r in enumerate(self.roles_)}
        nb = len(self.roles_) + (1 if self.product else 0)
        D = len(self.descriptors)
        out = np.full((len(df), nb, D), np.nan, np.float32)
        counts = np.zeros((len(df), len(self.roles_)), np.float32)
        for row, (comps, prod) in enumerate(zip(df["components"], df["product"])):
            per_role: dict[int, list] = {}
            for role, smi in comps:
                if role in idx:
                    per_role.setdefault(idx[role], []).append(self._desc(smi))
                    counts[row, idx[role]] += 1
            for b, vs in per_role.items():
                out[row, b] = np.nanmean(np.stack(vs), axis=0) if vs else np.nan
            if self.product:
                out[row, -1] = np.nanmean(np.stack([self._desc(s) for s in prod.split(".")]), axis=0)
        return np.concatenate([out.reshape(len(df), -1), counts], axis=1)


@FEATURIZERS.register("lookup")
class LookupFeaturizer(Featurizer):
    """Join a user table of physical descriptors onto components, mean-pooled per role.

    This is the hook for the 'condition embedding' ideas: Kamlet-Taft solvent parameters,
    base pKaH, ligand %Vbur / Tolman angle, Boltzmann-averaged 3D descriptors, ...
    Table = CSV with a SMILES column (`key_col`) plus numeric columns.
    Components missing from the table give NaN (LightGBM handles it; other models impute).
    """

    def __init__(self, table: str, roles, key_col: str = "smiles", columns=None, prefix="lk"):
        tab = pd.read_csv(table)
        keys = [canonical_smiles(s)[0] for s in tab[key_col]]
        tab = tab.assign(_key=keys).dropna(subset=["_key"]).drop_duplicates("_key").set_index("_key")
        cols = list(columns) if columns else [c for c in tab.columns
                                              if c not in (key_col,) and tab[c].dtype.kind in "fi"]
        self.table = tab[cols].astype(float)
        self.roles, self.cols, self.prefix = list(roles), cols, prefix
        self.names_ = [f"{prefix}_{r}_{c}" for r in self.roles for c in cols]

    def transform(self, df):
        out = np.full((len(df), len(self.roles), len(self.cols)), np.nan, np.float32)
        ridx = {r: i for i, r in enumerate(self.roles)}
        for row, comps in enumerate(df["components"]):
            per_role: dict[int, list] = {}
            for role, smi in comps:
                if role in ridx and smi in self.table.index:
                    per_role.setdefault(ridx[role], []).append(self.table.loc[smi].to_numpy())
            for b, vs in per_role.items():
                out[row, b] = np.mean(vs, axis=0)
        return out.reshape(len(df), -1)
