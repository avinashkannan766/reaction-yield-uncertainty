"""Molecular graphs for the D-MPNN. Mapping-free: each component is its own graph.

(No atom-mapping is needed, so this works on any dataset. A CGR variant can be added later as
another featurizer for atom-mapped data such as USPTO.)
"""

from __future__ import annotations

import numpy as np
import torch
from rdkit import Chem

from ..chem import cached_mol
from .base import roles_in

_ATOMS = [6, 7, 8, 9, 15, 16, 17, 35, 53, 5, 14, 11, 19, 3, 12, 20, 26, 28, 29, 30, 46, 45, 44, 78]
_HYB = [Chem.HybridizationType.SP, Chem.HybridizationType.SP2, Chem.HybridizationType.SP3,
        Chem.HybridizationType.SP3D, Chem.HybridizationType.SP3D2]
_BOND = [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC]
_STEREO = [Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE,
           Chem.BondStereo.STEREOANY]


def _onehot(value, choices, with_other=True):
    v = [float(value == c) for c in choices]
    if with_other:
        v.append(float(value not in choices))
    return v


def atom_features(a) -> list[float]:
    return (_onehot(a.GetAtomicNum(), _ATOMS) + _onehot(a.GetTotalDegree(), [0, 1, 2, 3, 4, 5])
            + _onehot(a.GetFormalCharge(), [-2, -1, 0, 1, 2]) + _onehot(a.GetTotalNumHs(), [0, 1, 2, 3, 4])
            + _onehot(a.GetHybridization(), _HYB) + [float(a.GetIsAromatic()), a.GetMass() / 100.0]
            + _onehot(int(a.GetChiralTag()), [0, 1, 2], with_other=True))


def bond_features(b) -> list[float]:
    return (_onehot(b.GetBondType(), _BOND, with_other=False)
            + [float(b.GetIsConjugated()), float(b.IsInRing())]
            + _onehot(b.GetStereo(), _STEREO, with_other=True))


_probe = Chem.MolFromSmiles("CC=O")
ATOM_DIM = len(atom_features(_probe.GetAtomWithIdx(0)))
BOND_DIM = len(bond_features(_probe.GetBondWithIdx(0)))


class MolTable:
    """Growing store of molecular graphs, indexed by id. Shared by train/val/test so a molecule
    is converted once; unseen test molecules are simply appended."""

    def __init__(self):
        self.index: dict[str, int] = {}
        self.graphs: list[dict] = []

    def add(self, smi: str) -> int:
        if smi in self.index:
            return self.index[smi]
        mol = cached_mol(smi)
        if mol is None:
            raise ValueError(f"Cannot build graph for {smi}")
        x = torch.tensor([atom_features(a) for a in mol.GetAtoms()], dtype=torch.float32)
        src, dst, ea = [], [], []
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            f = bond_features(b)
            src += [i, j]
            dst += [j, i]
            ea += [f, f]
        E = len(src)
        # edges are appended in (i->j, j->i) pairs, so the reverse of edge e is e ^ 1
        rev = torch.tensor([e ^ 1 for e in range(E)], dtype=torch.long)
        self.graphs.append({
            "x": x,
            "ei": torch.tensor([src, dst], dtype=torch.long).reshape(2, E),
            "ea": torch.tensor(ea, dtype=torch.float32).reshape(E, BOND_DIM),
            "rev": rev,
        })
        self.index[smi] = len(self.graphs) - 1
        return self.index[smi]

    def collate(self, ids: list[int]):
        xs, eis, eas, revs, batch = [], [], [], [], []
        a_off = e_off = 0
        for k, i in enumerate(ids):
            g = self.graphs[i]
            n, e = g["x"].shape[0], g["ei"].shape[1]
            xs.append(g["x"])
            eis.append(g["ei"] + a_off)
            eas.append(g["ea"])
            revs.append(g["rev"] + e_off)
            batch.append(torch.full((n,), k, dtype=torch.long))
            a_off += n
            e_off += e
        return (torch.cat(xs), torch.cat(eis, 1), torch.cat(eas), torch.cat(revs),
                torch.cat(batch), len(ids), torch.tensor(ids, dtype=torch.long))


class GraphData:
    """A set of reactions as (role_id, mol_id) pairs + optional tabular side features."""

    def __init__(self, table: MolTable, items: list[list[tuple[int, int]]], extra=None,
                 n_roles: int = 0):
        self.table, self.items, self.extra, self.n_roles = table, items, extra, n_roles

    def __len__(self):
        return len(self.items)

    def subset_extra(self, extra):
        return GraphData(self.table, self.items, extra, self.n_roles)

    def batch(self, idxs):
        mol_ids = sorted({m for i in idxs for _, m in self.items[i]})
        pos = {m: j for j, m in enumerate(mol_ids)}
        p_rxn, p_role, p_mol = [], [], []
        for b, i in enumerate(idxs):
            for r, m in self.items[i]:
                p_rxn.append(b)
                p_role.append(r)
                p_mol.append(pos[m])
        extra = None if self.extra is None else torch.from_numpy(self.extra[idxs]).float()
        return (self.table.collate(mol_ids), torch.tensor(p_rxn), torch.tensor(p_role),
                torch.tensor(p_mol), len(idxs), extra)


class GraphFeaturizer:
    def __init__(self, roles="auto", product: bool = True):
        self.roles_cfg, self.product = roles, product
        self.table = MolTable()

    def fit(self, df):
        roles = roles_in(df) if self.roles_cfg == "auto" else list(self.roles_cfg)
        self.role_names_ = roles + (["product"] if self.product else [])
        self.role_idx_ = {r: i for i, r in enumerate(self.role_names_)}
        return self

    def transform(self, df, extra=None) -> GraphData:
        items = []
        for comps, prod in zip(df["components"], df["product"]):
            row = [(self.role_idx_[r], self.table.add(s)) for r, s in comps if r in self.role_idx_]
            if self.product:
                row += [(self.role_idx_["product"], self.table.add(s)) for s in prod.split(".")]
            items.append(row)
        return GraphData(self.table, items, extra, len(self.role_names_))
