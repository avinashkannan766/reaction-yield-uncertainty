"""Dataset adapters: the ONLY place that knows a dataset's raw format.

Every adapter returns a DataFrame in the standard schema:

    reaction_id : str
    components  : list[(role, smiles)]   reactants / reagents / ligands / bases ...
    product     : str                    '.'-joined product SMILES
    yield       : float                  as in the source (cleaning fixes the scale)
    <extras>    : any extra columns the adapter was told to keep (e.g. temperature)

To support a new dataset: write a small Adapter (or just use the generic 'csv'
adapter via config) -- nothing downstream changes.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .chem import split_reaction_smiles
from .registry import Registry

ADAPTERS = Registry("adapter")


def _is_missing(v) -> bool:
    return v is None or (isinstance(v, float) and np.isnan(v))


def _flatten(value, role: str) -> list[tuple[str, str]]:
    """Turn str / list / dict cell values into [(role, smiles), ...].

    dict keys become role names, so role-aware sources keep their roles.
    """
    if _is_missing(value):
        return []
    if isinstance(value, str):
        return [(role, s) for s in value.split(".") if s]
    if isinstance(value, dict):
        out = []
        for k, v in value.items():
            out += _flatten(v, str(k))
        return out
    if isinstance(value, (list, tuple, set, np.ndarray)):
        out = []
        for v in value:
            out += _flatten(v, role)
        return out
    raise TypeError(f"Cannot interpret cell of type {type(value)}: {value!r}")


def _from_reaction_smiles(rxn: str):
    reactants, reagents, products = split_reaction_smiles(rxn)
    comps = [("reactant", s) for s in reactants] + [("reagent", s) for s in reagents]
    return comps, ".".join(products)


class Adapter:
    def load(self) -> pd.DataFrame:  # pragma: no cover - interface
        raise NotImplementedError


@ADAPTERS.register("csv")
class CSVAdapter(Adapter):
    """Generic tabular source. Pick ONE of three ways to describe the reaction:

    * reaction_smiles_col                   : 'A.B>reagents>P' strings
    * role_cols {role: column}  + product_col: HTE-style, one column per component
    * reactants_col [+ reagents_col] + product_col
    """

    def __init__(self, path, yield_col="yield", id_col=None, reaction_smiles_col=None,
                 role_cols=None, reactants_col=None, reagents_col=None, product_col=None,
                 keep_cols=(), sep=","):
        self.path = Path(path)
        self.yield_col = yield_col
        self.id_col = id_col
        self.reaction_smiles_col = reaction_smiles_col
        self.role_cols = role_cols or {}
        self.reactants_col = reactants_col
        self.reagents_col = reagents_col
        self.product_col = product_col
        self.keep_cols = list(keep_cols)
        self.sep = sep

    def load(self) -> pd.DataFrame:
        raw = pd.read_csv(self.path, sep=self.sep)
        comps, products = [], []
        for _, row in raw.iterrows():
            if self.reaction_smiles_col:
                c, p = _from_reaction_smiles(row[self.reaction_smiles_col])
            else:
                c = []
                if self.role_cols:
                    for role, col in self.role_cols.items():
                        c += _flatten(row[col], role)
                else:
                    c += _flatten(row[self.reactants_col], "reactant")
                    if self.reagents_col:
                        c += _flatten(row[self.reagents_col], "reagent")
                p = row[self.product_col] if self.product_col else ""
            comps.append(c)
            products.append(p if isinstance(p, str) else "")
        out = pd.DataFrame({
            "reaction_id": (raw[self.id_col].astype(str) if self.id_col
                            else pd.Series(range(len(raw))).astype(str)),
            "components": comps,
            "product": products,
            "yield": raw[self.yield_col].to_numpy(),
        })
        for col in self.keep_cols:
            out[col] = raw[col].to_numpy()
        return out


@ADAPTERS.register("tdc_pickle")
class TDCPickleAdapter(Adapter):
    """TDC 'Yields' pickle (e.g. buchwald-hartwig.pkl from the TDC dataverse).

    We have not been able to hard-code the exact column layout, so this adapter
    is deliberately tolerant: the reactant cell may be a reaction string, a
    '.'-joined string, a list, or a dict {role: smiles}. If the layout is not
    recognised it raises an error that prints the columns and a sample row, so
    the fix is a one-line config change (reactant_col / product_col / yield_col).

    NOTE: .pkl files are Python pickles and can execute code when loaded. Only
    load the official TDC file you downloaded yourself.
    """

    YIELD_NAMES = ("Yield", "yield", "Y")
    REACTANT_NAMES = ("Reactant", "Reactants", "reactants", "Reaction", "rxn")
    PRODUCT_NAMES = ("Product", "Products", "product")
    ID_NAMES = ("Reaction_ID", "reaction_id", "ID", "id")
    ROLE_EXTRAS = {"Catalyst": "catalyst", "Reagent": "reagent", "Reagents": "reagent",
                   "Solvent": "solvent"}

    _MAPPED = re.compile(r":\d+\]")      # atom-map number, e.g. [Cl:5]

    def __init__(self, path, yield_col=None, reactant_col=None, product_col=None,
                 id_col=None, extra_role_cols=None, keep_cols=(), atom_map_roles=False,
                 unmapped_role="agent", role_renames=None):
        """atom_map_roles: for atom-mapped sources (TDC USPTO yields). Inside the dot-joined
        'reactant' string, a fragment WITH mapped atoms is a real reactant (its atoms end up in the
        product); a fragment without (AlCl3, bases, catalysts ...) becomes `unmapped_role`.
        role_renames: {source_key: role}, e.g. {catalyst: agent} - in USPTO the 'catalyst' key holds
        solvents / reagents / catalysts alike, so it is renamed rather than trusted."""
        self.atom_map_roles = atom_map_roles
        self.unmapped_role = unmapped_role
        self.role_renames = role_renames or {}
        self.path = Path(path)
        self.yield_col, self.reactant_col, self.product_col = yield_col, reactant_col, product_col
        self.id_col = id_col
        self.extra_role_cols = extra_role_cols or {}
        self.keep_cols = list(keep_cols)

    @staticmethod
    def _is_x_dict_layout(raw) -> bool:
        """Real TDC layout: columns [ID, X, Y]; X = {'product':..., 'catalyst':..., 'reactant':...}."""
        return ("X" in raw.columns and "Y" in raw.columns and len(raw) > 0
                and isinstance(raw["X"].iloc[0], dict))

    def _load_x_dict(self, raw) -> pd.DataFrame:
        comps, products = [], []
        for x in raw["X"]:
            c = []
            for key, value in x.items():
                if key == "product":
                    continue
                if self.atom_map_roles and key == "reactant" and isinstance(value, str):
                    c += [("reactant" if self._MAPPED.search(f) else self.unmapped_role, f)
                          for f in value.split(".") if f]
                    continue
                c += _flatten(value, self.role_renames.get(key, key))  # dict keys become roles
            p = x.get("product", "")
            if isinstance(p, (list, tuple)):
                p = ".".join(p)
            comps.append(c)
            products.append(p if isinstance(p, str) else "")
        out = pd.DataFrame({
            "reaction_id": (raw["ID"].astype(str) if "ID" in raw.columns
                            else pd.Series(range(len(raw))).astype(str)),
            "components": comps,
            "product": products,
            "yield": raw["Y"].to_numpy(),
        })
        for col in self.keep_cols:
            out[col] = raw[col].to_numpy()
        return out

    @staticmethod
    def _pick(df, explicit, candidates):
        if explicit:
            return explicit
        return next((c for c in candidates if c in df.columns), None)

    def load(self) -> pd.DataFrame:
        raw = pd.read_pickle(self.path)
        if not isinstance(raw, pd.DataFrame):
            raise TypeError(f"Expected a DataFrame in {self.path}, got {type(raw)}")
        if self._is_x_dict_layout(raw):
            return self._load_x_dict(raw)
        ycol = self._pick(raw, self.yield_col, self.YIELD_NAMES)
        rcol = self._pick(raw, self.reactant_col, self.REACTANT_NAMES)
        pcol = self._pick(raw, self.product_col, self.PRODUCT_NAMES)
        icol = self._pick(raw, self.id_col, self.ID_NAMES)
        if ycol is None or rcol is None:
            raise ValueError(
                "Could not recognise the TDC layout.\n"
                f"  columns : {list(raw.columns)}\n"
                f"  row 0   : {raw.iloc[0].to_dict() if len(raw) else '<empty>'}\n"
                "Set dataset.yield_col / reactant_col / product_col in the config.")
        role_cols = {**{c: r for c, r in self.ROLE_EXTRAS.items() if c in raw.columns},
                     **self.extra_role_cols}
        comps, products = [], []
        for _, row in raw.iterrows():
            cell = row[rcol]
            if isinstance(cell, str) and ">" in cell:
                c, p = _from_reaction_smiles(cell)
            else:
                c = _flatten(cell, "reactant")
                p = row[pcol] if pcol else ""
                if isinstance(p, (list, tuple)):
                    p = ".".join(p)
            for col, role in role_cols.items():
                c += _flatten(row[col], role)
            comps.append(c)
            products.append(p if isinstance(p, str) else "")
        out = pd.DataFrame({
            "reaction_id": (raw[icol].astype(str) if icol
                            else pd.Series(range(len(raw))).astype(str)),
            "components": comps,
            "product": products,
            "yield": raw[ycol].to_numpy(),
        })
        for col in self.keep_cols:
            out[col] = raw[col].to_numpy()
        return out


def build_adapter(dataset_cfg: dict) -> Adapter:
    """dataset: {name: <label>, adapter: <registered name>, params: {...}}"""
    return ADAPTERS.build({"name": dataset_cfg["adapter"], **dataset_cfg.get("params", {})})
