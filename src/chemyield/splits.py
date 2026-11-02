"""Train / val / cal / test splits.

random : i.i.d. rows (optimistic - use only as the 'easy' reference)
group  : the TEST set is made of whole groups never seen in training. `group_by` says what a
         group is:
           scaffold         Murcko scaffold of the product
           cluster:<k>      k-means clusters of product fingerprints (for diverse data: every
                            scaffold / product is nearly unique, clusters are the unit of 'unseen')
           generic_scaffold Murcko scaffold with all atoms->C, bonds->single (coarser: topology only)
           product          identical product
           role:<name>      the set of components in that role (e.g. role:ligand -> unseen ligand)
           column:<name>    any column of the cleaned table
loo    : leave-one-group-out. Fold = seed % n_groups, so seeds [0..G-1] visit every group once
         (e.g. group_by role:ligand with 4 ligands -> seeds [0,1,2,3]). The whole group is the test set.
The remaining rows are split at random into train / val (early stopping) / cal (conformal) -
unless `val_by_group: true`, in which case val and cal are ALSO whole groups unseen in training
(needs enough groups; otherwise the split is skipped with a "too few usable groups" message).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit.Chem import rdFingerprintGenerator

from .chem import cached_mol, generic_scaffold, murcko_scaffold


_CLUSTER_CACHE: dict = {}


def _product_clusters(products: list[str], k: int) -> list[str]:
    """k-means (k clusters, fixed random_state => same groups for every seed) on 512-bit Morgan
    fingerprints of each product's largest fragment. For DIVERSE datasets, where every scaffold /
    product is nearly unique (median group size 1), whole clusters are the unit of 'unseen chemistry'."""
    key = (hash(tuple(products)), k)
    if key in _CLUSTER_CACHE:
        return _CLUSTER_CACHE[key]
    from sklearn.cluster import MiniBatchKMeans
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=512)
    uniq = sorted(set(products))
    X = np.zeros((len(uniq), 512), np.float32)
    ok = np.zeros(len(uniq), bool)
    for i, p in enumerate(uniq):
        mol = cached_mol(max(p.split("."), key=len))
        if mol is not None:
            try:
                X[i] = gen.GetFingerprintAsNumPy(mol)
                ok[i] = True
            except Exception:
                pass
    label = {p: "bad" for p in uniq}
    if ok.sum() >= 2:
        km = MiniBatchKMeans(n_clusters=int(min(k, ok.sum())), n_init=3, random_state=0,
                             batch_size=4096).fit(X[ok])
        for p, lab in zip(np.array(uniq)[ok], km.labels_):
            label[str(p)] = f"c{int(lab)}"
    out = [label[p] for p in products]
    _CLUSTER_CACHE[key] = out
    return out


def group_keys(df: pd.DataFrame, spec: str) -> list:
    if spec == "product":
        return df["product"].tolist()
    if spec == "scaffold":
        cache: dict[str, str] = {}
        out = []
        for p in df["product"]:
            if p not in cache:
                biggest = max(p.split("."), key=len)
                cache[p] = murcko_scaffold(biggest) or p
            out.append(cache[p])
        return out
    if spec.startswith("cluster:"):
        return _product_clusters(df["product"].tolist(), int(spec.split(":", 1)[1]))
    if spec == "generic_scaffold":
        cache2: dict[str, str] = {}
        out = []
        for p in df["product"]:
            if p not in cache2:
                biggest = max(p.split("."), key=len)
                cache2[p] = generic_scaffold(biggest) or "acyclic"
            out.append(cache2[p])
        return out
    if spec.startswith("role:"):
        role = spec.split(":", 1)[1]
        return [tuple(sorted(s for r, s in c if r == role)) for c in df["components"]]
    if spec.startswith("column:"):
        return df[spec.split(":", 1)[1]].tolist()
    raise ValueError(f"Unknown group_by '{spec}'")


def _groups(df, group_by):
    groups: dict = {}
    for i, k in enumerate(group_keys(df, group_by)):
        groups.setdefault(k, []).append(i)
    return groups


def held_out_label(df, group_by: str, seed: int) -> str:
    """Name of the group held out by split kind 'loo' for this seed (folds = seed % n_groups)."""
    groups = _groups(df, group_by)
    k = sorted(groups, key=str)[seed % len(groups)]
    return ".".join(k) if isinstance(k, tuple) else str(k)


def _take_groups(pool, groups, target, rng):
    """Take whole groups in random order until >= target rows. Groups far larger than the target
    are skipped so one big group cannot distort the requested fractions."""
    taken, rows = [], []
    for j in rng.permutation(len(pool)):
        if len(rows) >= target:
            break
        k = pool[j]
        if len(groups[k]) > 2 * max(target, 1):
            continue
        taken.append(k)
        rows += groups[k]
    return taken, rows


def _group_val_split(groups, test_idx, n, n_val, n_cal, rng, group_by):
    """Validation (early stopping) and calibration (conformal) sets that are ALSO whole groups the
    model never trains on - so they mimic 'a new group at test time' instead of 'more rows of a
    group already seen'. Random val/cal rows make early stopping favour interpolation and make the
    conformal intervals too narrow for new chemistry."""
    test_set = set(test_idx.tolist())
    rest = {k: [i for i in v if i not in test_set] for k, v in groups.items()}
    rest = {k: v for k, v in rest.items() if v}
    pool = sorted(rest, key=str)
    val_keys, val_rows = _take_groups(pool, rest, n_val, rng)
    pool = [k for k in pool if k not in set(val_keys)]
    cal_keys, cal_rows = _take_groups(pool, rest, n_cal, rng)
    held = set(val_keys) | set(cal_keys)
    train_rows = [i for k, v in rest.items() if k not in held for i in v]
    if len(val_rows) < 0.5 * n_val or len(cal_rows) < 0.5 * n_cal or len(train_rows) < 0.3 * n:
        raise ValueError(
            f"group_by='{group_by}' has too few usable groups ({len(groups)}) to hold out whole "
            "groups for validation and calibration as well as for test (val_by_group=true).")
    return {"val": np.sort(val_rows), "cal": np.sort(cal_rows),
            "train": np.sort(train_rows), "test": np.sort(test_idx)}


def make_split(df: pd.DataFrame, kind: str = "random", seed: int = 0, test_size: float = 0.2,
               val_size: float = 0.1, cal_size: float = 0.1, group_by: str | None = None,
               val_by_group: bool = False):
    n = len(df)
    rng = np.random.default_rng(seed)
    if kind == "random":
        perm = rng.permutation(n)
        test_idx = perm[: int(round(n * test_size))]
        rest = perm[len(test_idx):]
    elif kind == "group":
        if not group_by:
            raise ValueError("split.kind=group needs split.group_by")
        groups: dict = {}
        for i, k in enumerate(group_keys(df, group_by)):
            groups.setdefault(k, []).append(i)
        target = test_size * n
        # very large groups stay in train so one group cannot swallow the whole test set
        candidates = [k for k, v in groups.items() if len(v) <= target]
        order = rng.permutation(len(candidates))
        test_list: list[int] = []
        for j in order:
            if len(test_list) >= target:
                break
            test_list += groups[candidates[j]]
        if len(test_list) < 0.5 * target:
            raise ValueError(
                f"group_by='{group_by}' has too few usable groups ({len(groups)}) to build a "
                f"test set of {test_size:.0%}; pick a finer grouping.")
        test_idx = np.array(sorted(test_list))
        mask = np.ones(n, bool)
        mask[test_idx] = False
        rest = rng.permutation(np.flatnonzero(mask))
    elif kind == "loo":
        if not group_by:
            raise ValueError("split.kind=loo needs split.group_by")
        groups = _groups(df, group_by)
        if len(groups) < 2:
            raise ValueError(f"group_by='{group_by}' has too few usable groups ({len(groups)}) "
                             "for leave-one-out")
        held = sorted(groups, key=str)[seed % len(groups)]
        test_idx = np.array(sorted(groups[held]))
        mask = np.ones(n, bool)
        mask[test_idx] = False
        rest = rng.permutation(np.flatnonzero(mask))
    else:
        raise ValueError(f"Unknown split kind '{kind}'")
    n_val, n_cal = int(round(n * val_size)), int(round(n * cal_size))
    if val_by_group:
        if kind not in ("group", "loo"):
            raise ValueError("val_by_group needs split.kind 'group' or 'loo'")
        return _group_val_split(groups, np.sort(test_idx), n, n_val, n_cal, rng, group_by)
    return {
        "val": np.sort(rest[:n_val]),
        "cal": np.sort(rest[n_val:n_val + n_cal]),
        "train": np.sort(rest[n_val + n_cal:]),
        "test": np.sort(test_idx),
    }
