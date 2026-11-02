import numpy as np
import pandas as pd
import pytest
import torch

from chemyield.features import FeatureSet, GraphFeaturizer
from chemyield.features.graph import ATOM_DIM, BOND_DIM, MolTable
from chemyield.models.mpnn import DMPNNEncoder

SPECS = [{"name": "morgan", "n_bits": 64}, {"name": "rdkit_desc"},
         {"name": "conditions", "numeric_cols": ["temperature"], "onehot_roles": ["ligand"]}]


def test_shapes_and_determinism(clean_synth):
    fs = FeatureSet(SPECS).fit(clean_synth)
    X1, X2 = fs.transform(clean_synth), fs.transform(clean_synth)
    assert X1.shape == (len(clean_synth), len(fs.names)) and np.array_equal(X1, X2, equal_nan=True)


def test_unseen_ligand_is_all_zero_onehot(clean_synth):
    ligs = sorted({s for c in clean_synth["components"] for r, s in c if r == "ligand"})
    train = clean_synth[[all(s != ligs[0] for r, s in c if r == "ligand") for c in clean_synth["components"]]]
    held = clean_synth[[any(s == ligs[0] for r, s in c if r == "ligand") for c in clean_synth["components"]]]
    f = FeatureSet([{"name": "conditions", "onehot_roles": ["ligand"]}]).fit(train)
    assert f.transform(held).sum() == 0 and f.transform(train.head(5)).sum() == 5


def test_fingerprints_depend_on_role_not_just_identity(clean_synth):
    """Same molecule in a different role must land in a different block."""
    fs = FeatureSet([{"name": "morgan", "n_bits": 64, "product": False}]).fit(clean_synth)
    row = clean_synth.head(1).copy()
    swapped = row.copy()
    swapped["components"] = [[("amine" if r == "aryl_halide" else "aryl_halide" if r == "amine" else r, s)
                              for r, s in row["components"].iloc[0]]]
    assert not np.array_equal(fs.transform(row), fs.transform(swapped))


def test_lookup_featurizer(tmp_path, clean_synth):
    solvents = sorted({s for c in clean_synth["components"] for r, s in c if r == "solvent"})
    pd.DataFrame({"smiles": solvents, "dielectric": np.arange(len(solvents), dtype=float)}).to_csv(
        tmp_path / "t.csv", index=False)
    f = FeatureSet([{"name": "lookup", "table": str(tmp_path / "t.csv"), "roles": ["solvent", "base"]}])
    X = f.fit_transform(clean_synth)
    assert X.shape[1] == 2 and not np.isnan(X[:, 0]).any() and np.isnan(X[:, 1]).all()


def test_graph_edges_pair_up_with_reverse():
    t = MolTable()
    g = t.graphs[t.add("CC(=O)O")]
    ei, rev = g["ei"], g["rev"]
    assert g["x"].shape[1] == ATOM_DIM and g["ea"].shape[1] == BOND_DIM
    assert torch.equal(ei[0][rev], ei[1]) and torch.equal(ei[1][rev], ei[0])


def test_dmpnn_is_permutation_invariant_and_handles_bondless_atoms():
    torch.manual_seed(0)
    enc = DMPNNEncoder(hidden=16, depth=3).eval()
    a, b = MolTable(), MolTable()
    a.add("OCC(=O)N"); a.add("[Na+]")
    b.add("NC(=O)CO"); b.add("[Na+]")             # same molecule, different atom order
    ea = enc(*a.collate([0, 1]))
    eb = enc(*b.collate([0, 1]))
    assert torch.allclose(ea, eb, atol=1e-5)
    assert torch.isfinite(ea).all()


def test_graph_featurizer_roles(clean_synth):
    gf = GraphFeaturizer().fit(clean_synth)
    gd = gf.transform(clean_synth.head(10))
    assert gd.n_roles == len(gf.role_names_) == 6 and len(gd) == 10
    graph, p_rxn, p_role, p_mol, n, _ = gd.batch(np.arange(4))
    assert n == 4 and len(p_rxn) == len(p_role) == len(p_mol) and graph[5] >= 1
