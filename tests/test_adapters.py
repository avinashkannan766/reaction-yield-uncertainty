"""TDC adapter against fabricated pickles in each plausible layout (the real file was not
reachable from the dev sandbox, so tolerance to layout is tested explicitly)."""

import pandas as pd
import pytest

from chemyield.adapters import TDCPickleAdapter

A, B, P = "Brc1ccccc1", "Nc1ccccc1", "c1ccc(Nc2ccccc2)cc1"


def _load(tmp_path, df, **kw):
    df.to_pickle(tmp_path / "x.pkl")
    return TDCPickleAdapter(path=tmp_path / "x.pkl", **kw).load()


def test_reaction_smiles_cell(tmp_path):
    df = pd.DataFrame({"Reaction_ID": [1], "Reactant": [f"{A}.{B}>CCN>{P}"], "Yield": [55.0]})
    out = _load(tmp_path, df)
    assert out["components"][0] == [("reactant", A), ("reactant", B), ("reagent", "CCN")]
    assert out["product"][0] == P and out["reaction_id"][0] == "1"


def test_dotted_string_cell_with_product_column(tmp_path):
    out = _load(tmp_path, pd.DataFrame({"Reactant": [f"{A}.{B}"], "Product": [P], "Yield": [0.5]}))
    assert [s for _, s in out["components"][0]] == [A, B] and out["product"][0] == P


def test_dict_cell_keeps_roles(tmp_path):
    cell = {"aryl_halide": A, "amine": [B], "ligand": "CCN"}
    out = _load(tmp_path, pd.DataFrame({"Reactant": [cell], "Product": [P], "Yield": [10.0]}))
    assert sorted(out["components"][0]) == [("amine", B), ("aryl_halide", A), ("ligand", "CCN")]


def test_list_cell_and_explicit_columns(tmp_path):
    df = pd.DataFrame({"rs": [[A, B]], "pr": [[P]], "yy": [3.0]})
    out = _load(tmp_path, df, reactant_col="rs", product_col="pr", yield_col="yy")
    assert out["product"][0] == P and len(out["components"][0]) == 2


def test_unrecognised_layout_gives_actionable_error(tmp_path):
    with pytest.raises(ValueError, match="columns"):
        _load(tmp_path, pd.DataFrame({"foo": [1], "bar": [2]}))


# ---- the REAL TDC layout, taken from the first row of the actual buchwald-hartwig.pkl ----------
REAL_ROW = {
    "ID": "reactions_1",
    "X": {"product": "Cc1ccc(Nc2ccc(C(F)(F)F)cc2)cc1", "catalyst": "",
          "reactant": "FC(F)(F)c1ccc(Cl)cc1.Cc1ccc(N)cc1.O=S(=O)(O[Pd]1c2ccccc2-c2ccccc2N~1)C(F)(F)F."
                      "CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2CCCCC2)c(C(C)C)c1."
                      "CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C)C.c1ccc(-c2ccno2)cc1"},
    "Y": 0.10657811819999999,
}


def test_real_tdc_layout_end_to_end_with_config_role_rules(tmp_path):
    import yaml

    from chemyield.cleaning import clean_reactions

    pd.DataFrame([REAL_ROW]).to_pickle(tmp_path / "x.pkl")
    raw = TDCPickleAdapter(path=tmp_path / "x.pkl").load()
    assert raw["reaction_id"][0] == "reactions_1" and len(raw["components"][0]) == 6
    rules = yaml.safe_load(open("configs/tdc_buchwald.yaml"))["cleaning"]["role_rules"]
    df, rep = clean_reactions(raw, {"role_rules": rules})
    roles = {r: s for r, s in df["components"][0]}
    assert sorted(roles) == ["additive", "amine", "aryl_halide", "base", "catalyst", "ligand"]
    assert "Pd" in roles["catalyst"] and "P(" in roles["ligand"] and roles["amine"] == "Cc1ccc(N)cc1"
    assert abs(df["yield"][0] - 10.657811819999999) < 1e-6            # fraction -> percent
    assert rep["role_signatures"][0]["roles"] == {r: 1 for r in roles}


def test_uspto_layout_atom_map_roles(tmp_path):
    """Real TDC USPTO rows: 'catalyst' key holds a solvent; unmapped reactant fragments are agents."""
    import pandas as pd
    from chemyield.adapters import TDCPickleAdapter

    rows = [
        {"ID": "reactions_1", "X": {
            "catalyst": "C(Cl)Cl",
            "reactant": "[Cl-].[Al+3].[Cl-].[Cl-].[Cl:5][CH2:6][CH2:7][CH2:8][C:9](Cl)=[O:10]."
                        "[C:12]1([CH:18]([CH3:20])[CH3:19])[CH:17]=[CH:16][CH:15]=[CH:14][CH:13]=1",
            "product": "[Cl:5][CH2:6][CH2:7][CH2:8][C:9]([C:15]1[CH:16]=[CH:17][C:12]([CH:18]([CH3:20])[CH3:19])=[CH:13][CH:14]=1)=[O:10]"},
         "Y": 0.86},
        {"ID": "reactions_2", "X": {"catalyst": "", "reactant": "[Al+3].[Cl-].[Cl-].[Cl-]."
                                   "[Cl:5][CH2:6][CH2:7][CH2:8][C:9](Cl)=[O:10].[C:12]1([CH3:18])[CH:17]=[CH:16][CH:15]=[CH:14][CH:13]=1",
                                   "product": "[Cl:5][CH2:6][CH2:7][CH2:8][C:9]([C:15]1[CH:16]=[CH:17][C:12]([CH3:18])=[CH:13][CH:14]=1)=[O:10]"},
         "Y": 0.95},
    ]
    path = tmp_path / "u.pkl"
    pd.DataFrame(rows).to_pickle(path)
    df = TDCPickleAdapter(path, atom_map_roles=True, unmapped_role="agent",
                          role_renames={"catalyst": "agent"}).load()
    roles0 = [r for r, _ in df["components"][0]]
    assert roles0.count("reactant") == 2                 # the acid chloride and cumene (mapped)
    assert roles0.count("agent") == 5                    # Cl-, Al3+, Cl-, Cl-, plus the DCM solvent
    assert [r for r, _ in df["components"][1]].count("reactant") == 2
    # default behaviour (no flags) is unchanged: 'catalyst' stays a role name
    plain = TDCPickleAdapter(path).load()
    assert "catalyst" in [r for r, _ in plain["components"][0]]
