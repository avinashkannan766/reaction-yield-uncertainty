import pandas as pd

from chemyield.cleaning import clean_reactions, load_clean, save_clean


def _row(i, comps, prod, y):
    return {"reaction_id": str(i), "components": comps, "product": prod, "yield": y}


def test_cleaning_removes_each_kind_of_dirt():
    good = [("reactant", "Brc1ccccc1"), ("reactant", "Nc1ccccc1")]
    rows = [
        _row(0, good, "c1ccc(Nc2ccccc2)cc1", 50.0),                              # fine
        _row(1, good, "c1ccc(Nc2ccccc2)cc1", 70.0),                              # replicate of 0
        _row(2, [("reactant", "not_a_smiles(")], "CC", 10.0),                    # unparsable
        _row(3, good, "c1ccc(Nc2ccccc2)cc1", float("nan")),                      # missing yield
        _row(4, [("reactant", "Clc1ccccc1"), ("reactant", "CCN")], "CCNc1ccccc1", 250.0),  # impossible
        _row(5, [("reactant", "Clc1ccccc1"), ("reactant", "CCN")], "CCNc1ccccc1", 102.0),  # clip
        _row(6, [("reactant", "CCO")], "CCS", 40.0),                             # S from nowhere
        _row(7, [("reactant", "[CH3:1]Br"), ("reactant", "CCN")], "CCNC", 30.0),  # atom map stripped
    ]
    df, rep = clean_reactions(pd.DataFrame(rows))
    ids = set(df["reaction_id"])
    assert "2" not in ids and "3" not in ids and "4" not in ids
    assert "5" in ids and df.loc[df.reaction_id == "5", "yield"].iloc[0] == 100.0
    # replicates 0 and 1 merged to their mean, std recorded
    r0 = df[df.reaction_id == "0"].iloc[0]
    assert r0["yield"] == 60.0 and r0["n_replicates"] == 2 and r0["yield_std"] == 10.0
    # element mismatch is flagged (not dropped) by default
    assert bool(df.loc[df.reaction_id == "6", "element_mismatch"].iloc[0])
    assert all(":" not in s for comps in df["components"] for _, s in comps)
    assert rep["n_raw"] == 8 and rep["n_clean"] == len(df)


def test_fraction_yields_are_rescaled():
    comps = [("reactant", "Brc1ccccc1"), ("reactant", "Nc1ccccc1")]
    rows = [_row(i, comps + [("base", s)], "c1ccc(Nc2ccccc2)cc1", y)
            for i, (s, y) in enumerate([("CCN", 0.2), ("CCCN", 0.9), ("CCCCN", 0.5)])]
    df, _ = clean_reactions(pd.DataFrame(rows))
    assert df["yield"].max() == 90.0


def test_role_rules_split_an_unlabeled_bag():
    comps = [("reactant", "Brc1ccccc1"), ("reactant", "CC(C)(C)P(C(C)(C)C)C(C)(C)C")]
    cfg = {"role_rules": [{"role": "ligand", "smarts": "[P;X3]"}]}
    df, _ = clean_reactions(pd.DataFrame([_row(0, comps, "c1ccccc1", 10.0)]), cfg)
    roles = dict((s, r) for r, s in df["components"].iloc[0])
    assert roles["CC(C)(C)P(C(C)(C)C)C(C)(C)C"] == "ligand" and roles["Brc1ccccc1"] == "reactant"


def test_roundtrip_save_load(tmp_path, clean_synth):
    save_clean(clean_synth, tmp_path / "c.csv")
    back = load_clean(tmp_path / "c.csv")
    assert back["components"].iloc[0] == clean_synth["components"].iloc[0]
    assert len(back) == len(clean_synth)
