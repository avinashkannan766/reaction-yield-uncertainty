import pytest

from chemyield.adapters import CSVAdapter
from chemyield.cleaning import clean_reactions
from chemyield.synthetic import make_dataframe


@pytest.fixture(scope="session")
def raw_synth(tmp_path_factory):
    """Small SYNTHETIC reaction table passed through the real CSV adapter."""
    path = tmp_path_factory.mktemp("d") / "s.csv"
    make_dataframe().sample(1200, random_state=0).to_csv(path, index=False)
    return CSVAdapter(
        path=path, id_col="id", product_col="product", keep_cols=["temperature"],
        role_cols={"aryl_halide": "aryl_halide", "amine": "amine", "ligand": "ligand",
                   "base": "base", "solvent": "solvent"}).load()


@pytest.fixture(scope="session")
def clean_synth(raw_synth):
    df, _ = clean_reactions(raw_synth)
    return df
