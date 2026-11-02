"""SYNTHETIC Buchwald-Hartwig-like data for smoke tests and pipeline development.

Products are real (built with an RDKit reaction), but YIELDS ARE INVENTED from a hidden additive +
interaction rule plus noise. Results on this data say nothing about real chemistry - it only proves
the pipeline runs and that models can recover a known signal.

    python -m chemyield.synthetic data/raw/synthetic/synthetic_bh.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ARYL_HALIDES = ["Brc1ccccc1", "Brc1ccc(C)cc1", "Brc1ccc(OC)cc1", "Brc1ccc(C(F)(F)F)cc1",
                "Clc1ccccc1", "Clc1ccc(C)cc1", "Clc1ccc(OC)cc1", "Clc1ccc(C(F)(F)F)cc1",
                "Brc1cccnc1", "Clc1cccnc1", "Ic1ccccc1", "Ic1ccc(C)cc1", "Brc1ccccc1C",
                "Clc1ccccc1C"]
AMINES = ["Nc1ccccc1", "Cc1ccc(N)cc1", "COc1ccc(N)cc1", "C1CCNCC1", "C1COCCN1", "CCCCN",
          "CCNCC", "NCc1ccccc1"]
LIGANDS = ["c1ccc(cc1)P(c1ccccc1)c1ccccc1", "CC(C)(C)P(C(C)(C)C)C(C)(C)C",
           "CC(C)c1cc(C(C)C)c(c(c1)C(C)C)-c1ccccc1P(C1CCCCC1)C1CCCCC1",
           "C1CCC(CC1)P(C1CCCCC1)C1CCCCC1"]
BASES = ["C1CCC2=NCCCN2CC1", "CC(C)(C)O[Na]", "CCN(CC)CC"]
SOLVENTS = ["Cc1ccccc1", "C1CCOC1", "C1COCCO1", "CN(C)C=O"]

_RXN = AllChem.ReactionFromSmarts("[c:1][Cl,Br,I].[N;H2,H1;!$(NC=O):2]>>[c:1][N:2]")


def _product(halide: str, amine: str) -> str:
    ps = _RXN.RunReactants((Chem.MolFromSmiles(halide), Chem.MolFromSmiles(amine)))
    mol = ps[0][0]
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def make_dataframe(seed: int = 0, replicate_frac: float = 0.15) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    eff = {name: rng.normal(0, 1, len(lst)) for name, lst in
           [("halide", ARYL_HALIDES), ("amine", AMINES), ("ligand", LIGANDS), ("base", BASES),
            ("solvent", SOLVENTS)]}
    reactivity = {"I": 1.0, "Br": 0.6, "Cl": 0.0}
    bulky = np.array([0.0, 1.0, 1.0, 0.6])        # ligand 'bulk' helps hard (Cl) substrates
    rows = []
    for i, h in enumerate(ARYL_HALIDES):
        x = "I" if "I" in h[:2] else ("Br" if "Br" in h else "Cl")
        for j, a in enumerate(AMINES):
            prod = _product(h, a)
            for k, l in enumerate(LIGANDS):
                for m, b in enumerate(BASES):
                    for n, s in enumerate(SOLVENTS):
                        temp = int(rng.choice([60, 80, 100]))
                        base_score = (35 + 12 * eff["halide"][i] + 10 * eff["amine"][j]
                                      + 8 * eff["ligand"][k] + 6 * eff["base"][m]
                                      + 4 * eff["solvent"][n] + 25 * reactivity[x]
                                      + 22 * bulky[k] * (1 - reactivity[x]) + 0.15 * (temp - 80))
                        for _ in range(2 if rng.random() < replicate_frac else 1):
                            y = float(np.clip(base_score + rng.normal(0, 4), 0, 100))
                            rows.append({"aryl_halide": h, "amine": a, "ligand": l, "base": b,
                                         "solvent": s, "temperature": temp, "product": prod,
                                         "yield": round(y, 2)})
    df = pd.DataFrame(rows)
    df.insert(0, "id", [f"syn{i}" for i in range(len(df))])
    return df


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "data/raw/synthetic/synthetic_bh.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    d = make_dataframe()
    d.to_csv(out, index=False)
    print(f"wrote {len(d)} SYNTHETIC reactions to {out}")
