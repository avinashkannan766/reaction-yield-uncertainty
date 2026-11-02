"""RDKit helpers shared by cleaning, featurizers and splits."""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


def mol_from_smiles(smi):
    """Parse SMILES. Returns (mol, sanitized_ok).

    Organometallic catalysts often fail sanitization; we fall back to an
    unsanitized parse so they are kept (and counted) instead of silently dropped.
    """
    if not isinstance(smi, str) or not smi.strip():
        return None, False
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        return mol, True
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        return None, False
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
    except Exception:
        return None, False
    return mol, False


@lru_cache(maxsize=None)
def canonical_smiles(smi: str):
    """Canonical SMILES without atom maps. Returns (smiles | None, sanitized_ok)."""
    mol, ok = mol_from_smiles(smi)
    if mol is None:
        return None, False
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        return Chem.MolToSmiles(mol), ok
    except Exception:
        return None, False


@lru_cache(maxsize=None)
def cached_mol(smi: str):
    """Shared read-only Mol for featurizers. Do not mutate the result."""
    mol, _ = mol_from_smiles(smi)
    return mol


def split_reaction_smiles(rxn: str):
    """'A.B>C>P' -> (reactants, reagents, products) as lists of SMILES."""
    parts = rxn.split(">")
    if len(parts) != 3:
        raise ValueError(f"Not a reaction SMILES: {rxn!r}")
    return tuple([s for s in side.split(".") if s] for side in parts)


def murcko_scaffold(smi: str) -> str:
    mol = cached_mol(smi)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return ""


def generic_scaffold(smi: str) -> str:
    """Murcko scaffold with every atom -> C and every bond -> single: keeps only the ring/chain
    topology, so many chemically different molecules share a group (coarser than murcko_scaffold)."""
    mol = cached_mol(smi)
    if mol is None:
        return ""
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaf))
    except Exception:
        return ""


def element_set(smiles_list) -> set[int]:
    """Atomic numbers present (hydrogen excluded)."""
    out: set[int] = set()
    for smi in smiles_list:
        mol = cached_mol(smi)
        if mol is not None:
            out |= {a.GetAtomicNum() for a in mol.GetAtoms()}
    out.discard(1)
    return out
