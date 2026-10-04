"""Fixtures shared by the Python tests."""

from collections.abc import Callable
from pathlib import Path

import pytest

import molframe

ROOT = Path(__file__).resolve().parents[2]
CCD = ROOT / "crates" / "molframe-chem" / "data" / "CCD-amino-acids.cif"
ROLES = (
    ("N", "protein_nitrogen"),
    ("CA", "protein_alpha_carbon"),
    ("C", "protein_carbonyl_carbon"),
    ("O", "protein_carbonyl_oxygen"),
    ("CB", "protein_beta_carbon"),
)


def _with_roles(structure: molframe.Structure) -> molframe.Structure:
    roles = molframe.chemistry.polymer_atom_roles()
    rules = [
        molframe.chemistry.PolymerRoleRule(name, roles[role], component_kind=1)
        for name, role in ROLES
    ]
    return molframe.chemistry.apply_polymer_role_profile(
        structure, CCD, rules, profile_id="backbone", version="wwPDB-2026-10-03"
    ).structure


@pytest.fixture
def with_roles() -> Callable[[molframe.Structure], molframe.Structure]:
    """Return a function that annotates a protein's backbone and beta-carbon atom roles."""
    return _with_roles
