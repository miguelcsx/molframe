"""Fixtures shared by the Python tests."""

import os
import shutil
import tempfile
from collections.abc import Callable
from pathlib import Path

import pytest

import molframe

ROOT = Path(__file__).resolve().parents[2]
CCD = ROOT / "crates" / "molframe-chem" / "data" / "CCD-amino-acids.cif"
DATA = ROOT / "data"
PROFILE = """\
ramachandran: top8000-rama
rotamer: top8000-rotamer
chi: chi-atoms
nucleic: nucleic-sanity
geometry: alphafold-peptide-link
heavy-roles: heavy-atom-roles
max-area: tien-2013-max-area
components: fixture-ccd
score-matrix: blosum62
vdw-radii: bondi
covalent-radii: pyykko-covalent
metals: macromolecular-metals
ionic-radii: shannon-ionic
bond-reach: macromolecular-reach
"""
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


# One carbon in chain A at the origin; assembly 1 is that chain and a copy moved 3 A along x.
DIMER = """\
data_dimer
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 1 A 0 0 0
2 H H1 GLY A 1 1 A 1.0 0 0
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 0 1 0 0 0 0 1 0
2 1 0 0 3 0 1 0 0 0 0 1 0
_pdbx_struct_assembly.id 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '1,2' 'A'
"""


@pytest.fixture
def dimer(tmp_path: Path) -> molframe.Structure:
    path = tmp_path / "dimer.cif"
    path.write_text(DIMER)
    return molframe.read(path)


def pytest_configure(config: pytest.Config) -> None:
    """Point the table cache at the repository's data files before any native call."""
    cache = Path(tempfile.mkdtemp(prefix="molframe-tables-"))
    for entry in DATA.iterdir():
        if entry.is_file() and entry.suffix != ".py":
            shutil.copy(entry, cache / entry.name)
    shutil.copy(CCD, cache / "fixture-ccd.cif")
    (cache / "profile.txt").write_text(PROFILE)
    os.environ["MOLFRAME_DATA_DIR"] = str(cache)
    os.environ.pop("MOLFRAME_PROFILE", None)
    config.stash_cache = cache  # pyright: ignore[reportAttributeAccessIssue]


def pytest_unconfigure(config: pytest.Config) -> None:
    """Remove the temporary table cache."""
    cache = getattr(config, "stash_cache", None)
    if cache is not None:
        shutil.rmtree(cache, ignore_errors=True)
