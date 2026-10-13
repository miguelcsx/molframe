"""Heavy-atom hydrogen bonds work on structures with no hydrogens and no dictionary."""

import math

import numpy as np
import pytest

import molframe

PDB = b"""ATOM      1  OG  SER A   1       0.000   0.000   0.000  1.00 20.00           O
ATOM      2  O   GLY A   2       2.800   0.000   0.000  1.00 20.00           O
ATOM      3  C   GLY A   2       3.400   1.070   0.000  1.00 20.00           C
END
"""


def test_a_serine_hydroxyl_donates_to_a_carbonyl_oxygen() -> None:
    structure = molframe.read(PDB, format="pdb")
    table = molframe.analysis.heavy_atom_hydrogen_bonds(structure)
    assert len(table) == 1
    assert (int(table["donor"][0]), int(table["acceptor"][0])) == (0, 1)
    assert table["distance"][0] == pytest.approx(2.8, abs=1e-5)
    expected = math.degrees(math.acos(-0.6 / math.hypot(0.6, 1.07)))
    assert table["acceptor_angle"][0] == pytest.approx(expected, abs=1e-3)
    assert set(table.names) == {"donor", "acceptor", "distance", "acceptor_angle"}


def test_the_distance_window_and_angle_are_explicit_options() -> None:
    structure = molframe.read(PDB, format="pdb")
    assert len(molframe.analysis.heavy_atom_hydrogen_bonds(structure, max_distance=2.5)) == 0
    assert len(molframe.analysis.heavy_atom_hydrogen_bonds(structure, min_distance=3.0)) == 0
    strict = molframe.analysis.heavy_atom_hydrogen_bonds(structure, min_acceptor_angle=150.0)
    assert len(strict) == 0
    loose = molframe.analysis.heavy_atom_hydrogen_bonds(structure, min_acceptor_angle=None)
    assert len(loose) == 1
    assert np.isfinite(loose["acceptor_angle"][0])


def test_an_invalid_window_is_refused() -> None:
    structure = molframe.read(PDB, format="pdb")
    with pytest.raises(ValueError, match="invalid hydrogen-bond geometry"):
        molframe.analysis.heavy_atom_hydrogen_bonds(structure, min_distance=5.0)


ROLES = """name: lab-roles
version: 4
source: bench notes
licence: CC0-1.0
amino: SER GLY
@amino O acceptor
SER OG donor
"""


def test_a_supplied_role_table_replaces_the_builtin_typing(tmp_path) -> None:
    structure = molframe.read(PDB, format="pdb")
    default = molframe.analysis.heavy_atom_hydrogen_bonds(structure)
    assert "hbplus-baker-hubbard-heavy-roles" in default.applied["roles"]
    path = tmp_path / "roles.txt"
    path.write_text(ROLES)
    for table in (
        molframe.analysis.heavy_atom_hydrogen_bonds(structure, role_table=ROLES),
        molframe.analysis.heavy_atom_hydrogen_bonds(structure, role_file=path),
    ):
        assert (int(table["donor"][0]), int(table["acceptor"][0])) == (0, 1)
        assert "lab-roles 4 (source: bench notes; licence: CC0-1.0)" in table.applied["roles"]
    acceptors_only = ROLES.replace("SER OG donor", "SER OG acceptor")
    assert (
        len(molframe.analysis.heavy_atom_hydrogen_bonds(structure, role_table=acceptors_only)) == 0
    )


def test_a_malformed_role_table_or_both_arguments_are_refused() -> None:
    structure = molframe.read(PDB, format="pdb")
    with pytest.raises(molframe.MolframeValueError, match="line"):
        molframe.analysis.heavy_atom_hydrogen_bonds(
            structure, role_table="name: lab\nSER OG fish\n"
        )
    with pytest.raises(molframe.MolframeValueError, match="not both"):
        molframe.analysis.heavy_atom_hydrogen_bonds(structure, role_table=ROLES, role_file="x")
