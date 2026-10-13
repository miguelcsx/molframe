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
