"""Native internal coordinates retain topology, radians and missing-value semantics."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import geometry

DATA = Path(__file__).resolve().parents[2] / "crates/molframe-bench/data/1ubq.cif"


def test_internal_coordinates_and_bat_rebuild_a_real_protein():
    structure = molframe.read(DATA)
    forest = geometry.internal_coordinates(structure)
    assert len(forest.seeds) + len(forest.atoms) == structure.atom_count
    assert forest.atoms
    rebuilt = np.asarray(forest.rebuild(), dtype=np.float32)
    assert np.allclose(rebuilt, structure.coordinates, atol=2e-3)
    shifted = structure.coordinates + np.array([3.0, -2.0, 1.0], dtype=np.float32)
    bat = forest.measure_bat([tuple(row) for row in shifted])
    assert len(bat.coordinates) == len(forest.atoms)
    assert np.allclose(forest.rebuild_bat(bat), shifted, atol=2e-3)
    assert all(0 <= row[1] <= np.pi for row in bat.coordinates)
    assert all(abs(row[2]) <= np.pi for row in bat.coordinates)
    with pytest.raises(molframe.MolframeError):
        forest.measure_bat([None] * structure.atom_count)
    with pytest.raises(molframe.MolframeError):
        forest.measure_bat([])
    with pytest.raises(molframe.MolframeError) as raised:
        geometry.internal_coordinates(structure, model=structure.model_count)
    assert raised.value.code == "MOLFRAME-E6003"


def test_place_atom_reproduces_a_measured_torsion_and_refuses_degeneracy():
    first, second, third = (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)
    result = geometry.place_atom(first, second, third, 2.0, np.pi / 2, np.pi / 3)
    assert result is not None
    arrays = [np.array([point], dtype=np.float32) for point in (first, second, third, result)]
    assert geometry.distances(arrays[2], arrays[3])[0] == pytest.approx(2.0)
    assert geometry.angles(*arrays[1:])[0] == pytest.approx(90.0, abs=1e-5)
    assert geometry.dihedrals(*arrays)[0] == pytest.approx(60.0, abs=1e-5)
    assert geometry.place_atom(first, first, third, 2.0, np.pi / 2, 0.0) is None
