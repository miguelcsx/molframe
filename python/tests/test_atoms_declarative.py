"""Declarative atom selection, nearest-atom search and extraction through the Python contract."""

from pathlib import Path

import numpy as np
import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-py" / "tests" / "data"


def _water_and_sulfur() -> molframe.Structure:
    coordinates = np.array(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0], [4.5, 0.0, 0.0]],
        dtype=np.float32,
    )
    return molframe.chemistry.Molecule(["C", "S", "O", "S"], coordinates).to_structure()


def test_selection_from_indices_sorts_deduplicates_and_checks_bounds():
    structure = _water_and_sulfur()
    selection = molframe.Selection.from_indices(structure, np.array([3, 1, 1], dtype=np.uint32))
    assert selection.indices.tolist() == [1, 3]
    assert len(selection) == 2
    with pytest.raises(ValueError, match="outside the structure"):
        molframe.Selection.from_indices(structure, np.array([4], dtype=np.uint32))


def test_selection_from_mask_selects_true_entries_and_checks_length():
    structure = _water_and_sulfur()
    mask = np.array([False, True, False, True])
    assert molframe.Selection.from_mask(structure, mask).indices.tolist() == [1, 3]
    with pytest.raises(ValueError, match="mask length"):
        molframe.Selection.from_mask(structure, np.array([True, False]))


def test_element_query_selects_atoms_of_one_symbol_case_insensitively():
    structure = _water_and_sulfur()
    sulfur = structure.select(molframe.sel.element("s"))
    assert sulfur.indices.tolist() == [1, 3]


def test_nearest_matches_a_brute_force_search_on_a_real_structure():
    structure = molframe.read(DATA / "basic.cif")
    rng = np.random.default_rng(7)
    points = rng.uniform(-2.0, 2.0, size=(40, 3)).astype(np.float32)
    points += structure.coordinates.mean(axis=0)
    hit = molframe.nearest(structure, points)
    coordinates = structure.coordinates.astype(np.float64)
    distances = np.linalg.norm(points[:, None, :] - coordinates[None, :, :], axis=2)
    assert hit.found.all()
    assert hit.index.tolist() == distances.argmin(axis=1).tolist()
    assert hit.distance == pytest.approx(distances.min(axis=1), abs=1e-4)


def test_nearest_restricted_to_a_selection_still_reports_structure_indices():
    structure = _water_and_sulfur()
    sulfur = molframe.Selection.from_indices(structure, np.array([1, 3], dtype=np.uint32))
    hit = molframe.nearest(structure, [[2.9, 0.0, 0.0]], among=sulfur)
    assert hit.index.tolist() == [1]
    # Atom 1 sits at x=1.5 and atom 3 at x=4.5, so the point at x=2.9 is 1.4 Å from atom 1.
    assert hit.distance[0] == pytest.approx(1.4, abs=1e-5)


def test_nearest_reports_points_beyond_max_distance_as_not_found():
    structure = _water_and_sulfur()
    hit = molframe.nearest(structure, [[0.0, 0.0, 0.5], [0.0, 50.0, 0.0]], max_distance=2.0)
    assert hit.found.tolist() == [True, False]
    assert hit.index[1] == np.iinfo(np.uint32).max
    assert np.isinf(hit.distance[1])


def test_extract_keeps_exactly_the_selected_atoms():
    structure = _water_and_sulfur()
    sulfur = molframe.Selection.from_indices(structure, np.array([1, 3], dtype=np.uint32))
    extracted = molframe.extract(structure, sulfur)
    assert extracted.atom_count == 2
    assert [atom.element for atom in extracted.atoms] == ["S", "S"]
    assert structure.atom_count == 4
