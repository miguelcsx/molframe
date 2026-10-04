"""CE alignment, CAD scores, contact-map overlap and sequence mapping, against closed forms."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import compare, sequence

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"
CCD = ROOT / "crates" / "molframe-chem" / "data" / "CCD-amino-acids.cif"


@pytest.fixture(scope="module")
def ubiquitin():
    return molframe.read(BENCH / "1ubq.cif")


@pytest.fixture(scope="module")
def polymer(ubiquitin):
    """Ubiquitin without its waters, which the amino-acid dictionary does not define."""
    editor = ubiquitin.edit()
    editor.clear_extensions()
    editor.delete(ubiquitin.select("not protein"))
    return editor.finish()


@pytest.fixture(scope="module")
def alpha_carbons(ubiquitin):
    atoms = ubiquitin.atoms
    names = [atoms[index].name for index in range(ubiquitin.atom_count)]
    keep = [index for index, name in enumerate(names) if name == "CA"]
    return np.ascontiguousarray(np.asarray(ubiquitin.coordinates)[keep])


def rotation(angle: float) -> np.ndarray:
    return np.array(
        [[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]]
    )


def test_ce_aligns_a_structure_with_a_rigidly_moved_copy_of_itself(alpha_carbons):
    moved = np.ascontiguousarray(
        (alpha_carbons.astype(np.float64) @ rotation(1.0).T + [12.0, -7.0, 3.0]).astype(np.float32)
    )
    best = compare.ce_align(alpha_carbons, moved)[0]
    assert best.rmsd == pytest.approx(0.0, abs=1e-3)
    assert len(best) >= len(alpha_carbons) - 8
    reference = np.asarray(best.reference_indices)
    mobile = np.asarray(best.mobile_indices)
    assert np.array_equal(reference, mobile)
    assert best.z_score is not None
    assert best.z_score > 3.0
    assert compare.ce_align(alpha_carbons, moved, significance=False)[0].z_score is None


def test_ce_recovers_the_offset_between_a_chain_and_its_trimmed_copy(alpha_carbons):
    trimmed = np.ascontiguousarray(alpha_carbons[5:70])
    best = compare.ce_align(alpha_carbons, trimmed)[0]
    reference = np.asarray(best.reference_indices)
    mobile = np.asarray(best.mobile_indices)
    assert best.rmsd == pytest.approx(0.0, abs=1e-3)
    assert np.all(reference - mobile == 5)
    assert len(best) >= 50
    with pytest.raises(molframe.MolframeError):
        compare.ce_align(alpha_carbons[:5].copy(), trimmed)


def two_residue_structure(separation: float):
    positions = np.array([[0.0, 0.0, 0.0], [separation, 0.0, 0.0]], dtype=np.float32)
    radii = np.array([1.5, 1.5], dtype=np.float32)
    residues = np.array([0, 1], dtype=np.uint32)
    return positions, radii, residues


def test_contact_area_between_two_atoms_is_a_pair_of_caps():
    probe, separation, radius = 1.4, 2.0, 1.5
    positions, radii, residues = two_residue_structure(separation)
    table = compare.contact_areas(positions, radii, residues, probe=probe, density=8.0)
    assert len(table) == 1
    assert (int(table["first"][0]), int(table["second"][0])) == (0, 1)
    rho = radius + probe
    # The contact area is the part of one expanded sphere that lies inside the other: a cap
    # of height rho - separation / 2.
    expected = 2.0 * np.pi * rho * (rho - separation / 2.0)
    assert float(table["area"][0]) == pytest.approx(expected, rel=0.05)


def test_cad_scores_one_for_the_same_areas_and_less_for_a_pulled_apart_model():
    reference = compare.contact_areas(*two_residue_structure(2.0), density=8.0)
    assert compare.cad_score(reference, reference).score == pytest.approx(1.0)
    apart = compare.contact_areas(*two_residue_structure(7.0), density=8.0)
    result = compare.cad_score(reference, apart)
    assert len(apart) == 0
    assert result.score == pytest.approx(0.0)
    assert result.lost_area == pytest.approx(result.reference_area)
    contacts = result.contacts()
    assert float(contacts["model_area"][0]) == 0.0
    assert float(contacts["lost_area"][0]) == pytest.approx(float(contacts["reference_area"][0]))
    local = result.local()
    assert np.allclose(local["score"], 0.0)
    assert set(np.asarray(local["residue"]).tolist()) == {0, 1}


def test_cad_caps_the_difference_at_the_reference_area():
    reference = {
        "first": np.array([0], dtype=np.uint32),
        "second": np.array([1], dtype=np.uint32),
        "area": np.array([10.0]),
    }
    bigger = {**reference, "area": np.array([100.0])}
    smaller = {**reference, "area": np.array([4.0])}
    assert compare.cad_score(reference, bigger).score == pytest.approx(0.0)
    assert compare.cad_score(reference, smaller).score == pytest.approx(1.0 - 6.0 / 10.0)
    with pytest.raises(molframe.MolframeError):
        compare.cad_score(reference, {**reference, "area": np.array([-1.0])})


def test_contact_similarity_is_the_jaccard_overlap_regardless_of_pair_order():
    first = np.array([[0, 1], [1, 2], [3, 4]], dtype=np.uint32)
    second = np.array([[1, 0], [2, 1], [5, 6], [7, 8]], dtype=np.uint32)
    result = compare.contact_similarity(first, second)
    assert (result.shared, result.union) == (2, 5)
    assert result.jaccard == pytest.approx(2.0 / 5.0)
    empty = np.zeros((0, 2), dtype=np.uint32)
    assert compare.contact_similarity(empty, empty).jaccard == 1.0
    with pytest.raises(molframe.MolframeValueError):
        compare.contact_similarity(np.zeros((2, 3), dtype=np.uint32), empty)


SCORING = sequence.Scoring(match_score=2, mismatch_score=-1, gap_open=-4, gap_extend=-1)


def test_chains_are_assigned_by_sequence_identity(polymer):
    result = compare.assign_chains(
        polymer,
        polymer,
        components=CCD,
        components_version="test",
        scoring=SCORING,
        min_identity=0.9,
    )
    assert len(result.primary) == 1
    reference, target, identity = result.primary[0]
    assert reference == target
    assert identity == pytest.approx(1.0)
    with pytest.raises(molframe.MolframeError):
        compare.assign_chains(
            polymer,
            polymer,
            components=CCD,
            components_version="test",
            scoring=SCORING,
            min_identity=1.5,
        )


def test_a_query_sequence_is_placed_on_the_residues_of_a_structure(polymer):
    ubq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    placed = compare.map_sequence_to_structure(
        ubq,
        polymer,
        components=CCD,
        components_version="test",
        scoring=SCORING,
    )
    assert len(placed) == len(ubq)
    positions = [position for position, _ in placed]
    residues = [residue for _, residue in placed]
    assert positions == list(range(len(ubq)))
    assert residues == sorted(residues)
    assert len(set(residues)) == len(residues)
