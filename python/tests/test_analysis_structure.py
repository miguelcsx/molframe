"""The governed structure analyses, each checked against a NumPy brute force."""

from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"
CCD = ROOT / "crates" / "molframe-chem" / "data" / "CCD-amino-acids.cif"


def residue_atoms(structure):
    return [
        [atom.index for atom in (residue.atoms[i] for i in range(len(residue.atoms)))]
        for residue in (structure.residues[k] for k in range(structure.residue_count))
    ]


@pytest.fixture(scope="module")
def crambin():
    return molframe.read(BENCH / "1crn.cif")


def test_contact_map_equals_a_brute_force_residue_minimum(crambin):
    xyz = np.asarray(crambin.coordinates, dtype=np.float64)
    groups = residue_atoms(crambin)
    cutoff, separation = 4.5, 3
    expected = {}
    for first, second in combinations(range(len(groups)), 2):
        if second - first < separation:
            continue
        gap = xyz[groups[first]][:, None, :] - xyz[groups[second]][None, :, :]
        closest = float(np.sqrt((gap**2).sum(axis=2)).min())
        if closest <= cutoff:
            expected[(first, second)] = closest
    result = analysis.contact_map(crambin, cutoff=cutoff, min_separation=separation)
    table = result.value
    found = {
        (int(a), int(b)): float(d)
        for a, b, d in zip(
            table["first_residue"], table["second_residue"], table["min_distance"], strict=True
        )
    }
    assert result.status == "complete"
    assert found.keys() == expected.keys()
    for pair, distance in expected.items():
        assert found[pair] == pytest.approx(distance, abs=1e-4)


def test_a_contact_map_states_its_cutoff_or_is_refused(crambin):
    with pytest.raises(TypeError):
        analysis.contact_map(crambin, min_separation=1)  # type: ignore[call-arg]


def test_chain_interface_equals_a_brute_force_boundary():
    hemoglobin = molframe.read(BENCH / "4hhb.cif")
    xyz = np.asarray(hemoglobin.coordinates, dtype=np.float64)
    cutoff = 4.5
    groups = residue_atoms(hemoglobin)

    def named(name):
        """Collect the residues of every chain labelled `name` or authored as `name`."""
        found = set()
        for k in range(hemoglobin.chain_count):
            chain = hemoglobin.chains[k]
            if name in {chain.label, chain.auth_label}:
                found.update(chain.residues[i].index for i in range(len(chain.residues)))
        return found

    expected = set()
    for left in named("A"):
        for right in named("B"):
            gap = xyz[groups[left]][:, None, :] - xyz[groups[right]][None, :, :]
            if np.sqrt((gap**2).sum(axis=2)).min() <= cutoff:
                expected.update({left, right})
    result = analysis.chain_interface(hemoglobin, first_chain="A", second_chain="B", cutoff=cutoff)
    assert expected
    assert {int(each) for each in result.value["residue"]} == expected


def test_native_contacts_equal_a_brute_force_pair_count(crambin):
    cutoff, scale = 6.0, 1.6
    same = analysis.native_contacts(crambin, crambin, cutoff=cutoff, tolerance=1.0)
    assert same.value["fraction"] == pytest.approx(1.0)
    xyz = np.asarray(crambin.coordinates, dtype=np.float64)
    distances = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(axis=2))
    upper = np.triu(np.ones_like(distances, dtype=bool), k=1)
    native = int(((distances <= cutoff) & upper).sum())
    assert same.value["native"] == same.value["kept"] == native
    scoped = crambin.edit().coordinates()
    scoped.positions()[:] *= np.float32(scale)
    expanded = analysis.native_contacts(crambin, scoped.commit(), cutoff=cutoff, tolerance=1.0)
    kept = int(((distances <= cutoff) & (distances * scale <= cutoff) & upper).sum())
    assert expanded.value["native"] == native
    assert expanded.value["kept"] == kept
    assert expanded.value["fraction"] == pytest.approx(kept / native)
    assert 0 < kept < native


def test_half_sphere_exposure_counts_alpha_carbons_within_the_radius(crambin, with_roles):
    annotated = with_roles(crambin)
    radius = 13.0
    result = analysis.half_sphere_exposure(annotated, radius=radius)
    table = result.value
    carbons = {}
    for index in range(annotated.residue_count):
        alpha = annotated.residues[index].atom("CA")
        if alpha is not None and annotated.residues[index].atom("CB") is not None:
            carbons[index] = np.array(alpha.coordinate, dtype=np.float64)
    assert len(table) > 0
    for residue, upper, lower in zip(table["residue"], table["upper"], table["lower"], strict=True):
        centre = carbons[int(residue)]
        within = sum(
            1
            for other, point in carbons.items()
            if other != int(residue) and np.linalg.norm(point - centre) <= radius
        )
        assert int(upper) + int(lower) == within


def test_pi_stacking_matches_a_ring_plane_brute_force():
    ubiquitin = molframe.read(BENCH / "1ubq.cif")
    annotated = molframe.chemistry.annotate(ubiquitin, CCD, version="wwPDB-2026-10-03")
    rings = {}
    for index in range(annotated.residue_count):
        residue = annotated.residues[index]
        if residue.name in {"PHE", "TYR"}:
            names = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
            points = np.array([residue.atom(name).coordinate for name in names], dtype=np.float64)
            centre = points.mean(axis=0)
            normal = np.linalg.svd(points - centre)[2][-1]
            rings[index] = (centre, normal)
    limit = 7.0
    expected = {}
    for first, second in combinations(sorted(rings), 2):
        (c1, n1), (c2, n2) = rings[first], rings[second]
        distance = float(np.linalg.norm(c1 - c2))
        if distance <= limit:
            cosine = abs(float(np.dot(n1, n2)))
            expected[(first, second)] = (distance, np.degrees(np.arccos(min(1.0, cosine))))
    result = analysis.pi_stacking(
        annotated,
        max_centre_distance=limit,
        max_parallel_angle=30.0,
        min_perpendicular_angle=60.0,
    )
    table = result.value
    found = {
        (int(a), int(b)): (float(d), float(g))
        for a, b, d, g in zip(
            table["first_residue"],
            table["second_residue"],
            table["centre_distance"],
            table["angle_degrees"],
            strict=True,
        )
        if int(a) in rings and int(b) in rings
    }
    classified = {
        pair: values for pair, values in expected.items() if values[1] <= 30.0 or values[1] >= 60.0
    }
    assert found.keys() == classified.keys()
    for pair, (distance, angle) in classified.items():
        assert found[pair][0] == pytest.approx(distance, abs=1e-3)
        assert found[pair][1] == pytest.approx(angle, abs=0.05)
    assert set(table["kind"]) <= {0, 1}


def test_cation_pi_and_water_bridges_return_their_documented_columns(crambin):
    annotated = molframe.chemistry.annotate(crambin, CCD, version="wwPDB-2026-10-03")
    table = analysis.cation_pi(annotated, max_distance=6.0, max_face_angle=30.0).value
    assert table.names == ["cation_residue", "ring_residue", "distance"]
    assert np.all(np.asarray(table["distance"]) <= 6.0 + 1e-4)
    bridges = analysis.water_bridges(annotated, max_distance=3.5, min_angle=120.0).value
    assert bridges.names == ["water", "first", "second"]


def test_nucleic_torsions_refuse_a_structure_without_roles(crambin):
    with pytest.raises(molframe.PolicyError, match="atom-role") as raised:
        analysis.nucleic_torsions(crambin)
    assert raised.value.code == "MOLFRAME-E6103"


def test_analysis_failures_carry_the_code_of_their_kind(crambin):
    with pytest.raises(molframe.MolframeError) as invalid:
        analysis.contact_map(crambin, cutoff=-1.0, min_separation=0)
    assert invalid.value.code is not None
    annotated = molframe.chemistry.annotate(crambin, CCD, version="wwPDB-2026-10-03")
    with pytest.raises(molframe.MolframeError) as options:
        analysis.pi_stacking(
            annotated,
            max_centre_distance=-1.0,
            max_parallel_angle=30.0,
            min_perpendicular_angle=60.0,
        )
    assert options.value.code == "MOLFRAME-E5101"
    with pytest.raises(molframe.MolframeError) as budget:
        analysis.contact_map(
            crambin,
            cutoff=4.5,
            min_separation=0,
            context=molframe.ExecutionContext(memory_budget=64),
        )
    assert budget.value.code in {"MOLFRAME-E7001", "MOLFRAME-E1902"}
