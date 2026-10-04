"""lDDT, DockQ and QS, checked against NumPy and against their own invariances."""

from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import compare

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"
CCD = ROOT / "crates" / "molframe-chem" / "data" / "CCD-amino-acids.cif"
rng = np.random.default_rng(7)


def reference_lddt(model, reference, radius, tolerances):
    m, r = model.astype(np.float64), reference.astype(np.float64)
    kept, total = 0.0, 0
    for i, j in combinations(range(len(r)), 2):
        d = np.linalg.norm(r[i] - r[j])
        if d <= radius:
            deviation = abs(np.linalg.norm(m[i] - m[j]) - d)
            kept += sum(deviation <= tolerance for tolerance in tolerances) / len(tolerances)
            total += 1
    return kept / total if total else 1.0


def test_lddt_equals_a_brute_force_local_distance_test():
    reference = (rng.normal(size=(60, 3)) * 6.0).astype(np.float32)
    model = reference + rng.normal(scale=0.8, size=reference.shape).astype(np.float32)
    for radius in (8.0, 15.0):
        found = compare.lddt(model, reference, inclusion_radius=radius)
        assert found == pytest.approx(
            reference_lddt(model, reference, radius, (0.5, 1.0, 2.0, 4.0)), abs=1e-9
        )
    custom = compare.lddt(model, reference, inclusion_radius=10.0, tolerances=[1.0, 3.0])
    assert custom == pytest.approx(reference_lddt(model, reference, 10.0, (1.0, 3.0)), abs=1e-9)
    assert compare.lddt(reference, reference) == pytest.approx(1.0)
    assert compare.lddt(model[:, :], reference[: len(model)]) <= 1.0


def test_lddt_is_invariant_to_rigid_motion_unlike_rmsd():
    reference = (rng.normal(size=(50, 3)) * 5.0).astype(np.float32)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    moved = (reference @ (q * np.sign(np.linalg.det(q))).T + np.float32(11.0)).astype(np.float32)
    assert compare.lddt(moved, reference) == pytest.approx(1.0, abs=1e-4)
    assert molframe.geometry.rmsd(moved, reference) > 1.0


@pytest.fixture(scope="module")
def hemoglobin():
    """Keep the alpha (A) and beta (B) chains of 4HHB, whose residues the dictionary knows."""
    whole = molframe.read(BENCH / "4hhb.cif")
    editor = whole.edit()
    editor.clear_extensions()
    editor.delete(whole.select("not (protein and (chain A or chain B))"))
    return editor.finish()


def renamed(structure, mapping):
    editor = structure.edit()
    editor.clear_extensions()
    for index in range(structure.chain_count):
        label = structure.chains[index].label
        if label in mapping:
            editor.rename_chain(index, mapping[label])
    return editor.finish()


def displaced(structure, chain, shift):
    scoped = structure.edit().coordinates()
    positions = scoped.positions()
    owned = [
        atom.index
        for k in range(structure.chain_count)
        if structure.chains[k].label == chain
        for r in range(len(structure.chains[k].residues))
        for atom in (
            structure.chains[k].residues[r].atoms[a]
            for a in range(len(structure.chains[k].residues[r].atoms))
        )
    ]
    positions[owned] += np.asarray(shift, dtype=np.float32)
    return scoped.commit()


def test_dockq_of_a_complex_with_itself_is_perfect(hemoglobin):
    score = compare.dockq(hemoglobin, hemoglobin, receptor="A", ligand="B")
    assert score.score == pytest.approx(1.0)
    assert score.fnat == pytest.approx(1.0)
    assert score.ligand_rmsd == pytest.approx(0.0, abs=1e-4)
    assert score.interface_rmsd == pytest.approx(0.0, abs=1e-4)


def test_dockq_falls_as_the_ligand_leaves_the_interface(hemoglobin):
    near = compare.dockq(
        displaced(hemoglobin, "B", (1.0, 0.0, 0.0)), hemoglobin, receptor="A", ligand="B"
    )
    far = compare.dockq(
        displaced(hemoglobin, "B", (20.0, 0.0, 0.0)), hemoglobin, receptor="A", ligand="B"
    )
    assert 1.0 > near.score > far.score > 0.0
    assert near.fnat > far.fnat
    assert far.fnat == pytest.approx(0.0)
    assert near.ligand_rmsd == pytest.approx(1.0, abs=0.05)
    assert far.ligand_rmsd == pytest.approx(20.0, abs=0.5)


def test_dockq_states_its_constants_and_refuses_nonsense(hemoglobin):
    tighter = compare.dockq(
        displaced(hemoglobin, "B", (1.0, 0.0, 0.0)),
        hemoglobin,
        receptor="A",
        ligand="B",
        contact_distance=3.5,
    )
    assert tighter.fnat < 1.0
    with pytest.raises(molframe.MolframeError):
        compare.dockq(hemoglobin, hemoglobin, receptor="A", ligand="B", contact_distance=-1.0)
    with pytest.raises(molframe.PolicyError):
        compare.dockq(hemoglobin, hemoglobin, receptor="A", ligand="B", chain_names="both")


def test_a_renamed_model_scores_the_same_through_the_mapping(hemoglobin):
    model = renamed(hemoglobin, {"A": "X", "B": "Y"})
    scoring = molframe.sequence.Scoring(
        match_score=2, mismatch_score=-1, gap_open=-2, gap_extend=-1
    )
    direct = compare.dockq(hemoglobin, hemoglobin, receptor="A", ligand="B")
    score, found = compare.mapped_dockq(
        model,
        hemoglobin,
        receptor="A",
        ligand="B",
        components=CCD,
        components_version="wwPDB-2026-10-03",
        scoring=scoring,
        min_identity=0.9,
        automorphism_limit=64,
    )
    assert score.score == pytest.approx(direct.score, abs=1e-9)
    assert score.fnat == pytest.approx(direct.fnat, abs=1e-9)
    assert {(native, target) for native, target, _ in found.chains} == {("A", "X"), ("B", "Y")}
    assert found.atoms > 1000
    assert found.chain_assignments_tried >= 1
    pair, _ = compare.mapped_qs_score(
        model,
        hemoglobin,
        first_chain="A",
        second_chain="B",
        components=CCD,
        components_version="wwPDB-2026-10-03",
        scoring=scoring,
        min_identity=0.9,
        automorphism_limit=64,
    )
    assert pair == pytest.approx(1.0)


def test_an_unresolvable_component_names_itself_in_the_error():
    whole = molframe.read(BENCH / "4hhb.cif")
    scoring = molframe.sequence.Scoring()
    with pytest.raises(molframe.MolframeError, match="component=") as raised:
        compare.mapped_qs_score(
            whole,
            whole,
            first_chain="A",
            second_chain="B",
            components=CCD,
            components_version="v",
            scoring=scoring,
            min_identity=0.9,
            automorphism_limit=8,
        )
    assert raised.value.findings[0].context["component"]


def test_qs_score_is_one_for_identity_and_less_for_a_displaced_chain(hemoglobin):
    assert compare.qs_score(
        hemoglobin, hemoglobin, first_chain="A", second_chain="B"
    ) == pytest.approx(1.0)
    moved = displaced(hemoglobin, "B", (20.0, 0.0, 0.0))
    assert compare.qs_score(moved, hemoglobin, first_chain="A", second_chain="B") < 0.5
