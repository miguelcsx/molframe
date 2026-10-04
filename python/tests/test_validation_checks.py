"""Validation checks, each compared with an independent computation."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import geometry, validation

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"
CCD = ROOT / "crates" / "molframe-chem" / "data" / "CCD-amino-acids.cif"


def pdb_line(  # noqa: PLR0913, PLR0917
    serial, name, residue, number, xyz, occupancy=1.0, bfactor=10.0, element="C", alt=" "
):
    x, y, z = xyz
    return (
        f"ATOM  {serial:5d} {name:<4s}{alt}{residue:>3s} A{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}\n"
    )


def read_text(text, **options: object):
    return molframe.read(
        text.encode(), name="t.pdb", options=molframe.ReadOptions(**options) if options else None
    )


@pytest.fixture(scope="module")
def ubiquitin():
    return molframe.read(BENCH / "1ubq.cif")


def test_quality_flags_name_each_kind_of_problem():
    text = (
        pdb_line(1, "N", "ALA", 1, (0, 0, 0), occupancy=0.0, element="N")
        + pdb_line(2, "CA", "ALA", 1, (1.4, 0, 0))
        + pdb_line(3, "C", "ALA", 1, (2.4, 1, 0), bfactor=-3.0)
        + pdb_line(4, "O", "ALA", 1, (3.6, 1, 0), element="O")
        + "END\n"
    )
    result = validation.quality_flags(read_text(text))
    table = result.value
    assert dict(zip(map(int, table["atom"]), map(int, table["issue"]), strict=True)) == {
        0: 0,
        2: 2,
    }
    assert result.status == "complete"


def test_an_occupancy_above_one_makes_the_structure_invalid_and_says_so():
    structure = read_text(pdb_line(1, "N", "ALA", 1, (0, 0, 0), occupancy=1.5, element="N"))
    with pytest.raises(molframe.MolframeError) as raised:
        validation.quality_flags(structure)
    assert raised.value.code == "MOLFRAME-E3009"
    assert "occupancy" in raised.value.message


def test_b_factor_distribution_equals_the_population_statistics(ubiquitin):
    values = np.asarray(ubiquitin.b_factors, dtype=np.float64)
    result = validation.b_factor_distribution(ubiquitin, outlier_standard_deviations=3.0)
    found = result.value
    assert found["assessed"] == len(values)
    assert found["mean"] == pytest.approx(values.mean())
    assert found["standard_deviation"] == pytest.approx(values.std())
    assert found["variance"] == pytest.approx(values.var())
    assert found["median"] == pytest.approx(np.median(values))
    assert (found["minimum"], found["maximum"]) == pytest.approx((values.min(), values.max()))
    z = (values - values.mean()) / values.std()
    expected = set(np.flatnonzero(np.abs(z) > 3.0).tolist())
    assert {int(a) for a in found["outliers"]["atom"]} == expected
    chosen = ubiquitin.select("name CA")
    subset = validation.b_factor_distribution(
        ubiquitin, outlier_standard_deviations=3.0, selection=chosen
    ).value
    assert subset["assessed"] == len(chosen)
    assert subset["mean"] == pytest.approx(values[chosen.indices].mean())


def test_altloc_occupancy_sums_flag_groups_that_do_not_total_one():
    text = (
        pdb_line(1, "CA", "ALA", 1, (0, 0, 0), occupancy=0.6, alt="A")
        + pdb_line(2, "CA", "ALA", 1, (0.1, 0, 0), occupancy=0.4, alt="B")
        + pdb_line(3, "CB", "ALA", 1, (1, 0, 0), occupancy=0.5, alt="A")
        + pdb_line(4, "CB", "ALA", 1, (1.1, 0, 0), occupancy=0.3, alt="B")
        + "END\n"
    )
    result = validation.altloc_occupancy_sums(
        read_text(text, mode="permissive"), expected_sum=1.0, tolerance=0.01
    )
    flagged = {
        record["atom_name"]: record["issue"]
        for record in result.value["records"]
        if record["issue"] is not None
    }
    assert flagged == {"CB": "sum_mismatch"}
    sums = {record["atom_name"]: record["occupancy_sum"] for record in result.value["records"]}
    assert sums["CA"] == pytest.approx(1.0)
    assert sums["CB"] == pytest.approx(0.8)


def test_valence_flags_an_atom_with_too_many_bonds():
    centre = pdb_line(1, "C1", "LIG", 1, (0, 0, 0))
    arms = [(1.5, 0, 0), (-1.5, 0, 0), (0, 1.5, 0), (0, -1.5, 0), (0, 0, 1.5)]
    atoms = "".join(
        pdb_line(i + 2, f"N{i}", "LIG", 1, xyz, element="N") for i, xyz in enumerate(arms)
    )
    conect = "CONECT    1    2    3    4    5    6\n"
    structure = read_text(centre + atoms + conect + "END\n")
    found = validation.valence(structure).value
    assert [int(a) for a in found["atom"]] == [0]
    assert (int(found["bonds"][0]), int(found["maximum"][0])) == (5, 4)


def test_bond_length_deviations_report_distances_the_coordinates_confirm(ubiquitin):
    result = validation.bond_length_deviations(ubiquitin, tolerance=0.05)
    table = result.value
    xyz = np.asarray(ubiquitin.coordinates, dtype=np.float64)
    a, b = table["atom_a"].astype(int), table["atom_b"].astype(int)
    assert np.allclose(np.linalg.norm(xyz[a] - xyz[b], axis=1), table["observed"], atol=1e-4)
    assert np.allclose(table["observed"] - table["expected"], table["deviation"], atol=1e-5)
    assert np.all(np.abs(table["deviation"]) > 0.05)
    looser = validation.bond_length_deviations(ubiquitin, tolerance=0.2).value
    assert len(looser) < len(table)
    assert {tuple(row) for row in zip(looser["atom_a"], looser["atom_b"], strict=True)} <= {
        tuple(row) for row in zip(table["atom_a"], table["atom_b"], strict=True)
    }


def reference_omegas(structure):
    """Compute every peptide omega from CA, C, N, CA of consecutive residues in NumPy."""
    omegas = {}
    for index in range(structure.residue_count - 1):
        here, after = structure.residues[index], structure.residues[index + 1]
        atoms = [here.atom("CA"), here.atom("C"), after.atom("N"), after.atom("CA")]
        if any(atom is None for atom in atoms) or here.chain.index != after.chain.index:
            continue
        points = [np.array([atom.coordinate], dtype=np.float32) for atom in atoms]
        omegas[index] = float(geometry.dihedrals(*points)[0])
    return omegas


def test_cis_peptides_equal_the_omega_computed_independently(ubiquitin, with_roles):
    annotated = with_roles(ubiquitin)
    omegas = reference_omegas(ubiquitin)
    assert len(omegas) > 60
    for threshold in (30.0, 170.0):
        found = validation.cis_peptides(annotated, threshold_degrees=threshold).value
        expected = {r for r, omega in omegas.items() if abs(omega) <= threshold}
        assert {int(r) for r in found["residue"]} == expected
        for residue, omega in zip(found["residue"], found["omega"], strict=True):
            assert omega == pytest.approx(omegas[int(residue)], abs=1e-3)
    assert len(validation.cis_peptides(annotated, threshold_degrees=30.0).value) == 0


def test_planarity_of_aromatic_rings_matches_a_plane_fit(ubiquitin):
    annotated = molframe.chemistry.annotate(ubiquitin, CCD, version="wwPDB-2026-10-03")
    rings = {}
    for index in range(annotated.residue_count):
        residue = annotated.residues[index]
        if residue.name in {"PHE", "TYR"}:
            points = np.array(
                [
                    residue.atom(name).coordinate
                    for name in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
                ],
                dtype=np.float32,
            )
            rings[index] = geometry.plane_deviation(points)
    assert rings
    tight = validation.planarity(annotated, max_deviation=0.0).value
    found = {int(r): float(d) for r, d in zip(tight["residue"], tight["deviation"], strict=True)}
    for residue, deviation in rings.items():
        assert found[residue] == pytest.approx(deviation, abs=1e-3)
    loose = validation.planarity(annotated, max_deviation=0.5).value
    assert len(loose) < len(tight)


def test_completeness_counts_modelled_residues_and_names_the_missing_ones(ubiquitin):
    complete = validation.completeness(ubiquitin).value
    chain_a = next(chain for chain in complete if chain["chain"] == "A")
    assert chain_a["observed"] == chain_a["canonical"] == 76
    assert chain_a["missing"] == []
    editor = ubiquitin.edit()
    editor.clear_extensions()
    editor.delete(ubiquitin.select("chain A and resid 70:73"))
    gapped = editor.finish()
    found = validation.completeness(gapped).value
    chain = next(entry for entry in found if entry["chain"] == "A")
    assert chain["observed"] == 72
    assert [position for position, _ in chain["missing"]] == [70, 71, 72, 73]
    assert [name for _, name in chain["missing"]] == ["VAL", "LEU", "ARG", "LEU"]


def test_ligand_geometry_assesses_heterogen_bonds_and_tightens_with_tolerance():
    hemoglobin = molframe.read(BENCH / "4hhb.cif")
    strict = validation.ligand_geometry(hemoglobin, tolerance=0.01).value
    loose = validation.ligand_geometry(hemoglobin, tolerance=0.5).value
    assert strict["intended"] > 0
    assert strict["assessed"] <= strict["intended"]
    assert len(loose["outliers"]) <= len(strict["outliers"])


def test_checks_state_their_tolerances_and_classify_their_failures(ubiquitin):
    with pytest.raises(TypeError):
        validation.bond_length_deviations(ubiquitin)  # type: ignore[call-arg]
    with pytest.raises(molframe.MolframeError) as raised:
        validation.b_factor_distribution(ubiquitin, outlier_standard_deviations=-1.0)
    assert raised.value.code == "MOLFRAME-E5101"
