"""Per-residue and streamed SASA: thin conversions of the Rust results."""

import threading
import time
from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis, trajectory

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"

HEADER = (
    "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
    "_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
)
SMALL = HEADER + (
    "ATOM 1 N N ALA A 1 0 0 0\nATOM 2 C CA ALA A 1 1.5 0 0\n"
    "ATOM 3 C C GLY A 2 3.0 0 0\nATOM 4 O O GLY A 2 4.2 0.5 0\n"
    "HETATM 5 O O HOH B 3 2.0 2 0\n"
)
WITH_HYDROGEN = SMALL.replace("ATOM 3 C C GLY", "ATOM 6 H H ALA A 1 0.5 1 0\nATOM 3 C C GLY")


def _small() -> molframe.Structure:
    return molframe.read(SMALL.encode(), name="small.cif")


def test_a_residue_area_is_the_sum_of_its_atoms_areas():
    structure = _small()
    atoms = analysis.sasa(structure).value
    residues = analysis.residue_sasa(structure).value
    assert len(residues) == 2
    assert residues.chain == ["A", "A"]
    assert residues.name == ["ALA", "GLY"]
    assert residues.number.tolist() == [1, 2]
    assert residues.atoms.tolist() == [2, 2]
    assert residues.area[0] == atoms[0] + atoms[1]
    assert residues.area[1] == atoms[2] + atoms[3]
    assert residues.skipped_atoms == 0
    assert residues.missing_radius == 0
    assert not residues.area.flags.writeable


def test_hydrogens_are_not_measured_and_do_not_occlude():
    heavy = analysis.residue_sasa(_small()).value
    with_hydrogen = analysis.residue_sasa(molframe.read(WITH_HYDROGEN.encode(), name="h.cif"))
    assert with_hydrogen.value.skipped_atoms == 1
    assert np.array_equal(with_hydrogen.value.area, heavy.area)


def test_residue_areas_total_the_atom_areas_of_a_protein():
    crambin = molframe.read(BENCH / "1crn.cif")
    atoms = analysis.sasa(crambin).value
    residues = analysis.residue_sasa(crambin).value
    assert np.isclose(residues.area.sum(), np.nansum(atoms), rtol=1e-9)
    assert int(residues.atoms.sum()) == int(np.isfinite(atoms).sum())


def test_the_applied_table_is_reported_and_a_conflict_is_refused():
    structure = _small()
    result = analysis.residue_sasa(structure)
    assert any("bondi" in note for note in result.assumptions)
    policy = molframe.AnalysisPolicy(vdw_radii="charmm")
    charmm = analysis.residue_sasa(structure, policy=policy)
    assert any("charmm" in note for note in charmm.assumptions)
    with pytest.raises(molframe.PolicyError, match="E7103"):
        analysis.residue_sasa(structure, radii="bondi", policy=policy)
    with pytest.raises(molframe.MolframeError, match="E7102"):
        analysis.residue_sasa(structure, radii="no-such-table")


DATA = Path(__file__).parent / "data" / "trajectory"


def _cif(positions: np.ndarray) -> molframe.Structure:
    """One carbon per residue at `positions`, written exactly."""
    rows = [
        f"ATOM {i + 1} C CA ALA A {i + 1} {x!r} {y!r} {z!r}"
        for i, (x, y, z) in enumerate(positions.astype(np.float64).tolist())
    ]
    return molframe.read((HEADER + "\n".join(rows) + "\n").encode(), name="frame.cif")


def test_a_streamed_total_equals_the_batch_total_of_every_frame():
    walk = DATA / "walk.xtc"
    positions = trajectory.read(walk).positions
    streamed = trajectory.sasa_stream(walk, _cif(positions[0]))
    series = streamed.value
    assert streamed.status == "complete"
    assert series.frames.tolist() == list(range(len(positions)))
    for frame, total in zip(positions, series.values, strict=True):
        assert total == pytest.approx(float(np.nansum(analysis.sasa(_cif(frame)).value)), rel=1e-9)
    assert not series.values.flags.writeable
    assert any("bondi" in note for note in streamed.assumptions)
    sliced = trajectory.sasa_stream(walk, _cif(positions[0]), frames=slice(2, None, 3)).value
    assert sliced.frames.tolist() == [2, 5, 8]
    assert np.array_equal(sliced.values, series.values[[2, 5, 8]])


def test_a_stream_follows_the_radii_rule_and_checks_the_atom_count():
    walk = DATA / "walk.xtc"
    topology = _cif(trajectory.read(walk, frames=slice(0, 1)).positions[0])
    policy = molframe.AnalysisPolicy(vdw_radii="charmm")
    charmm = trajectory.sasa_stream(walk, topology, policy=policy)
    assert any("charmm" in note for note in charmm.assumptions)
    with pytest.raises(molframe.PolicyError, match="E7103"):
        trajectory.sasa_stream(walk, topology, radii="bondi", policy=policy)
    with pytest.raises(molframe.MolframeError, match="E7102"):
        trajectory.sasa_stream(walk, topology, radii="no-such-table")
    with pytest.raises(molframe.MolframeError, match="E5102"):
        trajectory.sasa_stream(walk, _small())
    with pytest.raises(molframe.MolframeError, match="E5101"):
        trajectory.sasa_stream(walk, topology, method="lee-richards")


def test_residue_sasa_lets_other_python_threads_run():
    rows = []
    for residue in range(1, 30001):
        for slot, name in enumerate(("N", "CA", "C", "O")):
            x, y = (residue % 150) * 3.8 + slot, (residue // 150) * 3.8
            rows.append(f"ATOM {len(rows) + 1} {name[0]} {name} ALA A {residue} {x} {y} 0")
    big = molframe.read((HEADER + "\n".join(rows) + "\n").encode(), name="big.cif")
    analysis.residue_sasa(big)  # warm up
    ticks = 0
    stop = threading.Event()

    def tick() -> None:
        nonlocal ticks
        while not stop.is_set():
            ticks += 1
            time.sleep(0.001)

    ticker = threading.Thread(target=tick)
    ticker.start()
    try:
        start = time.perf_counter()
        for _ in range(3):
            analysis.residue_sasa(big)
        elapsed = time.perf_counter() - start
    finally:
        stop.set()
        ticker.join()
    assert elapsed > 0.05, "the calculation is too fast to observe"
    assert ticks >= 0.2 * elapsed * 1000
