"""Table-backed built-in validation: arrays, summaries and applied-table provenance."""

from pathlib import Path

import numpy as np
import pytest

import molframe

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"
STEM = ROOT / "crates" / "molframe-validate" / "tests" / "data" / "1ehz_stem.cif"


def test_ubiquitin_has_no_ramachandran_outliers_and_reports_its_tables():
    report = molframe.validation.protein_conformation(molframe.read(BENCH / "1ubq.cif"))
    assert report.summary["rama_outliers"] == 0
    assert report.summary["rama_evaluated"] > 60
    names = {table["name"] for table in report.tables}
    assert {"top8000-ramachandran", "top8000-rotamer"} <= names
    assert all(table["licence"] for table in report.tables)


def test_per_residue_arrays_share_one_length_and_chi_is_four_wide():
    structure = molframe.read(BENCH / "1ubq.cif")
    report = molframe.validation.protein_conformation(structure)
    phi = report["phi"]
    assert isinstance(phi, np.ndarray)
    assert phi.dtype == np.float32
    assert report["chi"].shape == (phi.shape[0], 4)
    assert report["rama_class"].shape == phi.shape
    assert report.arrays["phi"] is not None


def test_geometry_z_matches_the_bond_and_angle_counts_of_a_small_protein():
    report = molframe.validation.geometry_z(molframe.read(BENCH / "1crn.cif"))
    assert 300 < report.summary["bond_count"] < 400
    assert report.summary["bond_rmsz"] < 3
    assert report.tables[0]["name"].startswith("engh-huber")


def test_geometry_z_rejects_a_bad_cutoff():
    with pytest.raises(molframe.errors.MolframeValueError):
        molframe.validation.geometry_z(molframe.read(BENCH / "1crn.cif"), cutoff=-1.0)


def test_clashscore_records_its_policy():
    report = molframe.validation.clashscore(molframe.read(BENCH / "4hhb.cif"))
    assert report.summary["score"] > 0
    assert report.summary["hydrogen_aware"] == 0
    assert "radii=bondi" in report.tables[0]["source"]
    assert report.summary["tolerance"] == pytest.approx(0.4)
    with pytest.raises(molframe.errors.MolframeValueError):
        molframe.validation.clashscore(molframe.read(BENCH / "4hhb.cif"), hydrogens="maybe")  # pyright: ignore[reportArgumentType]


def test_nucleic_stem_is_mostly_c3_endo():
    report = molframe.validation.nucleic_conformation(molframe.read(STEM))
    summary = report.summary
    assert summary["c3_endo"] >= 6
    assert report["nu"].shape[1] == 5


def test_a_missing_user_table_is_a_value_error():
    with pytest.raises(molframe.errors.MolframeValueError):
        molframe.validation.protein_conformation(
            molframe.read(BENCH / "1crn.cif"), rama_table="/nonexistent/rama.bin"
        )
