"""Table-backed built-in validation: arrays, summaries and applied-table provenance."""

import os
import subprocess
import sys
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
    assert {"top8000-rama", "top8000-rotamer"} <= names
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


def test_geometry_z_reports_the_peptide_links_of_the_selected_table():
    report = molframe.validation.geometry_z(molframe.read(BENCH / "1crn.cif"))
    assert report.summary["bond_count"] == 45
    assert report.tables[0]["name"] == "alphafold-peptide-link"
    assert report.tables[0]["licence"] == "Apache-2.0"


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


def test_a_missing_user_file_is_an_io_error_with_its_code():
    with pytest.raises(molframe.errors.MolframeIOError, match="E7101"):
        molframe.validation.protein_conformation(
            molframe.read(BENCH / "1crn.cif"), rama_table="/nonexistent/rama.rama"
        )


def test_an_unresolvable_table_names_how_to_supply_it(tmp_path):
    script = (
        "import molframe, sys\n"
        f"s = molframe.read({str(BENCH / '1crn.cif')!r})\n"
        "try:\n"
        "    molframe.validation.protein_conformation(s)\n"
        "except molframe.errors.MolframeError as e:\n"
        "    print(e)\n"
    )
    env = {**os.environ, "MOLFRAME_DATA_DIR": str(tmp_path)}
    env.pop("MOLFRAME_PROFILE", None)
    out = subprocess.run(
        [sys.executable, "-I", "-c", script], env=env, capture_output=True, text=True, check=True
    ).stdout
    assert "E7102" in out
    assert "ramachandran" in out
    assert "molframe data fetch" in out
