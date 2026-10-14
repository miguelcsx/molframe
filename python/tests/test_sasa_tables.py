"""A caller-supplied radius table: parsed in Rust, applied, and reported on the result."""

import re
from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"
BONDI = (ROOT / "data" / "bondi.radii").read_text()


def _lab(tmp_path: Path, carbon: str | None = None) -> Path:
    text = (
        BONDI.replace("name: bondi", "name: lab")
        .replace("version: 1964", "version: 3")
        .replace("licence: MIT", "licence: CC0")
    )
    path = tmp_path / "lab.radii"
    if carbon is not None:
        text = re.sub(r"(?m)^C .*$", f"C {carbon}", text)
    path.write_text(text)
    return path


def test_a_copy_of_the_default_reproduces_it_and_names_itself(tmp_path: Path):
    crambin = molframe.read(BENCH / "1crn.cif")
    default = analysis.sasa(crambin)
    supplied = analysis.sasa(crambin, radii=str(_lab(tmp_path)))
    assert np.array_equal(supplied.value, default.value, equal_nan=True)
    assert any("lab 3" in note for note in supplied.assumptions)
    assert any("bondi" in note for note in default.assumptions)


def test_a_table_value_changes_the_area(tmp_path: Path):
    crambin = molframe.read(BENCH / "1crn.cif")
    wider = analysis.sasa(crambin, radii=str(_lab(tmp_path, "2.5")))
    assert not np.array_equal(wider.value, analysis.sasa(crambin).value, equal_nan=True)


def test_a_malformed_or_absent_table_is_refused_with_its_code(tmp_path: Path):
    crambin = molframe.read(BENCH / "1crn.cif")
    bad = tmp_path / "bad.radii"
    bad.write_text("format: 1\nname: lab\nversion: 1\nsource: s\nlicence: CC0\nC -1\n")
    with pytest.raises(molframe.MolframeError, match="E2201"):
        analysis.sasa(crambin, radii=str(bad))
    with pytest.raises(molframe.MolframeError, match="E7102"):
        analysis.sasa(crambin, radii="no-such-table")
