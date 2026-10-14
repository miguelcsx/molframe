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


def test_the_policy_picks_the_radius_table_and_the_result_says_so():
    crambin = molframe.read(BENCH / "1crn.cif")
    by_set = {}
    for name in ("bondi", "charmm", "alvarez", "amber_united"):
        policy = molframe.AnalysisPolicy(vdw_radii=name)
        result = analysis.sasa(crambin, policy=policy)
        by_set[name] = result.value
        assert any(name.replace("_", "-") in note for note in result.assumptions)
    assert not np.array_equal(by_set["bondi"], by_set["charmm"], equal_nan=True)
    explicit = analysis.sasa(crambin, radii="charmm")
    assert np.array_equal(explicit.value, by_set["charmm"], equal_nan=True)


def test_a_radii_argument_must_agree_with_the_policy():
    crambin = molframe.read(BENCH / "1crn.cif")
    policy = molframe.AnalysisPolicy(vdw_radii="charmm")
    agreed = analysis.sasa(crambin, radii="charmm", policy=policy)
    assert np.array_equal(agreed.value, analysis.sasa(crambin, policy=policy).value, equal_nan=True)
    with pytest.raises(molframe.PolicyError, match="E7103"):
        analysis.sasa(crambin, radii="bondi", policy=policy)
