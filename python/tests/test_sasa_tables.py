"""A caller-supplied radius table: parsed in Rust, applied, and reported on the result."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis

BENCH = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"
BASE_ONLY = "name: lab\nversion: 3\nsource: bench notes\nlicence: CC0\nbase: bondi\n"


def test_a_base_only_table_reproduces_the_default_and_names_itself(tmp_path: Path):
    crambin = molframe.read(BENCH / "1crn.cif")
    default = analysis.sasa(crambin)
    path = tmp_path / "radii.txt"
    path.write_text(BASE_ONLY)
    for supplied in (
        analysis.sasa(crambin, radii_table=BASE_ONLY),
        analysis.sasa(crambin, radii_file=path),
    ):
        assert np.array_equal(supplied.value, default.value, equal_nan=True)
        assert any(
            "lab 3 (source: bench notes; licence: CC0)" in note for note in supplied.assumptions
        )
    assert any("bondi" in note for note in default.assumptions)


def test_a_table_value_changes_the_area():
    crambin = molframe.read(BENCH / "1crn.cif")
    wider = analysis.sasa(crambin, radii_table=BASE_ONLY + "C 2.5\n")
    assert not np.array_equal(wider.value, analysis.sasa(crambin).value, equal_nan=True)


def test_a_malformed_table_or_both_arguments_are_refused():
    crambin = molframe.read(BENCH / "1crn.cif")
    with pytest.raises(molframe.MolframeValueError, match="line"):
        analysis.sasa(crambin, radii_table="name: lab\nC -1\n")
    with pytest.raises(molframe.MolframeValueError, match="not both"):
        analysis.sasa(crambin, radii_table=BASE_ONLY, radii_file="x.txt")
