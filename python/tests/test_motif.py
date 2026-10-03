"""Motif specifications are evaluated at every mapping of a structure."""

from pathlib import Path

import pytest

import molframe
from molframe import motif

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"

SPECIFICATION = """
[components.site]
role = "residue"
component_ids = ["SER"]
required_atoms = ["OG"]

[[constraints]]
kind = "distance"
name = "og_cb"
first = "site.OG"
second = "site.CB"
target = 1.43
tolerance = 0.2

[profile]
id = "example-1.0"
missing = "indeterminate"

[[profile.rules]]
metric = "og_cb"
comparison = "between"
values = [1.2, 1.7]
"""


def _write(tmp_path, text, name="spec.toml"):
    path = tmp_path / name
    path.write_text(text)
    return path


class TestMotif:
    def test_a_serine_motif_is_evaluated_at_every_mapping(self, tmp_path):
        structure = molframe.read(DATA / "1ubq.pdb")
        report = motif.evaluate(structure, _write(tmp_path, SPECIFICATION))

        assert len(report) > 1
        assert report.mapping_ambiguous
        first = report.evaluations[0]
        assert first.profile == "example-1.0"
        assert first.verdict == "pass"
        assert first.constraints[0].name == "og_cb"
        assert first.constraints[0].value == pytest.approx(1.43, abs=0.1)
        assert [e.mapping_index for e in report.evaluations] == list(range(len(report)))

    def test_a_bad_specification_or_limit_raises_value_error(self, tmp_path):
        structure = molframe.read(DATA / "1ubq.pdb")
        bad = _write(tmp_path, "not = [valid", "broken.toml")
        with pytest.raises(ValueError, match="invalid TOML functional specification"):
            motif.evaluate(structure, bad)
        good = _write(tmp_path, SPECIFICATION)
        with pytest.raises(ValueError, match="limits must be positive"):
            motif.evaluate(structure, good, limits=(0, 4))
        with pytest.raises(ValueError, match="functional specification I/O failed"):
            motif.evaluate(structure, tmp_path / "spec.yaml")
