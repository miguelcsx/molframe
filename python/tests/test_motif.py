"""Motif specifications are evaluated at every mapping of a structure."""

import hashlib
from pathlib import Path

import numpy as np
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


PAIR = """
[components.a]
role = "residue"
component_ids = ["SER"]
required_atoms = ["OG"]

[components.b]
role = "residue"
component_ids = ["LYS", "ARG"]
required_atoms = ["CA"]

[[constraints]]
kind = "distance"
name = "reach"
first = "a.OG"
second = "b.CA"
target = 8.0
tolerance = 3.0

[profile]
id = "pair-1.0"
missing = "indeterminate"

[[profile.rules]]
metric = "reach"
comparison = "between"
values = [5.0, 11.0]
"""


class TestLazyMatching:
    def test_default_limits_no_longer_reject_a_single_component_motif(self, tmp_path):
        structure = molframe.read(DATA / "1aon.bcif")
        spec = _write(tmp_path, SPECIFICATION)
        report = motif.evaluate(structure, spec)
        assert len(report) == motif.count(structure, spec, constrained=False)
        assert len(report) > 64

    def test_the_mapping_cap_message_says_mappings_and_how_to_stream(self, tmp_path):
        structure = molframe.read(DATA / "4hhb.pdb")
        spec = _write(tmp_path, SPECIFICATION)
        with pytest.raises(ValueError, match=r"limit of 3 mappings.*iter_matches"):
            motif.evaluate(structure, spec, limits=(3, 256))

    def test_count_matches_and_iteration_agree_with_evaluate(self, tmp_path):
        structure = molframe.read(DATA / "4hhb.pdb")
        spec = _write(tmp_path, PAIR)
        report = motif.evaluate(structure, spec)
        total = motif.count(structure, spec, constrained=False)
        assert total == len(report)
        everything = motif.matches(structure, spec, constrained=False)
        assert everything.residues.shape == (total, 2)
        assert everything.components == ["a", "b"]
        chunks = list(motif.iter_matches(structure, spec, chunk_size=7, constrained=False))
        assert all(c.shape[0] <= 7 for c in chunks)
        assert np.array_equal(np.concatenate(chunks), everything.residues)
        head = motif.first(structure, spec, 5, constrained=False)
        assert np.array_equal(head.residues, everything.residues[:5])

    def test_constrained_equals_the_satisfied_subset_of_the_full_enumeration(self, tmp_path):
        structure = molframe.read(DATA / "4hhb.pdb")
        spec = _write(tmp_path, PAIR)
        report = motif.evaluate(structure, spec)
        satisfied = [i for i, e in enumerate(report.evaluations) if e.constraints[0].satisfied]
        every = motif.matches(structure, spec, constrained=False).residues
        kept = motif.matches(structure, spec, constrained=True).residues
        assert kept.tolist() == every[satisfied].tolist()
        assert motif.count(structure, spec) == len(satisfied)
        constrained = motif.evaluate(structure, spec, constrained=True)
        assert len(constrained) == len(satisfied)

    def test_a_tiny_budget_is_refused_with_the_knob_named(self, tmp_path):
        structure = molframe.read(DATA / "4hhb.pdb")
        spec = _write(tmp_path, PAIR)
        context = molframe.ExecutionContext(memory_budget=64)
        with pytest.raises(molframe.errors.MemoryBudgetError, match="memory_budget"):
            motif.matches(structure, spec, constrained=False, context=context)
        with pytest.raises(molframe.errors.MemoryBudgetError, match="iter_matches"):
            motif.evaluate(structure, spec, context=context)
        # Counting and streaming a small chunk still work.
        assert motif.count(structure, spec, context=context) >= 0
        assert (
            sum(len(c) for c in motif.iter_matches(structure, spec, chunk_size=1, context=context))
            >= 0
        )

    def test_a_cancelled_context_stops_the_work(self, tmp_path):
        structure = molframe.read(DATA / "4hhb.pdb")
        spec = _write(tmp_path, PAIR)
        context = molframe.ExecutionContext()
        context.cancel()
        with pytest.raises(molframe.errors.Cancelled):
            motif.count(structure, spec, context=context)

    def test_a_zero_chunk_is_rejected(self, tmp_path):
        structure = molframe.read(DATA / "4hhb.pdb")
        with pytest.raises(ValueError, match="chunk_size"):
            motif.iter_matches(structure, _write(tmp_path, PAIR), chunk_size=0)


class TestAgainstBruteForce:
    def test_constrained_and_full_sets_equal_an_independent_reference(self, tmp_path):
        structure = molframe.read(DATA / "4hhb.pdb")
        spec = _write(tmp_path, PAIR)
        residues = list(structure.residues)
        sers = [r for r in residues if r.name == "SER" and r.atom("OG") is not None]
        others = [r for r in residues if r.name in ("LYS", "ARG") and r.atom("CA") is not None]

        def point(atom):
            return np.asarray(atom.coordinate, dtype=np.float64)

        full, near = set(), set()
        for a in sers:
            for b in others:
                full.add((a.index, b.index))
                if abs(np.linalg.norm(point(a.atom("OG")) - point(b.atom("CA"))) - 8.0) <= 3.0:
                    near.add((a.index, b.index))
        assert near
        assert len(near) < len(full)

        def digest(rows):
            return hashlib.sha256(repr(sorted(map(tuple, rows))).encode()).hexdigest()

        everything = motif.matches(structure, spec, constrained=False).residues
        kept = motif.matches(structure, spec, constrained=True).residues
        assert digest(everything.tolist()) == digest(full)
        assert digest(kept.tolist()) == digest(near)
