"""Governed results carry status, coverage, assumptions and provenance."""

import json
from pathlib import Path

import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-py" / "tests" / "data"

# One residue, two conformers of one atom: the altloc policy decides what is analysed.
ALTLOC = b"""data_alt
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
ATOM 1 N N . SER A 1 1 0.0 0.0 0.0 1.0 10.0 1 SER A N
ATOM 2 C CA . SER A 1 1 1.4 0.0 0.0 1.0 10.0 1 SER A CA
ATOM 3 O OG A SER A 1 1 2.0 1.0 0.0 0.6 10.0 1 SER A OG
ATOM 4 O OG B SER A 1 1 2.0 -1.0 0.0 0.4 10.0 1 SER A OG
"""


def test_a_governed_result_names_its_status_coverage_and_provenance():
    structure = molframe.read(DATA / "basic.cif")
    result = molframe.analysis.contacts(structure, 4.0)
    assert isinstance(result, molframe.Analysis)
    assert result.status == "complete"
    assert result.coverage.fraction == pytest.approx(1.0)
    assert result.profile == "molframe-default-1.0"
    record = json.loads(result.provenance)
    assert record["molframe_version"]
    assert len(result.value) == len(molframe.analysis.atom_contacts(structure, 4.0))


def test_an_altered_policy_is_recorded_and_changes_the_analysed_atoms():
    structure = molframe.read(ALTLOC, name="alt.cif")
    keep_all = molframe.AnalysisPolicy(altloc="keep_all")
    first = molframe.AnalysisPolicy(altloc="first")
    everything = molframe.analysis.contacts(structure, 5.0, policy=keep_all)
    one_conformer = molframe.analysis.contacts(structure, 5.0, policy=first)
    assert everything.profile is None
    assert one_conformer.profile is None
    assert json.loads(everything.provenance) != json.loads(one_conformer.provenance)
    assert len(one_conformer.value) < len(everything.value)


def test_validation_returns_the_same_envelope():
    structure = molframe.read(DATA / "basic.cif")
    result = molframe.validation.clashes(structure, tolerance=0.0)
    assert result.status in {"complete", "partial", "ambiguous"}
    assert result.coverage.intended >= result.coverage.used
