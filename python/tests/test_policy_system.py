"""Assembly, crystal contacts and hydrogens take part in governed analyses."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis

BENCH = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"


@pytest.fixture(scope="module")
def crambin() -> molframe.Structure:
    return molframe.read(BENCH / "1crn.cif")


def test_every_decision_is_nameable_and_read_back():
    policy = molframe.AnalysisPolicy(
        assembly="biological:1",
        symmetry="biological_assembly",
        model="index:2",
        atom_equivalence="none",
        alignment="explicit:chain A",
        periodic="minimum_image",
        contact_def="surface:1.4",
        float_tolerance=(1e-6, 1e-9),
    )
    assert policy.assembly == "biological:1"
    assert policy.model == "index:2"
    assert policy.alignment == "explicit:chain A"
    assert policy.contact_def == "surface:1.4"
    assert policy.float_tolerance == (1e-6, 1e-9)
    assert policy != molframe.AnalysisPolicy()
    with pytest.raises(molframe.MolframeError):
        molframe.AnalysisPolicy(assembly="crystal:eight")
    with pytest.raises(molframe.MolframeValueError):
        molframe.AnalysisPolicy(float_tolerance=(-1.0, 0.0))


def test_an_assembly_the_unit_does_not_form_is_analysed_when_the_policy_names_it(dimer):
    unit = analysis.contacts(dimer, 4.0)
    assembly = analysis.contacts(
        dimer, 4.0, policy=molframe.AnalysisPolicy(assembly="biological:1")
    )
    # The deposited unit is one carbon and one hydrogen; its assembly holds two copies of both.
    assert unit.coverage.intended == 2
    assert assembly.coverage.intended == 4
    xyz_pairs = len(unit.value), len(assembly.value)
    assert xyz_pairs[1] > xyz_pairs[0]
    assert any("biological assembly 1" in line for line in assembly.assumptions)


def test_hydrogens_can_be_left_out_and_cannot_be_invented(dimer):
    explicit = analysis.contacts(dimer, 2.0)
    excluded = analysis.contacts(dimer, 2.0, policy=molframe.AnalysisPolicy(hydrogens="exclude"))
    assert len(explicit.value) == 1
    assert len(excluded.value) == 0
    assert excluded.coverage.intended == 1
    with pytest.raises(molframe.PolicyError):
        analysis.contacts(dimer, 2.0, policy=molframe.AnalysisPolicy(hydrogens="include_inferred"))


def test_a_system_that_cannot_exist_is_refused_with_the_reason(dimer):
    with pytest.raises(molframe.PolicyError) as missing:
        analysis.contacts(dimer, 4.0, policy=molframe.AnalysisPolicy(assembly="biological:9"))
    assert missing.value.code == "MOLFRAME-E6002"
    with pytest.raises(molframe.PolicyError) as contradiction:
        analysis.contacts(dimer, 4.0, policy=molframe.AnalysisPolicy(symmetry="crystallographic"))
    assert contradiction.value.code == "MOLFRAME-E6004"
    with pytest.raises(molframe.MolframeError) as no_operators:
        analysis.contacts(dimer, 4.0, policy=molframe.AnalysisPolicy(assembly="crystal:4"))
    assert no_operators.value.code == "MOLFRAME-E6016"


def test_crystal_contacts_of_crambin_add_the_neighbours_a_unit_alone_cannot_see(crambin):
    cutoff = 3.5
    unit = analysis.contacts(crambin, cutoff)
    crystal = analysis.contacts(
        crambin, cutoff, policy=molframe.AnalysisPolicy(assembly="crystal:4.0")
    )
    atoms = crambin.atom_count
    assert crystal.coverage.intended > atoms
    first = np.asarray(crystal.value.first, dtype=np.int64)
    second = np.asarray(crystal.value.second, dtype=np.int64)
    cross = int(((first < atoms) != (second < atoms)).sum())
    assert cross > 0
    # The unit's own contacts are unchanged by its neighbours.
    inside = int(((first < atoms) & (second < atoms)).sum())
    assert inside == len(unit.value)
    assert any("crystal contacts" in line for line in crystal.assumptions)


def test_a_result_says_which_input_atom_each_analysed_atom_is(crambin):
    atoms = crambin.atom_count
    unit = analysis.contacts(crambin, 3.5)
    assert np.array_equal(unit.atom_origin, np.arange(atoms))
    crystal = analysis.contacts(
        crambin, 3.5, policy=molframe.AnalysisPolicy(assembly="crystal:4.0")
    )
    origin = np.asarray(crystal.atom_origin)
    assert len(origin) == crystal.coverage.intended > atoms
    # The deposited unit comes first, unchanged; every copy is a copy of one of its atoms.
    assert np.array_equal(origin[:atoms], np.arange(atoms))
    assert origin.max() < atoms
    assert not crystal.atom_origin.flags.writeable  # type: ignore[union-attr]


def test_the_ceiling_on_candidate_crystal_images_belongs_to_the_context(crambin):
    crystal = molframe.AnalysisPolicy(assembly="crystal:4.0")
    tight = molframe.ExecutionContext(image_limit=10)
    assert tight.image_limit == 10
    with pytest.raises(molframe.PolicyError) as refused:
        analysis.contacts(crambin, 3.5, policy=crystal, context=tight)
    assert refused.value.code == "MOLFRAME-E6017"
    generous = molframe.ExecutionContext(image_limit=50_000_000)
    assert len(analysis.contacts(crambin, 3.5, policy=crystal, context=generous).value) > 0
