"""Policy audits: measured against quantities worked out by hand."""

from pathlib import Path

import pytest

import molframe
from molframe import analysis, audit

BENCH = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"


def contact_count(structure: molframe.Structure, cutoff: float):
    """Return the structure's atom contacts as an analysis that depends on the policy."""

    def analyse(policy: molframe.AnalysisPolicy):
        return analysis.contacts(structure, cutoff, policy=policy)

    return analyse


def test_the_variation_among_runs_splits_into_decisions_and_their_interaction(dimer):
    # Contacts within 2 A of the dimer: carbon 0 and hydrogen at 1 A, and a copy 3 A away.
    #   hydrogens explicit, unit      -> 1 (C0-H0)
    #   hydrogens explicit, assembly  -> 3 (C0-H0, C1-H1, H0-C1 at exactly 2 A)
    #   hydrogens excluded, either    -> 0 (the two carbons are 3 A apart)
    # So the assembly matters only while hydrogens are present: an interaction.
    space = (
        audit.PolicySpace()
        .vary("hydrogens", ["explicit_only", "exclude"])
        .vary("assembly", ["asymmetric_unit", "biological:1"])
    )
    assert space.cost == 4
    result = audit.run(
        space.plan(),
        contact_count(dimer, 2.0),
        metric="absolute",
        project=len,
    )
    assert result.metric == "absolute-error"
    assert result.indeterminate == []
    assert result.indeterminate_fraction == 0.0
    # Outcomes in plan order (hydrogens outer): 1, 3, 0, 0.  Mean 1, total SS 6.
    assert result.total_variation == pytest.approx(6.0)
    effects = {effect.field: effect for effect in result.effects or []}
    assert effects["hydrogens"].share == pytest.approx(4.0 / 6.0)
    assert effects["assembly"].share == pytest.approx(1.0 / 6.0)
    assert effects["assembly"].mean_change == pytest.approx(1.0)
    (interaction,) = result.interactions or []
    assert {interaction.first, interaction.second} == {"hydrogens", "assembly"}
    assert interaction.share == pytest.approx(1.0 / 6.0)
    assert result.higher_order == pytest.approx(0.0)
    # Only one of the four universes agrees with the first.
    assert result.agreement_with_first == pytest.approx(0.25)
    assert result.read is not None
    assert {"hydrogens", "assembly"} <= set(result.read)


def test_sets_of_contacts_are_compared_by_overlap(dimer):
    space = audit.PolicySpace().vary("assembly", ["asymmetric_unit", "biological:1"])
    result = audit.run(
        space.plan(),
        contact_count(dimer, 2.0),
        metric="set",
        project=lambda table: list(zip(table.first.tolist(), table.second.tolist(), strict=True)),
    )
    # The unit has one contact; the assembly's three include it: Jaccard distance 1 - 1/3.
    (effect,) = result.effects or []
    assert effect.mean_change == pytest.approx(2.0 / 3.0)


def test_a_conclusion_that_flips_is_counted_not_averaged(dimer):
    space = audit.PolicySpace().vary("hydrogens", ["explicit_only", "exclude"])
    result = audit.run(
        space.plan(),
        contact_count(dimer, 2.0),
        metric="flip",
        project=lambda table: len(table) > 0,
    )
    assert result.effects is not None
    assert result.effects[0].mean_change == pytest.approx(1.0)
    assert result.agreement_with_first == pytest.approx(0.5)


def test_a_decision_the_analysis_never_applied_is_refused_not_reported_stable(dimer):
    # Contacts are defined by a cutoff; they never read the radius set. A sweep over it would
    # come back perfectly stable, and that stability would be an artefact.
    space = audit.PolicySpace().vary("vdw_radii", ["bondi", "charmm"])
    with pytest.raises(molframe.PolicyError) as refused:
        audit.run(space.plan(), contact_count(dimer, 2.0), metric="absolute", project=len)
    assert refused.value.code == "MOLFRAME-E6103"
    assert "vdw_radii" in str(refused.value)


def test_decisions_that_contradict_each_other_cannot_be_planned():
    space = (
        audit.PolicySpace()
        .vary("assembly", ["asymmetric_unit", "biological:1"])
        .vary("symmetry", ["none", "crystallographic"])
    )
    with pytest.raises(molframe.PolicyError) as conflict:
        space.plan()
    assert conflict.value.code == "MOLFRAME-E6004"


def test_the_vocabulary_and_the_size_of_a_space_are_checked_up_front():
    with pytest.raises(molframe.MolframeError):
        audit.PolicySpace().vary("altloc", ["first", "nonsense"])
    with pytest.raises(molframe.MolframeValueError, match="not a policy decision"):
        audit.PolicySpace().vary("colour", ["red"])
    big = audit.PolicySpace(max_runs=3).vary("model", ["first", "index:1", "index:2", "index:3"])
    with pytest.raises(molframe.MolframeError) as limit:
        big.plan()
    assert limit.value.code == "MOLFRAME-E7001"


def test_the_metric_and_the_analysis_are_checked(dimer):
    plan = audit.PolicySpace().vary("hydrogens", ["explicit_only", "exclude"]).plan()
    with pytest.raises(molframe.MolframeValueError, match="metric"):
        audit.run(plan, contact_count(dimer, 2.0), metric="euclid")  # type: ignore[arg-type]
    with pytest.raises(molframe.MolframeTypeError, match="Analysis"):
        audit.run(plan, lambda _policy: 3, metric="absolute")  # type: ignore[arg-type, return-value]


def test_crambin_in_its_crystal_is_a_different_system_from_crambin_alone():
    crambin = molframe.read(BENCH / "1crn.cif")
    space = audit.PolicySpace().vary("assembly", ["asymmetric_unit", "crystal:4.0"])
    result = audit.run(space.plan(), contact_count(crambin, 3.5), metric="absolute", project=len)
    unit, crystal = (len(run.value) for run in result.runs)
    assert crystal > unit
    assert (result.effects or [])[0].mean_change == pytest.approx(crystal - unit)
    assert result.policies[1].assembly == "crystal:4"
