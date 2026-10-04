"""Policy audits: measured against quantities worked out by hand."""

import json
from pathlib import Path

import numpy as np
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
        project=lambda analysis: len(analysis.value),
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
        project=lambda analysis: list(
            zip(analysis.value.first.tolist(), analysis.value.second.tolist(), strict=True)
        ),
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
        project=lambda analysis: len(analysis.value) > 0,
    )
    assert result.effects is not None
    assert result.effects[0].mean_change == pytest.approx(1.0)
    assert result.agreement_with_first == pytest.approx(0.5)


def test_a_decision_the_analysis_never_applied_is_refused_not_reported_stable(dimer):
    # Contacts are defined by a cutoff; they never read the radius set. A sweep over it would
    # come back perfectly stable, and that stability would be an artefact.
    space = audit.PolicySpace().vary("vdw_radii", ["bondi", "charmm"])
    with pytest.raises(molframe.PolicyError) as refused:
        audit.run(
            space.plan(),
            contact_count(dimer, 2.0),
            metric="absolute",
            project=lambda analysis: len(analysis.value),
        )
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
    result = audit.run(
        space.plan(),
        contact_count(crambin, 3.5),
        metric="absolute",
        project=lambda analysis: len(analysis.value),
    )
    unit, crystal = (len(run.value) for run in result.runs)
    assert crystal > unit
    assert (result.effects or [])[0].mean_change == pytest.approx(crystal - unit)
    assert result.policies[1].assembly == "crystal:4"


def pair_set(table: object) -> set[tuple[int, int]]:
    first, second = table["first"].tolist(), table["second"].tolist()  # type: ignore[attr-defined]
    return set(zip(first, second, strict=True))


def test_contacts_by_definition_equal_a_numpy_distance_test_against_the_named_radii():
    crambin = molframe.read(BENCH / "1crn.cif")
    xyz = np.asarray(crambin.coordinates, dtype=np.float64)
    for radii, tolerance in (("bondi", 0.5), ("charmm", 0.0), ("alvarez", 0.3)):
        policy = molframe.AnalysisPolicy(vdw_radii=radii, contact_def=f"distance:{tolerance}")
        found = pair_set(analysis.contacts_by_definition(crambin, policy=policy).value)
        r = np.asarray(molframe.chemistry.vdw_radii(crambin, radii=radii), dtype=np.float64)
        gap = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
        reach = r[:, None] + r[None, :] + tolerance
        upper = np.triu(gap <= reach, k=1)
        expected = {(int(i), int(j)) for i, j in zip(*np.nonzero(upper), strict=True)}
        assert found == expected, (radii, tolerance)


def test_the_contact_definition_and_the_radii_are_decisions_an_audit_can_vary():
    crambin = molframe.read(BENCH / "1crn.cif")
    space = (
        audit.PolicySpace()
        .vary("contact_def", ["distance:0.0", "distance:0.5"])
        .vary("vdw_radii", ["bondi", "charmm"])
    )

    def analyse(policy: molframe.AnalysisPolicy):
        return analysis.contacts_by_definition(crambin, policy=policy)

    result = audit.run(
        space.plan(),
        analyse,
        metric="set",
        project=lambda analysis: sorted(pair_set(analysis.value)),
    )
    assert {"contact_def", "vdw_radii"} <= set(result.read)
    effects = {effect.field: effect for effect in result.effects or []}
    # Both decisions change which atoms count as touching, so neither is inert.
    assert effects["contact_def"].mean_change > 0.0
    assert effects["vdw_radii"].mean_change > 0.0
    # Widening the tolerance only ever adds contacts: the narrower set is contained.
    narrow, wide = (pair_set(run.value) for run in result.runs[:2])
    assert narrow < wide
    assert result.agreement_with_first is not None


def test_shapley_shares_of_the_dimer_audit_are_the_ones_worked_out_by_hand(dimer):
    # Outcomes in plan order: 1, 3, 0, 0.  Holding a decision fixed removes 4/6 of the
    # variation for hydrogens, 1/6 for the assembly and all of it for both; averaging each
    # decision's marginal gain over the two orders gives 3/4 and 1/4.
    space = (
        audit.PolicySpace()
        .vary(
            "hydrogens",
            ["explicit_only", "exclude"],
            rationale="structures differ in whether hydrogens were modelled",
            evidence="PDBbind protein files carry added hydrogens",
        )
        .vary("assembly", ["asymmetric_unit", "biological:1"])
    )
    plan = space.plan()
    assert plan.balanced
    assert plan.skipped == 0
    hydrogens, assembly = plan.decisions
    assert hydrogens.uncertainty == "interpretive"
    assert "modelled" in hydrogens.rationale
    assert "PDBbind" in hydrogens.evidence
    assert assembly.rationale == ""
    result = audit.run(
        plan,
        contact_count(dimer, 2.0),
        metric="absolute",
        project=lambda analysis: len(analysis.value),
    )
    shares = {share.name: share.share for share in result.shapley or []}
    assert shares["hydrogens"] == pytest.approx(0.75)
    assert shares["assembly"] == pytest.approx(0.25)
    assert sum(shares.values()) == pytest.approx(1.0)
    (interpretive,) = result.by_class or []
    assert (interpretive.name, interpretive.share) == ("interpretive", pytest.approx(1.0))
    assert result.balanced is True


def test_a_constrained_plan_attributes_by_shapley_when_the_product_is_not_whole(dimer):
    space = (
        audit.PolicySpace()
        .vary("assembly", ["asymmetric_unit", "biological:1"])
        .vary("symmetry", ["none", "crystallographic"])
    )
    plan = space.plan(constrained=True)
    assert (plan.cost, plan.skipped, plan.balanced) == (2, 2, False)
    result = audit.run(
        plan,
        contact_count(dimer, 2.0),
        metric="absolute",
        project=lambda analysis: len(analysis.value),
    )
    assert result.balanced is False
    assert sum(share.share for share in result.shapley or []) == pytest.approx(1.0)
    assert result.indeterminate == []


def test_a_forbidden_pair_is_dropped_from_a_constrained_plan_and_refuses_a_strict_one():
    space = (
        audit.PolicySpace()
        .vary("hydrogens", ["explicit_only", "exclude"])
        .vary("assembly", ["asymmetric_unit", "biological:1"])
        .forbid(("hydrogens", "exclude"), ("assembly", "biological:1"))
    )
    with pytest.raises(molframe.PolicyError) as strict:
        space.plan()
    assert strict.value.code == "MOLFRAME-E6004"
    plan = space.plan(constrained=True)
    assert (plan.cost, plan.skipped) == (3, 1)
    with pytest.raises(molframe.MolframeError):
        space.forbid(("hydrogens", "nonsense"), ("assembly", "biological:1"))


def references(node: object) -> set[str]:
    if isinstance(node, dict):
        found = {node["@id"]} if set(node) == {"@id"} else set()
        for value in node.values():
            found |= references(value)
        return found
    if isinstance(node, list):
        return set().union(*(references(item) for item in node)) if node else set()
    return set()


def test_an_audit_writes_a_certificate_that_is_a_closed_ro_crate_graph(dimer):
    space = (
        audit.PolicySpace()
        .vary(
            "hydrogens",
            ["explicit_only", "exclude"],
            rationale="the deposited entry has none, the prepared file has them added",
            evidence="PDBbind v2020 protein files",
        )
        .vary("assembly", ["asymmetric_unit", "biological:1"])
    )
    result = audit.run(
        space.plan(),
        contact_count(dimer, 2.0),
        metric="absolute",
        project=lambda analysis: len(analysis.value),
    )
    document = json.loads(result.certificate)
    assert document["@context"] == "https://w3id.org/ro/crate/1.1/context"
    graph = {node["@id"]: node for node in document["@graph"]}
    assert graph["ro-crate-metadata.json"]["about"] == {"@id": "./"}
    assert len([key for key in graph if key.startswith("#universe-")]) == 4
    # Every local reference points at an entity the document defines; the specification
    # the crate conforms to is an external IRI.
    local = {ref for ref in references(document["@graph"]) if not ref.startswith("https://")}
    assert local <= set(graph)
    decision = graph["#decision-hydrogens"]
    assert (
        decision["description"] == "the deposited entry has none, the prepared file has them added"
    )
    classes = {item["name"]: item["value"] for item in decision["additionalProperty"]}
    assert classes["uncertainty_class"] == "interpretive"
    assert classes["evidence"] == "PDBbind v2020 protein files"
    assert sorted(decision["value"]) == ["Exclude", "ExplicitOnly"]
    measured = {item["name"]: item["value"] for item in graph["#findings"]["variableMeasured"]}
    assert measured["metric"] == "absolute-error"
    assert measured["universes"] == 4
    assert measured["main_effect.hydrogens.share"] == pytest.approx(4.0 / 6.0)
    assert measured["interaction.hydrogens.assembly"] == pytest.approx(1.0 / 6.0)
    assert measured["shapley.hydrogens"] == pytest.approx(0.75)
    assert graph["#algorithm"]["description"] is not None
    # A structure read from memory carries no digest, and the certificate says so.
    assert measured["inputs_without_sha256"] >= 1
    assert (
        result.certificate
        == audit.run(
            space.plan(),
            contact_count(dimer, 2.0),
            metric="absolute",
            project=lambda analysis: len(analysis.value),
        ).certificate
    )
