"""Explicit polymer links refine component chemistry through the native policy."""

from pathlib import Path

import numpy as np
import pytest

import molframe

ROOT = Path(__file__).resolve().parents[2]
CCD = ROOT / "crates/molframe-chem/data/CCD-amino-acids.cif"


def test_explicit_peptide_policy_keeps_correct_roles_without_copying_coordinates(tmp_path):
    path = tmp_path / "two-residues.cif"
    path.write_text(
        "data_peptide\nloop_\n_atom_site.id\n_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "1 N N GLY A 1 0 0 0\n2 C CA GLY A 1 1.4 0 0\n3 C C GLY A 1 2.8 0 0\n"
        "4 N N GLY A 2 4.1 0 0\n5 C CA GLY A 2 5.5 0 0\n6 C C GLY A 2 6.9 0 0\n"
    )
    structure = molframe.read(path)
    rule = molframe.chemistry.PolymerLinkRule(1, "C", 1, "N")
    policy = molframe.chemistry.PolymerLinkPolicy(1.6, [rule])
    assert (policy.angstrom, policy.rules[0].left_atom, policy.rules[0].right_atom) == (
        pytest.approx(1.6),
        "C",
        "N",
    )
    assert (rule.left_kind, rule.right_kind) == (1, 1)
    with pytest.warns(UserWarning, match="component finding"):
        isolated = molframe.chemistry.annotate(structure, CCD)
    with pytest.warns(UserWarning, match="component finding"):
        linked = molframe.chemistry.annotate(structure, CCD, polymer_links=policy)
    nitrogens = [index for index, atom in enumerate(linked.atoms) if atom.name == "N"]
    isolated_roles = np.asarray(isolated.annotations["hbond_acceptor"].values)
    linked_roles = np.asarray(linked.annotations["hbond_acceptor"].values)
    assert len(nitrogens) == 2
    # The facade read already infers this short peptide bond; disabling new
    # polymer links preserves that connectivity and its refined roles.
    assert np.array_equal(isolated_roles[nitrogens], linked_roles[nitrogens])
    assert linked_roles[nitrogens[0]]
    assert not any(linked_roles[nitrogens[1:]])
    assert np.shares_memory(linked.coordinates, structure.coordinates)


def test_deposited_protein_links_refine_roles_without_an_inference_policy():
    structure = molframe.read(ROOT / "crates/molframe-bench/data/1ubq.cif")
    with pytest.warns(UserWarning, match="component finding"):
        annotated = molframe.chemistry.annotate(structure, CCD)
    nitrogens = [index for index, atom in enumerate(annotated.atoms) if atom.name == "N"]
    roles = np.asarray(annotated.annotations["hbond_acceptor"].values)
    assert len(nitrogens) == 76
    assert roles[nitrogens[0]]
    assert not any(roles[nitrogens[1:]])


@pytest.mark.parametrize(("distance", "rules"), [(0.0, True), (float("nan"), True), (1.6, False)])
def test_invalid_link_policies_reach_the_native_validator(distance, rules):
    structure = molframe.read(ROOT / "crates/molframe-bench/data/1ubq.cif")
    attachments = [molframe.chemistry.PolymerLinkRule(1, "C", 1, "N")] if rules else []
    policy = molframe.chemistry.PolymerLinkPolicy(distance, attachments)
    with pytest.raises(molframe.ConversionError, match="polymer-link"):
        molframe.chemistry.annotate(structure, CCD, polymer_links=policy)


def test_attachment_rules_reject_unknown_component_kind_codes():
    with pytest.raises(ValueError, match="component kind"):
        molframe.chemistry.PolymerLinkRule(999, "C", 1, "N")
