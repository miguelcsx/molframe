"""The analysis policy names the decisions a selection would otherwise make."""

import warnings

import pytest

import molframe

CIF = b"""data_t
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0 1 GLY B N
ATOM 2 C CA GLY A 1 1 1.4 0.0 0.0 1 GLY B CA
"""


@pytest.fixture(scope="module")
def structure():
    return molframe.read(CIF, name="t.cif")


def _count(structure, query, policy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return len(structure.select(query, policy=policy))


def test_the_identifier_namespace_changes_what_a_chain_selector_means(structure):
    auth = molframe.AnalysisPolicy(identifiers="auth")
    label = molframe.AnalysisPolicy(identifiers="label")
    assert _count(structure, "chain B", auth) == 2
    assert _count(structure, "chain B", label) == 0
    assert _count(structure, "chain A", label) == 2
    assert _count(structure, "chain A", auth) == 0
    # The default policy is the auth namespace, and a query object honours it too.
    assert len(structure.select("chain B")) == 2
    query = molframe.Query("chain A")
    assert len(query.select(structure, policy=label)) == 2


def test_a_policy_is_a_named_value_with_a_stable_fingerprint():
    default = molframe.AnalysisPolicy()
    assert default.profile == "molframe-default-1.0"
    assert (default.identifiers, default.altloc) == ("auth", "conformer_consistent")
    changed = molframe.AnalysisPolicy(identifiers="label", altloc="label:A")
    assert changed.profile is None
    assert changed.altloc == "label:A"
    assert changed.fingerprint != default.fingerprint
    assert changed.fingerprint == molframe.AnalysisPolicy(identifiers="label", altloc="label:A").fingerprint
    assert changed == molframe.AnalysisPolicy(identifiers="label", altloc="label:A")


@pytest.mark.parametrize(
    "keyword",
    [
        {"identifiers": "both"},
        {"altloc": "label:"},
        {"altloc": "newest"},
        {"missing_atoms": "guess"},
        {"hydrogens": "all"},
        {"symmetry": "p1"},
        {"vdw_radii": "uff"},
        {"precision": "f16"},
    ],
)
def test_unknown_choices_are_rejected_with_the_allowed_names(keyword):
    with pytest.raises(ValueError):
        molframe.AnalysisPolicy(**keyword)
