"""Regression tests for selection order, chain naming and chain interfaces."""

import numpy as np
import pytest

import molframe

# Label A is author H and label B is author L; a second pair of chains reuses
# the letters swapped, so the bare names "A" and "B" are ambiguous between the
# label and author namespaces.
SWAPPED = b"""data_sw
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA GLY A B 1 0 0 0
ATOM 2 C CA GLY A B 2 1 0 0
ATOM 3 C CA GLY B A 1 3 0 0
ATOM 4 C CA GLY B A 2 50 0 0
"""


@pytest.fixture(scope="module")
def swapped():
    return molframe.read(SWAPPED, name="swapped.cif")


def atoms(structure, query):
    return sorted(int(i) for i in structure.select(query).indices)


def test_conjunction_with_within_does_not_depend_on_operand_order(swapped):
    forward = atoms(swapped, "label_chain A and within 2.5 of label_chain B")
    reverse = atoms(swapped, "within 2.5 of label_chain B and label_chain A")
    assert forward == reverse == [1]
    assert atoms(swapped, "label_chain A and byres within 2.5 of label_chain B") == atoms(
        swapped, "byres within 2.5 of label_chain B and label_chain A"
    )


def test_chain_values_can_name_their_namespace(swapped):
    assert atoms(swapped, "chain label:A") == atoms(swapped, "label_chain A") == [0, 1]
    assert atoms(swapped, "chain auth:A") == atoms(swapped, "auth_chain A") == [2, 3]
    assert atoms(swapped, "chain label:A auth:A") == [0, 1, 2, 3]


def test_chain_interface_accepts_namespaced_chain_names(swapped):
    with pytest.raises(ValueError, match=r"label:A.*auth:A|auth:A.*label:A"):
        molframe.analysis.chain_interface(swapped, first_chain="A", second_chain="B", cutoff=4.0)
    by_label = molframe.analysis.chain_interface(
        swapped, first_chain="label:A", second_chain="label:B", cutoff=2.5
    )
    by_auth = molframe.analysis.chain_interface(
        swapped, first_chain="auth:B", second_chain="auth:A", cutoff=2.5
    )
    assert np.array_equal(by_label.value["residue"], by_auth.value["residue"])
    assert len(by_label.value) == 2


def test_salt_bridges_without_chemistry_is_indeterminate_not_empty(swapped):
    result = molframe.analysis.salt_bridges(swapped)
    assert result.status == "indeterminate"
    assert result.indeterminacy is not None
    assert "annotate" in result.indeterminacy
    with pytest.raises(molframe.IndeterminateError):
        result.value  # noqa: B018


_TITRATABLE = b"""data_t
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 20 20 20
ATOM 2 O OD2 ASP A 2 0 0 0
ATOM 3 N NZ LYS A 3 3 0 0
"""


def test_salt_bridges_find_protein_charges_without_a_charge_annotation():
    structure = molframe.read(_TITRATABLE, name="titratable.cif")
    result = molframe.analysis.salt_bridges(structure, max_distance=4.0)
    assert result.status == "complete"
    table = result.value
    assert [(int(a), int(c)) for a, c in zip(table["anion"], table["cation"], strict=True)] == [
        (1, 2)
    ]
    assert any("Asp/Glu" in assumption for assumption in result.assumptions)


_WATER_AND_IRON = b"""data_w
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA GLY A 1 0 0 0
HETATM 2 FE FE HEM B . 30 0 0
HETATM 3 O O HOH C . 0 3 0
"""


def test_sasa_names_the_policy_field_behind_a_partial_result():
    structure = molframe.read(_WATER_AND_IRON, name="fe.cif")
    result = molframe.analysis.sasa(structure, hetero="all", metals="none")
    assert result.status == "partial"
    assert any("W6004" in str(warning) and "Fe" in str(warning) for warning in result.warnings)
    assert any(a.startswith("VdwRadii") for a in result.assumptions)
    assert any("1 water" in a for a in result.assumptions)
    assert np.isnan(result.value[1])
