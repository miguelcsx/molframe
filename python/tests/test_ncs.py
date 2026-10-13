"""Non-crystallographic symmetry operators, chain pairs and expansion."""

import molframe

ENTRY = b"""data_demo
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 1 A 0 0 0
2 C CA GLY A 2 2 A 3.8 0 0
3 C CA GLY B 1 1 B 10 0 0
4 C CA GLY B 2 2 B 13.8 0 0
loop_
_struct_ncs_oper.id
_struct_ncs_oper.code
_struct_ncs_oper.matrix[1][1]
_struct_ncs_oper.matrix[1][2]
_struct_ncs_oper.matrix[1][3]
_struct_ncs_oper.vector[1]
_struct_ncs_oper.matrix[2][1]
_struct_ncs_oper.matrix[2][2]
_struct_ncs_oper.matrix[2][3]
_struct_ncs_oper.vector[2]
_struct_ncs_oper.matrix[3][1]
_struct_ncs_oper.matrix[3][2]
_struct_ncs_oper.matrix[3][3]
_struct_ncs_oper.vector[3]
1 given 1 0 0 10 0 1 0 0 0 0 1 0
"""


def test_operators_expose_matrix_and_translation_as_arrays():
    structure = molframe.read(ENTRY, name="entry.cif")
    (operator,) = molframe.crystal.ncs_operators(structure)
    assert operator.id == "1"
    assert operator.code == "given"
    assert operator.rotation.shape == (3, 3)
    assert operator.translation.shape == (3,)


def test_a_given_operator_binds_the_chain_it_maps_onto():
    structure = molframe.read(ENTRY, name="entry.cif")
    (pair,) = molframe.crystal.ncs_chain_pairs(structure)
    assert (pair.operator, pair.source, pair.target) == ("1", "A", "B")
    assert pair.rmsd is not None
    assert pair.rmsd < 1e-6
    assert pair.atoms == 2
    assert not pair.declared


def test_expansion_leaves_a_structure_without_generate_operators_alone():
    structure = molframe.read(ENTRY, name="entry.cif")
    expanded = molframe.crystal.ncs_structure(structure)
    assert expanded.atom_count == structure.atom_count == 4
