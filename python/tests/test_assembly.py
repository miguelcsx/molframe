"""Biological assemblies as placements of the deposited structure."""

import pytest

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
1 C CA GLY A 1 1 A 1 0 0
2 C CA GLY B 1 1 B 0 2 0
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
I 1 0 0 0 0 1 0 0 0 0 1 0
Z 0 -1 0 5 1 0 0 0 0 0 1 0
_pdbx_struct_assembly.id 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 I A,B
1 Z A
"""


def test_an_assembly_is_a_list_of_placements_with_their_chains():
    structure = molframe.read(ENTRY, name="entry.cif")
    assert molframe.crystal.assemblies(structure) == ["1"]
    identity, rotated = molframe.crystal.assembly(structure, "1")
    assert identity.chains == ["A", "B"]
    assert identity.matrix == [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    assert rotated.chains == ["A"]
    # column-major: a quarter turn about z, then the 5 Å translation in x
    assert rotated.matrix == pytest.approx([0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 5, 0, 0, 1])


def test_unknown_assemblies_and_structures_without_any_are_errors():
    structure = molframe.read(ENTRY, name="entry.cif")
    with pytest.raises(ValueError, match="E6002"):
        molframe.crystal.assembly(structure, "9")
    plain = molframe.read(
        b"ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 10.00           N\nEND\n",
        name="p.pdb",
    )
    assert molframe.crystal.assemblies(plain) == []
    with pytest.raises(ValueError, match="the structure declares no biological assemblies"):
        molframe.crystal.assembly(plain, "1")


def test_an_assembly_expands_to_one_structure_with_every_copy_placed():
    structure = molframe.read(ENTRY, name="entry.cif")
    whole = molframe.crystal.assembly_structure(structure, "1")
    assert whole.atom_count == 3
    assert whole.chain_count == 3
    # the copy of A: a quarter turn about z of (1, 0, 0), then +5 in x
    assert whole.coordinates[2].tolist() == pytest.approx([5, 1, 0])


def test_streamed_pieces_concatenate_to_the_whole_expansion():
    structure = molframe.read(ENTRY, name="entry.cif")
    whole = molframe.crystal.assembly_structure(structure, "1")
    for bound in (1, 2, 1000):
        pieces = list(molframe.crystal.assembly_chunks(structure, "1", bound))
        assert sum(piece.atom_count for piece in pieces) == whole.atom_count
        merged = [row for piece in pieces for row in piece.coordinates.tolist()]
        assert merged == whole.coordinates.tolist()
    # placements are never split: the first placement holds both source chains
    sizes = [p.atom_count for p in molframe.crystal.assembly_chunks(structure, "1", 1)]
    assert sizes == [2, 1]


def test_expansion_refuses_unknown_assemblies():
    structure = molframe.read(ENTRY, name="entry.cif")
    with pytest.raises(ValueError, match="E6002"):
        molframe.crystal.assembly_structure(structure, "9")
    with pytest.raises(ValueError, match="E6002"):
        molframe.crystal.assembly_chunks(structure, "9")
