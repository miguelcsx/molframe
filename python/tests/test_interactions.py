"""Interaction analyses need a component dictionary to know charges and roles."""

import numpy
import pytest

import molframe

CCD = """data_ASP
_chem_comp.id ASP
_chem_comp.name 'ASPARTIC ACID'
_chem_comp.type 'L-peptide linking'
_chem_comp.one_letter_code D
_chem_comp.formula 'C4 H7 N O4'
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
CG C 0 N N N
OD1 O 0 N N N
OD2 O -1 N N N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
CG OD1 DOUB N N
CG OD2 SING N N
#
data_LYS
_chem_comp.id LYS
_chem_comp.name LYSINE
_chem_comp.type 'L-peptide linking'
_chem_comp.one_letter_code K
_chem_comp.formula 'C6 H15 N2 O2'
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
CE C 0 N N N
NZ N 1 N N N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
CE NZ SING N N
#
"""

STRUCTURE = b"""data_sb
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
ATOM 1 C CG ASP A 1 1 0.0 0.0 0.0 1 ASP A CG
ATOM 2 O OD1 ASP A 1 1 0.6 1.1 0.0 1 ASP A OD1
ATOM 3 O OD2 ASP A 1 1 0.6 -1.1 0.0 1 ASP A OD2
ATOM 4 N NZ LYS A 1 2 2.4 -1.0 0.0 2 LYS A NZ
ATOM 5 C CE LYS A 1 2 3.5 -1.2 0.5 2 LYS A CE
"""


@pytest.fixture(scope="module")
def annotated(tmp_path_factory):
    path = tmp_path_factory.mktemp("ccd") / "components.cif"
    path.write_text(CCD)
    structure = molframe.read(STRUCTURE, name="sb.cif")
    return structure, molframe.chemistry.annotate(structure, path, version="test-1")


def test_without_a_dictionary_the_chemistry_analyses_say_what_is_missing(annotated):
    plain, _ = annotated
    with pytest.raises(ValueError, match="MissingChemistry"):
        molframe.analysis.hydrogen_bonds(plain)


def test_annotating_keeps_the_atoms_and_their_coordinates(annotated):
    plain, chemical = annotated
    assert chemical.atom_count == plain.atom_count
    assert numpy.array_equal(chemical.coordinates, plain.coordinates)


def test_a_salt_bridge_is_found_between_the_charged_atoms(annotated):
    _, chemical = annotated
    result = molframe.analysis.salt_bridges(chemical, max_distance=4.0)
    table = result.value
    assert len(table) == 1
    assert (int(table["anion"][0]), int(table["cation"][0])) == (2, 3)
    assert table["distance"][0] == pytest.approx(numpy.hypot(1.8, 0.1), abs=1e-3)
    assert not table["distance"].flags.writeable
    assert set(table.names) == {"anion", "cation", "distance"}
    assert "anion" in table and "missing" not in table
    with pytest.raises(KeyError):
        table["missing"]
    # Outside the cutoff there is nothing to report, and that is a complete answer.
    far = molframe.analysis.salt_bridges(chemical, max_distance=1.0)
    assert len(far.value) == 0 and far.status == "complete"
