"""Interaction analyses need a component dictionary to know charges and roles."""

import numpy as np
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


def test_without_a_dictionary_the_chemistry_analyses_say_what_is_missing():
    # The hydrogen keeps the input from being indeterminate for lack of hydrogens first.
    with_hydrogen = STRUCTURE + b"ATOM 6 H HZ1 LYS A 1 2 2.9 -1.5 0.5 2 LYS A HZ1\n"
    plain = molframe.read(with_hydrogen, name="sb-h.cif")
    with pytest.raises(ValueError, match="CCD donor/acceptor annotations"):
        molframe.analysis.hydrogen_bonds(plain)


def test_annotating_keeps_the_atoms_and_their_coordinates(annotated):
    plain, chemical = annotated
    assert chemical.atom_count == plain.atom_count
    assert np.array_equal(chemical.coordinates, plain.coordinates)


def test_a_salt_bridge_is_found_between_the_charged_atoms(annotated):
    _, chemical = annotated
    result = molframe.analysis.salt_bridges(chemical, max_distance=4.0)
    table = result.value
    assert len(table) == 1
    assert (int(table["anion"][0]), int(table["cation"][0])) == (2, 3)
    assert table["distance"][0] == pytest.approx(np.hypot(1.8, 0.1), abs=1e-3)
    assert not table["distance"].flags.writeable
    assert set(table.names) == {"anion", "cation", "distance"}
    assert "anion" in table
    assert "missing" not in table
    with pytest.raises(KeyError):
        table["missing"]
    # Outside the cutoff there is nothing to report, and that is a complete answer.
    far = molframe.analysis.salt_bridges(chemical, max_distance=1.0)
    assert len(far.value) == 0
    assert far.status == "complete"


def test_a_hydrogen_bond_geometry_without_hydrogens_has_no_answer_not_an_empty_one(annotated):
    _, chemical = annotated
    result = molframe.analysis.hydrogen_bonds(chemical)
    # An empty table would say "there are no hydrogen bonds"; the input cannot say that.
    assert result.status == "indeterminate"
    assert result.indeterminacy is not None
    assert "no hydrogen" in result.indeterminacy
    with pytest.raises(molframe.IndeterminateError):
        result.value  # noqa: B018


def test_excluding_hydrogens_is_refused_for_an_analysis_that_measures_them(annotated):
    _, chemical = annotated
    with pytest.raises(molframe.PolicyError) as refused:
        molframe.analysis.hydrogen_bonds(
            chemical, policy=molframe.AnalysisPolicy(hydrogens="exclude")
        )
    assert refused.value.code == "MOLFRAME-E6103"
    assert "hydrogen" in str(refused.value)


def test_an_analysis_says_what_it_estimates(annotated):
    _, chemical = annotated
    estimand = molframe.analysis.salt_bridges(chemical).estimand
    assert estimand is not None
    assert "charge" in estimand
