import math
from pathlib import Path

import pytest

import molframe


def _structure():
    return molframe.read(
        b"HETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        b"HETATM    2  O   LIG A   1       1.400   0.000   0.000  1.00  0.00           O\nEND\n",
        name="ligand.pdb",
    )


def test_complete_file_charges_do_not_require_a_dictionary():
    structure = molframe.read(
        b"ATOM      1 C    LIG A   1       0.000   0.000   0.000  0.1250 1.7000\n",
        name="charged.pqr",
    )
    charges = molframe.chemistry.partial_charges(structure)
    assert charges.values == [0.125]
    assert charges.source == "file"
    assert charges.dictionary_version is None
    assert charges.parameter_profile is None


def test_ccd_charge_projection_and_provenance_are_public(tmp_path):
    dictionary = tmp_path / "ligand.cif"
    dictionary.write_text(
        "data_LIG\n_chem_comp.id LIG\n_chem_comp.name methanol\n"
        "_chem_comp.type NON-POLYMER\nloop_\n_chem_comp_atom.atom_id\n"
        "_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n"
        "_chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_leaving_atom_flag\n"
        "_chem_comp_atom.pdbx_stereo_config\nC C 0 N N N\nO O 0 N N N\nH H 0 N N N\n"
        "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n"
        "_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n"
        "_chem_comp_bond.pdbx_stereo_config\nC O SING N N\nO H SING N N\n"
    )
    structure = _structure()
    charges = molframe.chemistry.partial_charges(structure, dictionary, version="fixture")
    assert charges.source == "peoe"
    assert charges.dictionary_version == "fixture"
    assert charges.parameter_profile == "gasteiger-marsili"
    assert charges.values[0] > 0
    assert charges.values[1] < 0
    assert sum(charges.values) == pytest.approx(0, abs=1e-12)
    with pytest.raises(AttributeError):
        charges.source = "file"
    with pytest.raises(ValueError, match="UnknownComponent"):
        molframe.chemistry.partial_charges(structure)


@pytest.mark.parametrize("radius", [5.0, 10.0])
def test_public_contact_potential_matches_a_unit_charge(radius):
    grid = molframe.analysis.GridSpec(
        [[1, 0, 0, radius], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], [1, 1, 1]
    )
    field = molframe.analysis.contact_potential(_structure(), [1.0, 0.0], grid)
    expected = 1.602176634e-19**2 / (
        4 * math.pi * 8.8541878128e-12 * 1e-10 * 1.380649e-23 * 298 * 4 * radius**2
    )
    assert field.values[0] == pytest.approx(expected, rel=1e-4)
    assert field.spec.dimensions == [1, 1, 1]
    assert field.spec.voxel_to_world == grid.voxel_to_world
    with pytest.raises(ValueError, match="dense coordinate and charge"):
        molframe.analysis.contact_potential(_structure(), [1.0], grid)


def test_real_4hhb_assembly_has_no_cross_instance_covalent_links():
    path = Path(__file__).resolve().parents[2] / "crates/molframe-bench/data/4hhb.cif"
    structure = molframe.read(path)
    assert molframe.crystal.assembly_covalent_links(structure, "1") == []
    with pytest.raises(ValueError, match="E6002"):
        molframe.crystal.assembly_covalent_links(structure, "missing")
    with pytest.raises(ValueError, match="E6003"):
        molframe.crystal.assembly_covalent_links(structure, "1", model=1)
