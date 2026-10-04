"""Small molecules: MDL records, their data fields and substructure queries."""

import numpy as np
import pytest

import molframe
from molframe import chemistry

SDF = """ethanol
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
> <ID>
EtOH-1

> <NOTE>
first line
second line

$$$$
acetate
  test

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -1.1000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
M  CHG  1   4  -1
M  END
$$$$
"""


def test_every_record_of_an_sdf_file_is_read_with_its_data_fields():
    first, second = chemistry.read_sdf(SDF)
    assert (first.name, first.atom_count, first.bond_count) == ("ethanol", 3, 2)
    assert first.elements == ["C", "C", "O"]
    assert first.version == "v2000"
    assert first.properties == {"ID": "EtOH-1", "NOTE": "first line\nsecond line"}
    assert np.allclose(first.coordinates[2], [2.0, 1.2, 0.0])
    assert first.bonds.tolist() == [[0, 1, 1], [1, 2, 1]]
    assert second.bonds.tolist() == [[0, 1, 1], [1, 2, 2], [1, 3, 1]]
    assert second.formal_charges == [0, 0, 0, -1]
    assert second.properties == {}


def test_what_is_written_is_what_is_read_again():
    records = chemistry.read_sdf(SDF)
    again = chemistry.read_sdf(chemistry.write_sdf(records))
    assert len(again) == 2
    for original, copy in zip(records, again, strict=True):
        assert copy.name == original.name
        assert copy.elements == original.elements
        assert np.allclose(copy.coordinates, original.coordinates, atol=1e-4)
        assert copy.bonds.tolist() == original.bonds.tolist()
        assert copy.formal_charges == original.formal_charges
        assert copy.properties == original.properties
    block = records[1].to_mol()
    assert chemistry.read_mol(block).formal_charges == [0, 0, 0, -1]
    assert records[0].to_sdf().rstrip().endswith("$$$$")


def test_a_molecule_is_built_from_arrays_and_becomes_a_structure_with_its_bonds():
    built = chemistry.Molecule(
        ["C", "C", "O"],
        np.array([[0, 0, 0], [1.5, 0, 0], [2.0, 1.2, 0]], dtype=np.float32),
        np.array([[0, 1, 1], [1, 2, 1]], dtype=np.uint32),
        name="built",
        properties={"SOURCE": "test"},
    )
    assert built.properties == {"SOURCE": "test"}
    structure = built.to_structure()
    assert structure.atom_count == 3
    assert len(structure.bonds) == 2
    assert np.allclose(structure.coordinates, built.coordinates)
    assert chemistry.read_sdf(built.to_sdf())[0].properties == {"SOURCE": "test"}


def test_a_structure_gives_back_its_graph():
    record = chemistry.read_sdf(SDF)[0]
    graph = chemistry.molecule(record.to_structure(), name="again")
    assert graph.name == "again"
    assert graph.elements == record.elements
    assert graph.bonds.tolist() == record.bonds.tolist()


def test_a_record_that_cannot_be_read_or_built_is_refused_by_name():
    with pytest.raises(molframe.ParseError) as malformed:
        chemistry.read_mol("not a mol block")
    assert malformed.value.code == "MOLFRAME-E1201"
    with pytest.raises(molframe.MolframeValueError, match="not a chemical element"):
        chemistry.Molecule(["Xx"], np.zeros((1, 3), dtype=np.float32))
    with pytest.raises(molframe.MolframeValueError, match="one row"):
        chemistry.Molecule(["C", "C"], np.zeros((1, 3), dtype=np.float32))
    with pytest.raises(molframe.MolframeError) as endpoint:
        chemistry.Molecule(
            ["C"],
            np.zeros((1, 3), dtype=np.float32),
            np.array([[0, 5, 1]], dtype=np.uint32),
        ).to_mol()
    assert endpoint.value.code == "MOLFRAME-E3006"


def test_smarts_maps_a_pattern_onto_the_atoms_it_matches():
    ethanol, acetate = (record.to_structure() for record in chemistry.read_sdf(SDF))
    assert chemistry.smarts(ethanol, "CO") == [[1, 2]]
    assert chemistry.smarts(ethanol, "C") == [[0], [1]]
    assert chemistry.smarts(ethanol, "[#8]") == [[2]]
    # Both oxygens of acetate bond to the carboxyl carbon; one of the bonds is double.
    assert chemistry.smarts(acetate, "C=O") == [[1, 2]]
    assert chemistry.smarts(acetate, "C-O") == [[1, 3]]
    assert chemistry.smarts(acetate, "C(~O)~O") == [[1, 2, 3], [1, 3, 2]]


def test_a_pattern_that_is_not_smarts_names_its_position():
    structure = chemistry.read_sdf(SDF)[0].to_structure()
    with pytest.raises(molframe.ParseError) as bad:
        chemistry.smarts(structure, "C(")
    assert bad.value.code == "MOLFRAME-E1301"
