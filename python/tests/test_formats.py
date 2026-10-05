"""Structure writers through the Python contract."""

from pathlib import Path

import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-py" / "tests" / "data"


@pytest.fixture(scope="module")
def structure():
    return molframe.read(DATA / "basic.cif")


def test_text_writers_round_trip_the_atom_count(structure):
    for text in (molframe.formats.to_mmcif(structure), molframe.formats.to_pdb(structure)):
        again = molframe.read(
            text.encode(), name="again.cif" if text.startswith("data_") else "again.pdb"
        )
        assert again.atom_count == structure.atom_count


def test_bcif_bytes_round_trip_and_are_deterministic(structure):
    first = molframe.formats.to_bcif(structure)
    assert first == molframe.formats.to_bcif(structure)
    assert molframe.read(first, name="x.bcif").atom_count == structure.atom_count


def test_write_picks_the_format_from_the_extension(tmp_path, structure):
    for name in ("out.cif", "out.pdb", "out.bcif"):
        target = tmp_path / name
        molframe.formats.write(structure, target)
        assert molframe.read(target).atom_count == structure.atom_count
    forced = tmp_path / "out.dat"
    molframe.formats.write(structure, forced, format="pdb")
    assert molframe.read(forced, name="x.pdb").atom_count == structure.atom_count


def test_unknown_formats_and_missing_extensions_are_rejected(tmp_path, structure):
    with pytest.raises(molframe.MolframeError) as unknown:
        molframe.formats.write(structure, tmp_path / "out.xyz")
    assert unknown.value.code == "MOLFRAME-E1001"
    with pytest.raises(molframe.MolframeError) as unnamed:
        molframe.formats.write(structure, tmp_path / "noextension")
    assert unnamed.value.code == "MOLFRAME-E1001"
    with pytest.raises(molframe.PolicyError) as forced:
        molframe.formats.write(structure, tmp_path / "out.dat", format="nonsense")
    assert forced.value.code == "MOLFRAME-E6104"


def test_explicit_transport_identifiers_keep_author_chain_insertion_and_coordinates():
    source = (
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   1A      1.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  O   HOH     2       4.000   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    structure = molframe.read(source.encode(), format="pdb")
    encoded = molframe.formats.to_bcif(
        structure, block_id="transport", generated_connections=True, transport=True
    )
    assert encoded == molframe.formats.to_bcif(
        structure, block_id="transport", generated_connections=True, transport=True
    )
    restored = molframe.read(encoded, format="bcif")
    assert restored.coordinates.tolist() == structure.coordinates.tolist()
    assert [chain.auth_label for chain in restored.chains] == [
        chain.auth_label for chain in structure.chains
    ]
    assert [(residue.auth_number, residue.insertion_code) for residue in restored.residues] == [
        (residue.auth_number, residue.insertion_code) for residue in structure.residues
    ]
    assert len({chain.label for chain in restored.chains}) == restored.chain_count
