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
        again = molframe.read(text.encode(), name="again.cif" if text.startswith("data_") else "again.pdb")
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
    with pytest.raises(ValueError):
        molframe.formats.write(structure, tmp_path / "out.xyz")
    with pytest.raises(ValueError):
        molframe.formats.write(structure, tmp_path / "noextension")
