"""Motion candidates preserve original atom addresses and explicit chemistry."""

from pathlib import Path

import numpy as np
import pytest

import molframe

CCD = Path(__file__).resolve().parents[2] / "crates/molframe-chem/data/CCD-amino-acids.cif"
PDB = b"""ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  H   ALA A   1       1.000   0.000   0.000  1.00 20.00           H
ATOM      3  O   ALA A   2       6.000   0.000   0.000  1.00 20.00           O
END
"""


@pytest.fixture
def chemical() -> molframe.Structure:
    """Use CCD roles and a real donor-hydrogen covalent edge."""
    source = molframe.read(PDB, format="pdb")
    with pytest.warns(UserWarning, match="component finding"):
        return molframe.chemistry.annotate(source, CCD)


def test_candidates_include_a_contact_inside_the_segment_only(
    chemical: molframe.Structure,
) -> None:
    """Both endpoints are outside the distance cutoff; the middle crosses it."""
    start = chemical.coordinates
    end = np.array(start, copy=True)
    end[2, 0] = -6.0
    table = molframe.analysis.hydrogen_bond_candidates(chemical, start, end)
    assert len(table) == 1
    assert tuple(int(table[name][0]) for name in ("donor", "hydrogen", "acceptor")) == (
        0,
        1,
        2,
    )
    assert not table["donor"].flags.writeable
    np.testing.assert_array_equal(start, chemical.coordinates)


def test_common_translation_keeps_an_outside_pair_outside(
    chemical: molframe.Structure,
) -> None:
    """Display-independent Cartesian motion does not invent a contact."""
    start = chemical.coordinates
    end = np.array(start, copy=True) + 1000.0
    assert len(molframe.analysis.hydrogen_bond_candidates(chemical, start, end)) == 0


def test_missing_roles_and_malformed_coordinate_columns_are_rejected(
    chemical: molframe.Structure,
) -> None:
    """No heavy-atom or malformed-buffer fallback enters the analysis."""
    source = molframe.read(PDB, format="pdb")
    with pytest.raises(ValueError, match="CCD donor/acceptor"):
        molframe.analysis.hydrogen_bond_candidates(source, source.coordinates, source.coordinates)
    with pytest.raises(ValueError, match="shape"):
        molframe.analysis.hydrogen_bond_candidates(
            chemical, np.zeros((3, 2), dtype=np.float32), chemical.coordinates
        )
    end = np.array(chemical.coordinates, copy=True)
    end[2, 0] = np.nan
    with pytest.raises(ValueError, match="geometry options"):
        molframe.analysis.hydrogen_bond_candidates(chemical, chemical.coordinates, end)


@pytest.mark.parametrize("distance", [0.0, -1.0, float("inf")])
def test_invalid_distance_envelopes_reach_the_native_validator(distance: float) -> None:
    """Python construction shares native candidate policy validation."""
    with pytest.raises(ValueError, match="geometry options"):
        molframe.analysis.HydrogenBondCandidateOptions(max_distance=distance)
