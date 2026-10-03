"""Direct X-ray structure factors against values computed with Gemmi 0.7.5."""

import numpy as np
import pytest

import molframe

# A tiny P 21 21 21 model with one anisotropic atom and one half-occupied atom.
PDB = b"""CRYST1   12.000   14.000   16.000  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   ALA A   1       1.500   2.250   3.125  1.00 12.00           N
ANISOU    1  N   ALA A   1      180    210    170     20    -10     15       N
ATOM      2  CA  ALA A   1       2.750   2.900   3.800  1.00 15.00           C
ATOM      3  C   ALA A   1       3.900   1.950   4.100  0.70 14.00           C
ATOM      4  O   ALA A   1       4.100   0.800   3.900  1.00 18.00           O
ATOM      5  S   CYS A   2       5.250   6.500   7.750  1.00 20.00           S
END
"""

# REFERENCE was produced by gemmi: StructureFactorCalculatorX(cell), then
# calculate_sf_from_model(model, hkl).
REFERENCE = {
    (1, 0, 0): 0.0,  # a screw-axis absence
    (2, 0, 0): -1.2190625415629484 + 0j,
    (0, 1, 3): -31.396970818024215j,
    (3, 2, 1): -17.27280307252613 + 20.648885613620507j,
    (-2, 5, 4): -14.170811514306425 + 3.3211426639433776j,
    (4, -1, 2): -9.13728370946391 - 6.471406855570705j,
}


@pytest.fixture(scope="module")
def structure():
    return molframe.read(PDB, name="fixture.pdb")


def test_values_agree_with_the_reference_implementation(structure):
    hkl = np.array(list(REFERENCE), dtype=np.int32)
    computed = molframe.crystal.structure_factors(structure, hkl)
    assert computed.dtype == np.complex128
    expected = np.array(list(REFERENCE.values()), dtype=np.complex128)
    # Gemmi keeps its scattering coefficients in single precision.
    assert np.allclose(computed, expected, rtol=1e-4, atol=1e-4)


def test_friedel_mates_have_equal_amplitudes_and_opposite_phases(structure):
    hkl = np.array([[3, 2, 1], [-3, -2, -1]], dtype=np.int32)
    plus, minus = molframe.crystal.structure_factors(structure, hkl)
    assert abs(plus) == pytest.approx(abs(minus), rel=1e-12)
    assert plus.conjugate() == pytest.approx(minus, abs=1e-9)


def test_inputs_without_a_cell_or_with_a_bad_shape_are_rejected(structure):
    with pytest.raises(ValueError, match="hkl must have shape"):
        molframe.crystal.structure_factors(structure, np.zeros((3, 2), dtype=np.int32))
    no_cell = molframe.read(
        b"ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 10.00           N\nEND\n",
        name="n.pdb",
    )
    with pytest.raises(ValueError, match="E5004"):
        molframe.crystal.structure_factors(no_cell, np.array([[1, 0, 0]], dtype=np.int32))
