"""Array parameters reject a wrong dtype or rank with one clear, typed error."""

import numpy as np
import pytest

import molframe
from molframe import compare, geometry, spatial, surface
from molframe.errors import MolframeTypeError

F32 = np.zeros((4, 3), dtype=np.float32)
F64 = np.zeros((4, 3), dtype=np.float64)
RADII32 = np.ones(4, dtype=np.float32)


def test_a_float64_array_names_both_dtypes_and_the_conversion():
    with pytest.raises(MolframeTypeError) as caught:
        spatial.neighbor_pairs(F64, 3.0)
    message = str(caught.value)
    assert "expected a float32 array" in message
    assert "received dtype float64" in message
    assert "astype(numpy.float32)" in message


def test_the_offending_parameter_is_named_on_the_exception():
    with pytest.raises(MolframeTypeError) as caught:
        spatial.neighbor_pairs(F64, 3.0)
    assert any("coordinates" in note for note in caught.value.__notes__)


def test_the_error_belongs_to_the_typed_family_and_is_a_type_error():
    assert issubclass(MolframeTypeError, TypeError)
    assert issubclass(MolframeTypeError, molframe.errors.MolframeError)


@pytest.mark.parametrize(
    "call",
    [
        pytest.param(lambda: spatial.neighbor_pairs(F64, 3.0), id="spatial"),
        pytest.param(lambda: geometry.centre_of_mass(F64), id="measures"),
        pytest.param(lambda: compare.tm_score(F64, F32), id="compare-model"),
        pytest.param(lambda: compare.tm_score(F32, F64), id="compare-reference"),
        pytest.param(lambda: surface.sasa(F64, RADII32), id="surface-coordinates"),
        pytest.param(
            lambda: surface.sasa(F32, RADII32.astype(np.float64)),
            id="surface-radii",
        ),
        pytest.param(lambda: compare.tm_scores([F32, F64]), id="readonly-list"),
        pytest.param(
            lambda: compare.weighted_rmsd(F32, F32, np.ones(4, dtype=np.float32)),
            id="float64-expected",
        ),
    ],
)
def test_every_array_parameter_raises_the_typed_dtype_error(call):
    with pytest.raises(MolframeTypeError, match=r"expected a float(32|64) array"):
        call()


def test_a_float64_expected_parameter_states_float64():
    with pytest.raises(MolframeTypeError, match="expected a float64 array"):
        compare.weighted_rmsd(F32, F32, np.ones(4, dtype=np.float32))


def test_a_list_is_not_an_array():
    with pytest.raises(MolframeTypeError, match="received list"):
        spatial.neighbor_pairs([[0.0, 0.0, 0.0]], 3.0)


def test_a_wrong_rank_names_both_ranks():
    with pytest.raises(MolframeTypeError, match=r"2-dimensional.*1 dimensions"):
        spatial.neighbor_pairs(np.zeros(3, dtype=np.float32), 3.0)


def test_a_correct_float32_array_still_works():
    first, second, _ = spatial.neighbor_pairs(F32, 3.0)
    assert len(first) == len(second) == 6


def test_converting_once_with_astype_is_enough():
    first, _, _ = spatial.neighbor_pairs(F64.astype(np.float32), 3.0)
    assert len(first) == 6
