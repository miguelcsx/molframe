"""Resolution bins and reflection normalisation."""

import math

import numpy as np
import pytest

from molframe import crystal


def _lattice(limit):
    return np.array(
        [
            (h, k, i)
            for h in range(-limit, limit + 1)
            for k in range(-limit, limit + 1)
            for i in range(-limit, limit + 1)
            if (h, k, i) != (0, 0, 0)
        ],
        dtype=np.int32,
    )


def _p1():
    return crystal.SpaceGroup(1)


@pytest.fixture
def cell():
    return crystal.UnitCell([20.0, 22.0, 24.0], [90.0, 90.0, 90.0])


@pytest.fixture
def hkl():
    return _lattice(5)


class TestReflectionStatistics:
    def test_shells_cover_every_reflection_and_end_open(self, cell, hkl):
        for method in ("equal_count", "dstar", "dstar2", "dstar3"):
            bins = crystal.ResolutionBins(cell, hkl, bins=6, method=method)
            assert len(bins) == 6
            assert math.isinf(bins.limits[-1])
            spacings = np.array(
                [cell.reciprocal_spacing_squared(tuple(map(int, row))) for row in hkl]
            )
            indices = bins.indices(spacings)
            assert indices.min() == 0
            assert indices.max() == 5
            assert bins.d_max(0) > bins.d_min(5)

    def test_constant_amplitudes_normalize_to_one_in_p1(self, cell, hkl):
        bins = crystal.ResolutionBins(cell, hkl, bins=5)
        amplitudes = np.full(len(hkl), 7.0)
        multipliers = crystal.normalizers(cell, _p1(), hkl, amplitudes, bins)
        np.testing.assert_allclose(multipliers * 7.0, 1.0, rtol=1e-12)

    def test_missing_amplitudes_give_nan(self, cell, hkl):
        bins = crystal.ResolutionBins(cell, hkl, bins=5)
        amplitudes = np.full(len(hkl), 2.0)
        amplitudes[3] = np.nan
        multipliers = crystal.normalizers(cell, _p1(), hkl, amplitudes, bins)
        assert math.isnan(multipliers[3])
        assert not np.isnan(np.delete(multipliers, 3)).any()

    def test_bad_input_raises_value_error(self, cell, hkl):
        with pytest.raises(ValueError, match="the number of resolution bins must be positive"):
            crystal.ResolutionBins(cell, hkl, bins=0)
        with pytest.raises(
            ValueError, match="method must be one of equal_count, dstar, dstar2, dstar3"
        ):
            crystal.ResolutionBins(cell, hkl, method="spiral")
        with pytest.raises(ValueError, match="hkl must have shape"):
            crystal.ResolutionBins(cell, np.zeros((3, 2), dtype=np.int32))
        bins = crystal.ResolutionBins(cell, hkl, bins=3)
        with pytest.raises(ValueError, match="reflection columns have different lengths"):
            crystal.normalizers(cell, _p1(), hkl, np.ones(2), bins)
        with pytest.raises(ValueError, match="shell index out of range"):
            bins.d_min(3)
