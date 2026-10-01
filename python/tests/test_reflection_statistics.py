import math
import unittest

import numpy as np

from molframe import crystal


def lattice(limit):
    return np.array(
        [
            (h, k, l)
            for h in range(-limit, limit + 1)
            for k in range(-limit, limit + 1)
            for l in range(-limit, limit + 1)
            if (h, k, l) != (0, 0, 0)
        ],
        dtype=np.int32,
    )


def p1():
    return crystal.SpaceGroup(1)


class ReflectionStatisticsTests(unittest.TestCase):
    def setUp(self):
        self.cell = crystal.UnitCell([20.0, 22.0, 24.0], [90.0, 90.0, 90.0])
        self.hkl = lattice(5)

    def test_shells_cover_every_reflection_and_end_open(self):
        for method in ("equal_count", "dstar", "dstar2", "dstar3"):
            bins = crystal.ResolutionBins(self.cell, self.hkl, bins=6, method=method)
            self.assertEqual(len(bins), 6)
            self.assertTrue(math.isinf(bins.limits[-1]))
            spacings = np.array(
                [self.cell.reciprocal_spacing_squared(tuple(map(int, row))) for row in self.hkl]
            )
            indices = bins.indices(spacings)
            self.assertEqual(indices.min(), 0)
            self.assertEqual(indices.max(), 5)
            self.assertGreater(bins.d_max(0), bins.d_min(5))

    def test_constant_amplitudes_normalize_to_one_in_p1(self):
        bins = crystal.ResolutionBins(self.cell, self.hkl, bins=5)
        amplitudes = np.full(len(self.hkl), 7.0)
        multipliers = crystal.normalizers(self.cell, p1(), self.hkl, amplitudes, bins)
        np.testing.assert_allclose(multipliers * 7.0, 1.0, rtol=1e-12)

    def test_missing_amplitudes_give_nan(self):
        bins = crystal.ResolutionBins(self.cell, self.hkl, bins=5)
        amplitudes = np.full(len(self.hkl), 2.0)
        amplitudes[3] = np.nan
        multipliers = crystal.normalizers(self.cell, p1(), self.hkl, amplitudes, bins)
        self.assertTrue(math.isnan(multipliers[3]))
        self.assertFalse(np.isnan(np.delete(multipliers, 3)).any())

    def test_bad_input_raises_value_error(self):
        with self.assertRaises(ValueError):
            crystal.ResolutionBins(self.cell, self.hkl, bins=0)
        with self.assertRaises(ValueError):
            crystal.ResolutionBins(self.cell, self.hkl, method="spiral")
        with self.assertRaises(ValueError):
            crystal.ResolutionBins(self.cell, np.zeros((3, 2), dtype=np.int32))
        bins = crystal.ResolutionBins(self.cell, self.hkl, bins=3)
        with self.assertRaises(ValueError):
            crystal.normalizers(self.cell, p1(), self.hkl, np.ones(2), bins)
        with self.assertRaises(ValueError):
            bins.d_min(3)


if __name__ == "__main__":
    unittest.main()
