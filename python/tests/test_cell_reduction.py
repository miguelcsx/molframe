import unittest

from molframe import crystal


class CellReductionTests(unittest.TestCase):
    def test_a_reduced_cell_is_returned_with_the_identity_basis(self):
        reduced = crystal.reduce_cell([10.0, 12.0, 15.0], [90.0, 90.0, 90.0])
        self.assertEqual(reduced.change_of_basis, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertTrue(reduced.converged)
        self.assertEqual(tuple(round(v, 9) for v in reduced.lengths), (10.0, 12.0, 15.0))

    def test_a_skewed_cell_reduces_to_shorter_edges(self):
        reduced = crystal.reduce_cell([8.9, 12.1, 57.2], [99.8, 79.5, 90.9])
        self.assertLess(reduced.lengths[2], 57.2)
        self.assertLessEqual(reduced.lengths[0], reduced.lengths[1])
        self.assertLessEqual(reduced.lengths[1], reduced.lengths[2])

    def test_a_bad_cell_or_tolerance_raises_value_error(self):
        with self.assertRaises(ValueError):
            crystal.reduce_cell([0.0, 1.0, 1.0], [90.0, 90.0, 90.0])
        with self.assertRaises(ValueError):
            crystal.reduce_cell([1.0, 1.0, 1.0], [90.0, 90.0, 90.0], epsilon=-1.0)


if __name__ == "__main__":
    unittest.main()
