"""Cross-set searches preserve local row identities without within-set pairs."""

import numpy as np
import pytest

import molframe


def test_cross_pairs_match_brute_force_with_local_indices():
    first = np.array([[0, 0, 0], [1, 0, 0], [10, 0, 0]], dtype=np.float32)
    second = np.array([[0.5, 0, 0], [20, 0, 0]], dtype=np.float32)
    left, right, distance = molframe.spatial.cross_pairs(first, second, 1.0)
    assert left.tolist() == [0, 1]
    assert right.tolist() == [0, 0]
    np.testing.assert_allclose(distance, [0.5, 0.5])
    empty = np.empty((0, 3), dtype=np.float32)
    assert len(molframe.spatial.cross_pairs(empty, second, 1.0)[0]) == 0
    with pytest.raises(ValueError, match="cutoff"):
        molframe.spatial.cross_pairs(first, second, -1)
