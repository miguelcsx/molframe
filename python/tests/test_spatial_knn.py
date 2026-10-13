"""Parallel k-nearest-neighbour search and chunked pair lists against scipy."""

import numpy as np
import pytest

import molframe

cKDTree = pytest.importorskip("scipy.spatial").cKDTree  # noqa: N816

NONE = np.iinfo(np.uint32).max
BACKENDS = ["auto", "cell", "kd_tree"]


def layouts():
    rng = np.random.default_rng(20260101)
    base = rng.random((4000, 3))
    clustered = np.concatenate([rng.normal(0, 0.4, (2500, 3)), rng.random((1500, 3)) * 80.0])
    outliers = rng.random((3000, 3)) * 30.0
    outliers[:3] = [[1e6, 0, 0], [0, -1e6, 0], [0, 0, 1e6]]
    duplicates = np.repeat(rng.random((400, 3)) * 20.0, 5, axis=0)
    line = np.zeros((3000, 3))
    line[:, 0] = rng.random(3000) * 300.0
    return {
        "uniform": base * 40.0,
        "needle": base * np.array([1.0, 1.0, 600.0]),
        "slab": base * np.array([60.0, 60.0, 2.0]),
        "clustered": clustered,
        "offset": base * 40.0 + 1.0e5,
        "collinear": line,
        "outliers": outliers,
        "duplicates": duplicates,
    }


LAYOUTS = layouts()


def f32(points):
    return np.ascontiguousarray(points, dtype=np.float32)


def near_boundary(points, pairs, cutoff):
    """Pairs whose float64 distance sits at the float32 rounding boundary."""
    if len(pairs) == 0:
        return set()
    gap = np.linalg.norm(points[pairs[:, 0]] - points[pairs[:, 1]], axis=1)
    return {tuple(p) for p in pairs[np.abs(gap - cutoff) <= 1e-4 * max(cutoff, 1.0)]}


@pytest.mark.parametrize("name", sorted(LAYOUTS))
@pytest.mark.parametrize("backend", BACKENDS)
def test_neighbor_pairs_match_scipy(name, backend):
    points = f32(LAYOUTS[name])
    wide = points.astype(np.float64)
    cutoff = 3.0
    first, second, distance = molframe.spatial.neighbor_pairs(points, cutoff, backend=backend)
    ours = {(int(a), int(b)) for a, b in zip(first, second, strict=True)}
    theirs = {tuple(p) for p in cKDTree(wide).query_pairs(cutoff, output_type="ndarray")}
    assert len(ours) == len(first)
    difference = ours ^ theirs
    allowed = near_boundary(wide, np.array(sorted(difference)).reshape(-1, 2), cutoff)
    assert difference <= allowed
    assert np.all(distance <= cutoff * (1 + 1e-5))
    assert np.all(first < second)


@pytest.mark.parametrize("name", sorted(LAYOUTS))
@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("k", [1, 8, 16])
def test_knn_matches_scipy(name, backend, k):
    points = f32(LAYOUTS[name])
    wide = points.astype(np.float64)
    indices, distances = molframe.spatial.knn(points, k, backend=backend)
    assert indices.shape == distances.shape == (len(points), k)
    assert indices.dtype == np.uint32
    assert distances.dtype == np.float32
    expected = cKDTree(wide).query(wide, k + 1)[0][:, 1:]
    np.testing.assert_allclose(distances, expected, rtol=2e-5, atol=2e-5)
    own = np.arange(len(points))[:, None]
    assert not np.any(indices == own)
    # Rows are ordered by distance, and the distance of every returned index is real.
    assert np.all(np.diff(distances, axis=1) >= 0)
    check = np.linalg.norm(wide[indices.astype(np.int64)] - wide[:, None, :], axis=2)
    np.testing.assert_allclose(distances, check, rtol=2e-5, atol=2e-5)


def test_knn_is_independent_of_backend_and_worker_count():
    points = f32(LAYOUTS["clustered"])
    reference = molframe.spatial.knn(points, 12, backend="brute_force")
    for backend in BACKENDS:
        for workers in (1, 4):
            context = molframe.ExecutionContext(workers=workers)
            got = molframe.spatial.knn(points, 12, backend=backend, context=context)
            np.testing.assert_array_equal(got[0], reference[0])
            np.testing.assert_array_equal(got[1], reference[1])


def test_knn_query_set_max_distance_and_padding():
    points = f32(LAYOUTS["uniform"])
    query = f32(LAYOUTS["uniform"][:50] + 0.01)
    indices, distances = molframe.spatial.knn(points, 6, query=query, max_distance=1.5)
    assert indices.shape == (50, 6)
    tree = cKDTree(points.astype(np.float64))
    expected, _ = tree.query(query.astype(np.float64), 6, distance_upper_bound=1.5 * (1 + 1e-6))
    missing = ~np.isfinite(expected)
    assert np.array_equal(indices == NONE, missing)
    assert np.all(np.isinf(distances[missing]))
    np.testing.assert_allclose(distances[~missing], expected[~missing], rtol=2e-5, atol=2e-5)
    short = molframe.spatial.knn(points[:4], 6)
    assert np.all(short[0][:, 3:] == NONE)
    assert np.all(np.isinf(short[1][:, 3:]))
    with pytest.raises(ValueError, match="finite"):
        molframe.spatial.knn(points, 4, max_distance=float("nan"))


@pytest.mark.parametrize("angles", [(90.0, 90.0, 90.0), (70.0, 80.0, 65.0)])
def test_periodic_knn_is_minimum_image(angles):
    rng = np.random.default_rng(5)
    lengths = (18.0, 21.0, 24.0)
    cell = molframe.crystal.UnitCell(list(lengths), list(angles))
    fractional = rng.random((1500, 3))
    matrix = np.array(cell.matrix if hasattr(cell, "matrix") else _matrix(lengths, angles))
    points = f32(fractional @ matrix)
    indices, distances = molframe.spatial.knn(
        points, 8, cell=cell, max_distance=5.0, backend="kd_tree"
    )
    brute = molframe.spatial.knn(points, 8, cell=cell, max_distance=5.0, backend="brute_force")
    np.testing.assert_array_equal(indices, brute[0])
    np.testing.assert_array_equal(distances, brute[1])
    assert np.isfinite(distances).any()


def _matrix(lengths, angles):
    a, b, c = lengths
    alpha, beta, gamma = np.radians(angles)
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c * c - cx * cx - cy * cy)
    return [[a, 0, 0], [b * np.cos(gamma), b * np.sin(gamma), 0], [cx, cy, cz]]


@pytest.mark.parametrize("name", ["uniform", "clustered", "duplicates"])
def test_neighbor_pairs_chunks_concatenate_to_neighbor_pairs(name):
    points = f32(LAYOUTS[name])
    whole = molframe.spatial.neighbor_pairs(points, 2.5)
    chunks = list(molframe.spatial.neighbor_pairs_chunks(points, 2.5, chunk_size=777))
    assert all(len(chunk[0]) <= 777 for chunk in chunks)
    for column in range(3):
        np.testing.assert_array_equal(np.concatenate([c[column] for c in chunks]), whole[column])
    with pytest.raises(ValueError, match="positive"):
        molframe.spatial.neighbor_pairs_chunks(points, 2.5, chunk_size=0)


def test_pair_columns_own_their_memory_without_copies():
    first, second, distance = molframe.spatial.neighbor_pairs(f32(LAYOUTS["uniform"]), 2.0)
    for column in (first, second, distance):
        assert column.flags.c_contiguous
        assert column.flags.owndata is False
