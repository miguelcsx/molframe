"""Ensemble analyses of a trajectory: each result is checked against NumPy."""

import math

import numpy as np
import pytest

import molframe
from molframe import trajectory


def kabsch_rmsd(first: np.ndarray, second: np.ndarray) -> float:
    """RMSD after the best rigid fit, from the SVD of the covariance."""
    a = first.astype(np.float64) - first.astype(np.float64).mean(axis=0)
    b = second.astype(np.float64) - second.astype(np.float64).mean(axis=0)
    u, singular, vt = np.linalg.svd(a.T @ b)
    sign = np.sign(np.linalg.det(u @ vt))
    flip = np.array([1.0, 1.0, sign])
    sum_squares = (a**2).sum() + (b**2).sum() - 2.0 * (singular * flip).sum()
    return math.sqrt(max(sum_squares, 0.0) / len(a))


@pytest.fixture(scope="module")
def frames():
    return np.random.default_rng(7).normal(size=(12, 6, 3)).astype(np.float32)


def test_pairwise_rmsd_equals_the_kabsch_rmsd_of_every_pair(frames):
    result = trajectory.pairwise_rmsd(frames)
    assert result.estimand is None or isinstance(result.estimand, str)
    matrix = np.asarray(result.value)
    assert matrix.shape == (12, 12)
    expected = np.array([[kabsch_rmsd(a, b) for b in frames] for a in frames])
    assert np.allclose(matrix, expected, atol=2e-5)
    assert np.allclose(matrix, matrix.T)
    assert np.all(np.diag(matrix) == 0.0)


def test_a_rotated_and_shifted_copy_is_at_distance_zero(frames):
    angle = 0.7
    rotation = np.array(
        [[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, 1]]
    )
    moved = (frames[0].astype(np.float64) @ rotation.T + [5.0, -2.0, 1.0]).astype(np.float32)
    ensemble = np.stack([frames[0], moved])
    assert float(np.asarray(trajectory.pairwise_rmsd(ensemble).value)[0, 1]) == pytest.approx(
        0.0, abs=1e-5
    )


def test_pca_eigenvalues_are_those_of_the_covariance_of_the_coordinates(frames):
    result = trajectory.pca(frames, components=4)
    flat = frames.reshape(len(frames), -1).astype(np.float64)
    centred = flat - flat.mean(axis=0)
    # Normalised by the number of frames less one.
    expected = np.sort(np.linalg.eigvalsh(centred.T @ centred / (len(frames) - 1)))[::-1][:4]
    pca = result.value
    assert np.allclose(pca.eigenvalues, expected, rtol=1e-6)
    assert np.allclose(pca.mean, flat.mean(axis=0), atol=1e-6)
    assert np.asarray(pca.components).shape == (4, 18)
    assert np.asarray(pca.projections).shape == (12, 4)
    # The axes are orthonormal and the projections are the centred data on them.
    axes = np.asarray(pca.components)
    assert np.allclose(axes @ axes.T, np.eye(4), atol=1e-6)
    assert np.allclose(np.asarray(pca.projections), centred @ axes.T, atol=1e-5)
    # The variance of each projection is its eigenvalue.
    assert np.allclose(
        np.var(np.asarray(pca.projections), axis=0, ddof=1), pca.eigenvalues, rtol=1e-5
    )


def test_pca_fits_frames_first_when_asked_and_refuses_an_unknown_fit(frames):
    angle = 1.1
    rotation = np.array(
        [[1, 0, 0], [0, math.cos(angle), -math.sin(angle)], [0, math.sin(angle), math.cos(angle)]]
    )
    rotated = np.stack([frames[0], (frames[0] @ rotation.T).astype(np.float32)] * 3)
    fitted = trajectory.pca(rotated, components=1, fit="mean").value
    assert fitted.eigenvalues[0] == pytest.approx(0.0, abs=1e-6)
    unfitted = trajectory.pca(rotated, components=1).value
    assert unfitted.eigenvalues[0] > 0.1
    by_reference = trajectory.pca(rotated, components=1, reference=frames[0]).value
    assert by_reference.eigenvalues[0] == pytest.approx(0.0, abs=1e-6)
    with pytest.raises(molframe.MolframeValueError, match="not a fit"):
        trajectory.pca(rotated, components=1, fit="best")  # type: ignore[arg-type]


def test_the_mean_structure_of_copies_of_one_shape_is_that_shape(frames):
    ensemble = np.stack([frames[0]] * 4)
    mean = np.asarray(trajectory.mean_structure(ensemble).value, dtype=np.float64)
    assert kabsch_rmsd(mean, frames[0]) == pytest.approx(0.0, abs=1e-5)


def test_kmeans_recovers_separated_groups_and_their_means():
    rng = np.random.default_rng(3)
    first = rng.normal(0.0, 0.1, size=(15, 2))
    second = rng.normal(10.0, 0.1, size=(15, 2))
    observations = np.ascontiguousarray(np.vstack([first, second]))
    result = trajectory.kmeans(observations, initial_centres=[0, 15]).value
    labels = np.asarray(result.labels)
    assert set(labels[:15]) == {0}
    assert set(labels[15:]) == {1}
    assert np.allclose(np.asarray(result.centres)[0], first.mean(axis=0))
    assert np.allclose(np.asarray(result.centres)[1], second.mean(axis=0))
    expected = ((first - first.mean(axis=0)) ** 2).sum() + (
        (second - second.mean(axis=0)) ** 2
    ).sum()
    assert result.inertia == pytest.approx(expected)
    with pytest.raises(molframe.MolframeError) as invalid:
        trajectory.kmeans(observations, initial_centres=[0, 0])
    assert invalid.value.code == "MOLFRAME-E5101"


def line_distances(points: list[float]) -> np.ndarray:
    values = np.array(points)
    return np.abs(values[:, None] - values[None, :])


def test_agglomerative_clustering_groups_a_line_by_its_gaps():
    distances = line_distances([0.0, 0.1, 0.2, 5.0, 5.1, 9.0])
    for linkage in ("single", "complete", "average"):
        result = trajectory.agglomerative_clustering(distances, clusters=3, linkage=linkage).value
        assert sorted(sorted(group) for group in result.members) == [[0, 1, 2], [3, 4], [5]]
        assert sorted(result.medoids) == [1, 3, 5] or sorted(result.medoids) == [1, 4, 5]
    with pytest.raises(molframe.MolframeError, match="not a linkage"):
        trajectory.agglomerative_clustering(distances, clusters=3, linkage="ward")


def test_dbscan_marks_the_isolated_point_as_noise():
    distances = line_distances([0.0, 0.1, 0.2, 5.0, 5.1, 5.2, 20.0])
    result = trajectory.dbscan(distances, epsilon=0.5, minimum_points=2).value
    labels = np.asarray(result.labels)
    assert labels[6] == -1
    assert labels[0] == labels[1] == labels[2] != labels[3]
    assert labels[3] == labels[4] == labels[5]
    assert len(result) == 2


def test_the_medoid_is_the_member_nearest_the_others():
    distances = line_distances([0.0, 1.0, 2.0, 10.0])
    assert trajectory.medoid(distances, [0, 1, 2]).value == 1


def test_diffusion_map_of_two_disconnected_clusters_separates_them_in_its_first_coordinate():
    distances = line_distances([0.0, 0.1, 0.2, 8.0, 8.1, 8.2])
    result = trajectory.diffusion_map(distances, epsilon=1.0, dimensions=1).value
    coordinates = np.asarray(result.coordinates)[:, 0]
    assert coordinates.shape == (6,)
    # The affinity between the clusters is exp(-64): numerically positive, so the graph is
    # one component, but the first non-constant eigenvalue is one to within that and its
    # eigenvector is the indicator of a cluster.
    assert result.eigenvalues[0] == pytest.approx(1.0, abs=1e-9)
    assert abs(coordinates[:3].mean() - coordinates[3:].mean()) > 0.1
    assert np.ptp(coordinates[:3]) < 1e-9
    assert np.ptp(coordinates[3:]) < 1e-9
    assert set(np.asarray(result.graph_components).tolist()) == {0}


def test_msd_of_uniform_motion_is_the_squared_distance_travelled():
    velocity = np.array([0.5, -0.25, 1.0])
    steps = np.arange(10)[:, None, None] * velocity[None, None, :]
    positions = np.ascontiguousarray(np.repeat(steps, 3, axis=1), dtype=np.float32)
    result = trajectory.msd(positions, maximum_lag=4).value
    lag = np.asarray(result.lag)
    assert lag.tolist() == [0, 1, 2, 3, 4]
    assert np.allclose(np.asarray(result.value), (lag**2) * float(velocity @ velocity), atol=1e-5)
    assert result.observations[0] == 30
    only = trajectory.msd(positions, maximum_lag=2, atoms=[1]).value
    assert only.observations[0] == 10
    with pytest.raises(molframe.MolframeError):
        trajectory.msd(positions, maximum_lag=10)


def test_group_variance_is_the_mean_of_each_atoms_variance_about_its_own_mean():
    rng = np.random.default_rng(5)
    positions = rng.normal(size=(20, 4, 3)).astype(np.float32)
    groups = [[0, 1], [2, 3]]
    found = trajectory.group_variance(positions, groups).value
    for group, result in zip(groups, found, strict=True):
        per_atom = positions[:, group, :].astype(np.float64).var(axis=0)
        axis = per_atom.mean(axis=0)
        assert np.allclose(result.variance_by_axis, axis, rtol=1e-5)
        assert result.rms_fluctuation == pytest.approx(math.sqrt(axis.sum()), rel=1e-5)
        assert result.atoms == group


def test_block_convergence_reports_block_and_cumulative_means():
    values = [1.0, 3.0, 5.0, 7.0, 9.0]
    result = trajectory.block_convergence(values, block_size=2).value
    assert np.asarray(result.start).tolist() == [0, 2, 4]
    assert np.asarray(result.end).tolist() == [2, 4, 5]
    assert np.allclose(np.asarray(result.block_mean), [2.0, 6.0, 9.0])
    assert np.allclose(np.asarray(result.cumulative_mean), [2.0, 4.0, 5.0])
    with pytest.raises(molframe.MolframeError):
        trajectory.block_convergence(values, block_size=2, remainder="reject")
    with pytest.raises(molframe.MolframeError, match="remainder"):
        trajectory.block_convergence(values, block_size=2, remainder="drop")


def test_harmonic_similarity_of_a_shifted_copy_is_the_gaussian_mean_term():
    rng = np.random.default_rng(11)
    first = rng.normal(size=(50, 1))
    second = first + 2.0
    regularization = 1e-3
    result = trajectory.harmonic_similarity(first, second, regularization=regularization).value
    variance = first.var() + regularization  # covariance normalised by the observation count
    assert result.covariance_term == pytest.approx(0.0, abs=1e-12)
    assert result.mean_term == pytest.approx(2.0**2 / (8.0 * variance), rel=1e-9)
    assert result.distance == pytest.approx(result.mean_term)
    assert result.similarity == pytest.approx(math.exp(-result.distance))
    same = trajectory.harmonic_similarity(first, first, regularization=regularization).value
    assert same.similarity == pytest.approx(1.0)


def test_population_similarity_is_one_for_equal_and_zero_for_disjoint_populations():
    assert trajectory.population_similarity([0, 1, 0, 1], [1, 0, 1, 0], clusters=2).value == (
        pytest.approx(1.0)
    )
    assert trajectory.population_similarity([0, 0], [1, 1], clusters=2).value == pytest.approx(0.0)
    with pytest.raises(molframe.MolframeError):
        trajectory.population_similarity([0, 5], [1, 1], clusters=2)


def test_path_similarity_of_a_path_with_itself_is_zero_and_of_a_shifted_one_is_the_shift():
    path = np.random.default_rng(2).normal(size=(5, 4, 3)).astype(np.float32)
    same = trajectory.path_similarity(path, path)
    assert same.hausdorff_distance == pytest.approx(0.0, abs=1e-5)
    assert same.discrete_frechet_distance == pytest.approx(0.0, abs=1e-5)
    shifted = path + np.float32([3.0, 0.0, 0.0])
    cartesian = trajectory.path_similarity(path, shifted, metric="cartesian_rmsd")
    assert cartesian.hausdorff_distance <= 3.0 + 1e-5
    assert cartesian.discrete_frechet_distance == pytest.approx(3.0, abs=1e-4)
    fitted = trajectory.path_similarity(path, shifted)
    assert fitted.hausdorff_distance == pytest.approx(0.0, abs=1e-4)
    with pytest.raises(molframe.MolframeError, match="path frame metric"):
        trajectory.path_similarity(path, path, metric="rmsd")
