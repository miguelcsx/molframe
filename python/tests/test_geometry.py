"""Geometry kernels, each checked against an independent NumPy formulation."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import geometry

BENCH = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"
rng = np.random.default_rng(20261003)


def cloud(n=40, scale=8.0):
    return (rng.normal(size=(n, 3)) * scale).astype(np.float32)


def test_distances_equal_the_row_norms():
    a, b = cloud(), cloud()
    expected = np.linalg.norm(a.astype(np.float64) - b.astype(np.float64), axis=1)
    assert np.allclose(geometry.distances(a, b), expected, atol=1e-6)
    with pytest.raises(molframe.MolframeError) as raised:
        geometry.distances(a, b[:-1])
    assert raised.value.code == "MOLFRAME-E5102"


def test_angles_equal_the_arccosine_of_the_normalised_dot_product():
    a, v, c = cloud(), cloud(), cloud()
    u, w = (a - v).astype(np.float64), (c - v).astype(np.float64)
    cosine = (u * w).sum(1) / (np.linalg.norm(u, axis=1) * np.linalg.norm(w, axis=1))
    assert np.allclose(geometry.angles(a, v, c), np.degrees(np.arccos(cosine)), atol=1e-4)
    coincident = np.zeros((1, 3), dtype=np.float32)
    assert np.isnan(geometry.angles(coincident, coincident, c[:1])[0])


def reference_dihedral(p0, p1, p2, p3):
    b0, b1, b2 = -(p1 - p0), p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - (b0 * b1).sum(1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(1, keepdims=True) * b1
    x = (v * w).sum(1)
    y = (np.cross(b1, v) * w).sum(1)
    return np.degrees(np.arctan2(y, x))


def test_dihedrals_equal_the_atan2_formulation_and_are_signed():
    quad = [cloud() for _ in range(4)]
    expected = reference_dihedral(*[q.astype(np.float64) for q in quad])
    found = geometry.dihedrals(*quad)
    assert np.allclose(found, expected, atol=1e-3)
    assert (found > 0).any()
    assert (found < 0).any()
    cis = [np.array([[0, 1, 0]], np.float32), np.array([[0, 0, 0]], np.float32)]
    cis += [np.array([[1, 0, 0]], np.float32), np.array([[1, 1, 0]], np.float32)]
    assert abs(geometry.dihedrals(*cis)[0]) < 1e-4
    trans = [*cis[:3], np.array([[1, -1, 0]], np.float32)]
    assert abs(abs(geometry.dihedrals(*trans)[0]) - 180.0) < 1e-4


def test_mass_weighted_moments_equal_their_definitions():
    points, masses = cloud(), rng.uniform(1.0, 16.0, 40)
    p = points.astype(np.float64)
    centre = (p * masses[:, None]).sum(0) / masses.sum()
    assert np.allclose(geometry.centre_of_mass(points, masses), centre)
    assert np.allclose(geometry.centre_of_mass(points), p.mean(0))
    offsets = p - centre
    rg = np.sqrt((masses * (offsets**2).sum(1)).sum() / masses.sum())
    assert geometry.radius_of_gyration(points, masses) == pytest.approx(rg)
    tensor = sum(
        m * ((o @ o) * np.eye(3) - np.outer(o, o)) for m, o in zip(masses, offsets, strict=True)
    )
    found = geometry.inertia_tensor(points, masses)
    assert np.allclose(found, tensor)
    moments, axes = geometry.principal_axes(points, masses)
    assert np.allclose(np.sort(moments), np.sort(np.linalg.eigvalsh(tensor)))
    assert np.allclose(found @ axes, axes * moments)
    assert geometry.centre_of_mass(np.empty((0, 3), np.float32)) is None


def test_shape_descriptors_distinguish_a_line_from_a_sphere():
    line = np.column_stack([np.linspace(0, 10, 30), np.zeros(30), np.zeros(30)]).astype(np.float32)
    ball = rng.normal(size=(4000, 3)).astype(np.float32)
    assert geometry.asphericity(line) == pytest.approx(1.0, abs=1e-6)
    assert geometry.asphericity(ball) < 0.02
    assert geometry.shape_parameter(line) == pytest.approx(2.0 / 2.0, abs=1.0)
    values, axes = geometry.gyration_axes(line)
    assert values[0] > 1.0
    assert np.allclose(np.abs(axes[:, 0]), [1, 0, 0], atol=1e-6)


def kabsch(mobile, reference):
    m, r = mobile.astype(np.float64), reference.astype(np.float64)
    mc, rc = m.mean(0), r.mean(0)
    u, _, vt = np.linalg.svd((m - mc).T @ (r - rc))
    sign = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1, 1, sign]) @ u.T
    moved = (m - mc) @ rotation.T + rc
    return rotation, rc - rotation @ mc, np.sqrt(((moved - r) ** 2).sum(1).mean())


def random_rotation():
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_superposition_recovers_a_known_rigid_motion_and_agrees_with_kabsch():
    reference = cloud()
    rotation = random_rotation()
    shift = np.array([4.0, -2.0, 7.5])
    mobile = (reference.astype(np.float64) @ rotation.T + shift).astype(np.float32)
    fit = geometry.superpose(mobile, reference)
    assert fit.rmsd < 1e-4
    assert np.allclose(fit.apply(mobile), reference, atol=1e-3)
    noisy = mobile + rng.normal(scale=0.3, size=mobile.shape).astype(np.float32)
    fit = geometry.superpose(noisy, reference)
    expected_rotation, expected_shift, expected_rmsd = kabsch(noisy, reference)
    assert np.allclose(fit.rotation, expected_rotation, atol=1e-6)
    assert np.allclose(fit.translation, expected_shift, atol=1e-4)
    assert fit.rmsd == pytest.approx(expected_rmsd, abs=1e-6)
    assert geometry.rmsd_after_fit(noisy, reference) == pytest.approx(expected_rmsd, abs=1e-6)
    assert geometry.rmsd(noisy, reference) > expected_rmsd
    assert np.linalg.det(fit.rotation) == pytest.approx(1.0)


def test_a_degenerate_set_cannot_be_superposed():
    line = np.column_stack([np.arange(5.0), np.zeros(5), np.zeros(5)]).astype(np.float32)
    with pytest.raises(molframe.MolframeError) as raised:
        geometry.superpose(line, line)
    assert raised.value.code is not None


def test_plane_fit_finds_the_normal_and_the_deviation():
    plane = np.array([[x, y, 0.0] for x in range(4) for y in range(4)], np.float32)
    centre, normal = geometry.best_fit_plane(plane)
    assert np.allclose(centre, [1.5, 1.5, 0.0])
    assert np.allclose(np.abs(normal), [0, 0, 1], atol=1e-9)
    assert geometry.plane_deviation(plane) == pytest.approx(0.0, abs=1e-9)
    lifted = plane.copy()
    lifted[:, 2] = 0.5 * (-1.0) ** (plane[:, 0] + plane[:, 1])
    assert geometry.plane_deviation(lifted) == pytest.approx(0.5, abs=1e-6)
    assert geometry.best_fit_plane(plane[:2]) is None


def test_rmsf_equals_the_per_atom_standard_deviation_about_the_mean():
    frames = (rng.normal(size=(30, 12, 3)) * 1.5).astype(np.float32)
    expected = np.sqrt(((frames - frames.mean(0)) ** 2).sum(2).mean(0))
    assert np.allclose(geometry.rmsf(frames), expected, atol=1e-5)
    with pytest.raises(molframe.MolframeError) as raised:
        geometry.rmsf(np.empty((0, 4, 3), np.float32))
    assert raised.value.code == "MOLFRAME-E5103"


def test_backbone_torsions_equal_the_dihedrals_of_the_named_atoms(with_roles):
    ubiquitin = molframe.read(BENCH / "1ubq.cif")
    table = geometry.backbone_torsions(with_roles(ubiquitin))
    assert table.names == ["residue", "phi", "psi", "omega"]

    def point(index, name):
        found = ubiquitin.residues[index].atom(name)
        return None if found is None else np.array([found.coordinate], dtype=np.float32)

    checked = 0
    for row, residue in enumerate(table["residue"]):
        i = int(residue)
        if 0 < i < ubiquitin.residue_count - 1:
            parts = [point(i - 1, "C"), point(i, "N"), point(i, "CA"), point(i, "C")]
            if all(part is not None for part in parts):
                assert table["phi"][row] == pytest.approx(
                    float(geometry.dihedrals(*parts)[0]), abs=1e-3
                )
                psi = [point(i, "N"), point(i, "CA"), point(i, "C"), point(i + 1, "N")]
                if all(part is not None for part in psi):
                    assert table["psi"][row] == pytest.approx(
                        float(geometry.dihedrals(*psi)[0]), abs=1e-3
                    )
                checked += 1
    assert checked > 60
    assert np.isnan(table["phi"][0])
    assert np.isnan(table["psi"][-1])
