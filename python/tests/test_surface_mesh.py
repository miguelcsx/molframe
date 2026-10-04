"""Surfaces and meshes: each quantity is checked against a closed form."""

import math

import numpy as np
import pytest

import molframe
from molframe import surface


def single(radius: float):
    return np.zeros((1, 3), dtype=np.float32), np.array([radius], dtype=np.float32)


def cube():
    vertices = np.array(
        [[x, y, z] for x in (0.0, 1.0) for y in (0.0, 1.0) for z in (0.0, 1.0)], dtype=np.float32
    )
    quads = [
        (0, 1, 3, 2),
        (4, 6, 7, 5),
        (0, 4, 5, 1),
        (2, 3, 7, 6),
        (0, 2, 6, 4),
        (1, 5, 7, 3),
    ]
    faces = [(a, b, c) for a, b, c, d in quads] + [(a, c, d) for a, b, c, d in quads]
    return vertices, np.array(faces, dtype=np.uint32)


def tetrahedron(offset: float = 0.0):
    vertices = np.array(
        [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=np.float32
    ) + np.float32(offset)
    faces = np.array([[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]], dtype=np.uint32)
    return vertices, faces


def test_the_accessible_surface_of_one_atom_is_a_sphere_of_radius_plus_probe():
    positions, radii = single(2.0)
    points = surface.surface_points(positions, radii, probe=1.4, samples=480)
    assert len(points) == 480
    centre = np.asarray(points.positions, dtype=np.float64)
    assert np.allclose(np.linalg.norm(centre, axis=1), 3.4, atol=1e-4)
    normals = np.asarray(points.normals, dtype=np.float64)
    assert np.allclose(normals, centre / 3.4, atol=1e-4)
    assert set(np.asarray(points.atoms).tolist()) == {0}
    dense = surface.surface_points_at_density(positions, radii, probe=1.4, density=2.0)
    expected = 4.0 * math.pi * 3.4**2 * 2.0
    assert len(dense) == pytest.approx(expected, rel=0.02)


def test_the_solvent_excluded_surface_of_one_atom_is_its_own_sphere():
    positions, radii = single(2.0)
    mesh = surface.solvent_excluded_surface(positions, radii, probe=1.4, resolution=0.2)
    assert mesh.is_manifold
    assert mesh.area == pytest.approx(4.0 * math.pi * 2.0**2, rel=0.06)
    assert mesh.component_count == 1


def test_two_touching_atoms_bury_the_two_caps_each_takes_from_the_other():
    radius, probe, separation = 1.5, 1.4, 2.0
    positions = np.array([[0.0, 0.0, 0.0], [separation, 0.0, 0.0]], dtype=np.float32)
    radii = np.array([radius, radius], dtype=np.float32)
    result = surface.buried_surface(
        positions, radii, np.array([True, False]), probe=probe, points=2000
    )
    rho = radius + probe
    # Each accessible sphere loses a cap of height rho - separation / 2 to the other.
    expected = 2.0 * (2.0 * math.pi * rho * (rho - separation / 2.0))
    assert result.buried == pytest.approx(expected, rel=0.03)
    assert result.buried == pytest.approx(
        result.first_alone + result.second_alone - result.together
    )
    by_role = surface.buried_solvent_excluded_surface(
        positions, radii, ["first", "second"], probe=probe, resolution=0.25
    )
    assert by_role.buried > 0.0
    assert by_role.together < by_role.first_alone + by_role.second_alone
    only_first = surface.buried_solvent_excluded_surface(
        positions, radii, ["first", "excluded"], probe=probe, resolution=0.25
    )
    assert only_first.buried == pytest.approx(0.0)


def test_an_unknown_role_is_refused_by_name():
    positions, radii = single(1.5)
    with pytest.raises(molframe.MolframeError, match="not a molecule role"):
        surface.buried_solvent_excluded_surface(positions, radii, ["both"])


def test_an_atoms_depth_is_its_distance_to_the_nearest_surface_point():
    atoms = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
    points = np.array([[0.0, 0.0, 5.0], [10.0, 0.0, 0.0]], dtype=np.float32)
    assert float(surface.atom_depths(atoms, points, cell_size=3.0)[0]) == pytest.approx(
        5.0, abs=1e-4
    )
    nothing = np.zeros((0, 3), dtype=np.float32)
    assert math.isinf(float(surface.atom_depths(atoms, nothing, cell_size=3.0)[0]))
    with pytest.raises(molframe.MolframeError):
        surface.atom_depths(atoms, points, cell_size=0.0)


def test_a_closed_mesh_reports_itself_and_its_area():
    mesh = surface.Mesh(*cube())
    assert mesh.is_manifold
    assert (mesh.boundary_edges, mesh.non_manifold_edges, mesh.degenerate_faces) == (0, 0, 0)
    assert mesh.area == pytest.approx(6.0)
    assert np.asarray(mesh.faces).shape == (12, 3)
    assert np.allclose(np.linalg.norm(mesh.vertex_normals, axis=1), 1.0, atol=1e-5)
    assert not mesh.vertices.flags.writeable


def test_a_mesh_names_its_faults_rather_than_repairing_them():
    vertices, faces = cube()
    broken = surface.Mesh(vertices, np.vstack([faces, [[0, 1, 99]]]).astype(np.uint32))
    assert not broken.is_manifold
    assert broken.degenerate_faces == 1
    with pytest.raises(molframe.MolframeError):
        broken.curvatures()
    with pytest.raises(molframe.MolframeError):
        broken.geodesic_distances(0)


def test_components_are_split_by_area_and_filtered():
    first, faces = tetrahedron()
    second = first + np.float32(20.0)
    mesh = surface.Mesh(
        np.vstack([first, second * np.float32(0.5)]),
        np.vstack([faces, faces + 4]).astype(np.uint32),
    )
    parts = mesh.components()
    assert len(parts) == 2
    areas = sorted(area for _, area in parts)
    # The second tetrahedron is half the size: a quarter of the area.
    assert areas[0] == pytest.approx(areas[1] / 4.0, rel=1e-5)
    kept = mesh.filter_components(minimum_area=areas[0] * 2.0)
    assert kept.area == pytest.approx(areas[1])
    assert mesh.filter_components(maximum_components=1).area == pytest.approx(areas[1])
    with pytest.raises(molframe.MolframeError):
        mesh.filter_components(minimum_area=-1.0)


def test_distances_along_the_edges_of_a_cube():
    mesh = surface.Mesh(*cube())
    distances = np.asarray(mesh.geodesic_distances(0), dtype=np.float64)
    assert distances[0] == 0.0
    # Along edges and face diagonals a neighbour is 1 or sqrt(2) away, and the far corner no
    # more than 1 + sqrt(2).
    assert distances.max() <= 1.0 + math.sqrt(2.0) + 1e-6
    near = set(mesh.patch(0, 1.0).tolist())
    assert 0 in near
    assert all(distances[v] <= 1.0 + 1e-9 for v in near)
    assert len(near) < len(distances)
    with pytest.raises(molframe.MolframeError):
        mesh.geodesic_distances(99)
    with pytest.raises(molframe.MolframeError):
        mesh.patch(0, -1.0)


def test_every_vertex_of_a_regular_tetrahedron_curves_alike():
    curvature = surface.Mesh(*tetrahedron()).curvatures()
    assert len(curvature) == 4
    assert set(curvature.quality) == {"interior"}
    gaussian = np.asarray(curvature.gaussian)
    assert np.all(gaussian > 0.0)
    assert np.allclose(gaussian, gaussian[0], rtol=1e-6)
    assert np.allclose(np.asarray(curvature.mean), np.asarray(curvature.mean)[0], rtol=1e-6)


def test_a_mesh_is_written_as_obj_once(tmp_path):
    mesh = surface.Mesh(*cube())
    destination = tmp_path / "cube.obj"
    mesh.write_obj(destination)
    lines = destination.read_text().splitlines()
    assert sum(line.startswith("v ") for line in lines) == 8
    assert sum(line.startswith("f ") for line in lines) == 12
    with pytest.raises(molframe.MolframeIOError):
        mesh.write_obj(destination)
