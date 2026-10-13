"""Periodic policy and unit-cell plumbing, checked against a NumPy image search."""

from itertools import product
from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "crates" / "molframe-bench" / "data"


def basis_of(cell):
    a, b, c = cell.lengths
    alpha, beta, gamma = np.radians(cell.angles)
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c * c - cx * cx - cy * cy)
    return np.array(
        [[a, 0, 0], [b * np.cos(gamma), b * np.sin(gamma), 0], [cx, cy, cz]], dtype=np.float64
    )


def image_distances(xyz, cell):
    """Minimum distance over the 27 nearest lattice images for every atom pair."""
    basis = basis_of(cell)
    shifts = np.array(list(product((-1, 0, 1), repeat=3)), dtype=np.float64) @ basis
    delta = xyz[:, None, :] - xyz[None, :, :]
    fractional = delta @ np.linalg.inv(basis)
    delta = (fractional - np.round(fractional)) @ basis
    best = np.full(delta.shape[:2], np.inf)
    for shift in shifts:
        best = np.minimum(best, np.linalg.norm(delta + shift, axis=2))
    return best


@pytest.fixture(scope="module")
def crambin():
    return molframe.read(BENCH / "1crn.cif")


def pair_count(distances, cutoff):
    return int(((distances <= cutoff) & np.triu(np.ones_like(distances, dtype=bool), 1)).sum())


def test_minimum_image_contacts_match_an_image_search_and_differ_from_open_space(crambin):
    xyz = np.asarray(crambin.coordinates, dtype=np.float64)
    cell = crambin.cell
    cutoff = 5.0
    periodic = pair_count(image_distances(xyz, cell), cutoff)
    flat = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
    assert periodic > pair_count(flat, cutoff)
    for mode in ("minimum_image", "pbc"):
        policy = molframe.AnalysisPolicy(periodic=mode)
        found = analysis.contacts(crambin, cutoff, policy=policy)
        assert len(found.value) == periodic
    plain = analysis.contacts(crambin, cutoff)
    assert len(plain.value) == pair_count(flat, cutoff)


def test_pbc_refuses_a_cutoff_larger_than_half_the_cell(crambin):
    policy = molframe.AnalysisPolicy(periodic="pbc")
    with pytest.raises(Exception, match=r"image|cutoff|cell"):
        analysis.contacts(crambin, 12.0, policy=policy)


def test_an_analysis_that_cannot_apply_the_cell_refuses_the_periodic_policy(crambin):
    policy = molframe.AnalysisPolicy(periodic="minimum_image")
    with pytest.raises(Exception, match="periodic"):
        analysis.half_sphere_exposure(crambin, radius=13.0, policy=policy)


def test_spatial_pairs_through_a_cell_match_the_image_search(crambin):
    xyz = np.asarray(crambin.coordinates, dtype=np.float32)
    cutoff = 5.0
    first, second, distance = molframe.spatial.neighbor_pairs(xyz, cutoff, cell=crambin.cell)
    reference = image_distances(xyz.astype(np.float64), crambin.cell)
    assert len(first) == pair_count(reference, cutoff)
    np.testing.assert_allclose(distance, reference[first, second], atol=2e-4)
    left, right, cross = molframe.spatial.cross_pairs(xyz[:20], xyz[20:], cutoff, cell=crambin.cell)
    expected = reference[:20, 20:]
    assert len(left) == int((expected <= cutoff).sum())
    np.testing.assert_allclose(cross, expected[left, right], atol=2e-4)


def test_orthorhombic_cell_wraps_two_atoms_across_the_boundary():
    cell = molframe.crystal.UnitCell([10.0, 10.0, 10.0], [90.0, 90.0, 90.0])
    xyz = np.array([[0.5, 5, 5], [9.5, 5, 5]], dtype=np.float32)
    _, _, distance = molframe.spatial.neighbor_pairs(xyz, 2.0, cell=cell)
    np.testing.assert_allclose(distance, [1.0], atol=1e-5)
    assert len(molframe.spatial.neighbor_pairs(xyz, 2.0)[0]) == 0
