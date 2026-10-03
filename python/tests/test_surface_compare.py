"""Surface and comparison kernels through the Python contract."""

import math

import numpy as np
import pytest

import molframe


def _helix(count=12):
    t = np.arange(count, dtype=np.float32)
    return np.stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * t], axis=1)


def test_a_lone_sphere_has_the_analytic_accessible_area():
    xyz = np.zeros((1, 3), dtype=np.float32)
    radii = np.array([1.7], dtype=np.float32)
    expected = 4.0 * math.pi * (1.7 + 1.4) ** 2
    assert molframe.surface.sasa(xyz, radii)[0] == pytest.approx(expected, rel=1e-3)
    assert molframe.surface.lee_richards(xyz, radii)[0] == pytest.approx(expected, rel=1e-2)


def test_buried_atoms_lose_area_and_bad_radii_are_rejected():
    xyz = np.array([[0, 0, 0], [1.5, 0, 0]], dtype=np.float32)
    radii = np.array([1.7, 1.7], dtype=np.float32)
    pair = molframe.surface.sasa(xyz, radii)
    lone = molframe.surface.sasa(xyz[:1], radii[:1])[0]
    assert pair.sum() < 2 * lone
    with pytest.raises(ValueError, match="expected 2 radii to match the positions, found 1"):
        molframe.surface.sasa(xyz, radii[:1])


def test_a_hollow_shell_of_atoms_encloses_a_cavity():
    golden = np.pi * (3.0 - math.sqrt(5.0))
    count = 120
    k = np.arange(count, dtype=np.float64)
    y = 1 - 2 * (k + 0.5) / count
    r = np.sqrt(1 - y * y)
    shell = 5.0 * np.stack([r * np.cos(golden * k), y, r * np.sin(golden * k)], 1)
    xyz = shell.astype(np.float32)
    radii = np.full(count, 2.4, dtype=np.float32)
    found = molframe.surface.cavities(xyz, radii, probe=1.0, resolution=0.5)
    assert found
    assert found[0][0] > 1.0
    assert molframe.surface.cavities(xyz[:1], radii[:1]) == []


def test_superposition_scores_are_rigid_motion_invariant():
    reference = _helix()
    angle = 0.7
    rotation = np.array(
        [[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, 1]],
        dtype=np.float32,
    )
    model = reference @ rotation.T + np.array([5, -3, 2], dtype=np.float32)
    assert molframe.compare.tm_score(model, reference) == pytest.approx(1.0, abs=1e-4)
    assert molframe.compare.gdt_ts(model, reference) == pytest.approx(1.0)
    assert molframe.compare.gdt_ha(model, reference) == pytest.approx(1.0)


def test_weighted_rmsd_follows_the_weights_and_rejects_mismatches():
    reference = np.zeros((2, 3), dtype=np.float32)
    model = np.array([[1, 0, 0], [3, 0, 0]], dtype=np.float32)
    even = molframe.compare.weighted_rmsd(model, reference, np.array([1.0, 1.0]))
    first = molframe.compare.weighted_rmsd(model, reference, np.array([1.0, 0.0]))
    assert even == pytest.approx(math.sqrt(5.0))
    assert first == pytest.approx(1.0)
    with pytest.raises(ValueError, match="weights has 1 entries but the coordinates have 2 points"):
        molframe.compare.weighted_rmsd(model, reference, np.array([1.0]))
    with pytest.raises(ValueError, match="superposition failed: TooFewPoints"):
        molframe.compare.tm_score(model[:1], reference[:1])
