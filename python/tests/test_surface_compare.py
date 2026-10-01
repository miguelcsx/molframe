"""Surface and comparison kernels through the Python contract."""

import math

import numpy
import pytest

import molframe


def _helix(count=12):
    t = numpy.arange(count, dtype=numpy.float32)
    return numpy.stack([2.3 * numpy.cos(t), 2.3 * numpy.sin(t), 1.5 * t], axis=1)


def test_a_lone_sphere_has_the_analytic_accessible_area():
    xyz = numpy.zeros((1, 3), dtype=numpy.float32)
    radii = numpy.array([1.7], dtype=numpy.float32)
    expected = 4.0 * math.pi * (1.7 + 1.4) ** 2
    assert molframe.surface.sasa(xyz, radii)[0] == pytest.approx(expected, rel=1e-3)
    assert molframe.surface.lee_richards(xyz, radii)[0] == pytest.approx(
        expected, rel=1e-2
    )


def test_buried_atoms_lose_area_and_bad_radii_are_rejected():
    xyz = numpy.array([[0, 0, 0], [1.5, 0, 0]], dtype=numpy.float32)
    radii = numpy.array([1.7, 1.7], dtype=numpy.float32)
    pair = molframe.surface.sasa(xyz, radii)
    lone = molframe.surface.sasa(xyz[:1], radii[:1])[0]
    assert pair.sum() < 2 * lone
    with pytest.raises(ValueError):
        molframe.surface.sasa(xyz, radii[:1])


def test_a_hollow_shell_of_atoms_encloses_a_cavity():
    golden = numpy.pi * (3.0 - math.sqrt(5.0))
    count = 120
    k = numpy.arange(count, dtype=numpy.float64)
    y = 1 - 2 * (k + 0.5) / count
    r = numpy.sqrt(1 - y * y)
    shell = 5.0 * numpy.stack([r * numpy.cos(golden * k), y, r * numpy.sin(golden * k)], 1)
    xyz = shell.astype(numpy.float32)
    radii = numpy.full(count, 2.4, dtype=numpy.float32)
    found = molframe.surface.cavities(xyz, radii, probe=1.0, resolution=0.5)
    assert found and found[0][0] > 1.0
    assert molframe.surface.cavities(xyz[:1], radii[:1]) == []


def test_superposition_scores_are_rigid_motion_invariant():
    reference = _helix()
    angle = 0.7
    rotation = numpy.array(
        [[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, 1]],
        dtype=numpy.float32,
    )
    model = reference @ rotation.T + numpy.array([5, -3, 2], dtype=numpy.float32)
    assert molframe.compare.tm_score(model, reference) == pytest.approx(1.0, abs=1e-4)
    assert molframe.compare.gdt_ts(model, reference) == pytest.approx(1.0)
    assert molframe.compare.gdt_ha(model, reference) == pytest.approx(1.0)


def test_weighted_rmsd_follows_the_weights_and_rejects_mismatches():
    reference = numpy.zeros((2, 3), dtype=numpy.float32)
    model = numpy.array([[1, 0, 0], [3, 0, 0]], dtype=numpy.float32)
    even = molframe.compare.weighted_rmsd(model, reference, numpy.array([1.0, 1.0]))
    first = molframe.compare.weighted_rmsd(model, reference, numpy.array([1.0, 0.0]))
    assert even == pytest.approx(math.sqrt(5.0))
    assert first == pytest.approx(1.0)
    with pytest.raises(ValueError):
        molframe.compare.weighted_rmsd(model, reference, numpy.array([1.0]))
    with pytest.raises(ValueError):
        molframe.compare.tm_score(model[:1], reference[:1])
