"""Element data, radii and validation through the Python contract."""

import math
from pathlib import Path

import numpy
import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-py" / "tests" / "data"


def test_element_properties_are_case_insensitive_and_complete():
    carbon = molframe.chemistry.element("c")
    assert (carbon.symbol, carbon.atomic_number, carbon.period) == ("C", 6, 2)
    assert carbon.atomic_weight == pytest.approx(12.011)
    assert molframe.chemistry.element("Fe").group == 8
    with pytest.raises(ValueError):
        molframe.chemistry.element("Zz")


def test_radius_sets_differ_and_unknown_sets_are_rejected():
    bondi = molframe.chemistry.vdw_radius("C", radii="bondi")
    assert bondi == pytest.approx(1.7)
    assert molframe.chemistry.vdw_radius("C", radii="alvarez") != bondi
    with pytest.raises(ValueError):
        molframe.chemistry.vdw_radius("C", radii="made_up")


def test_per_atom_radii_feed_the_surface_kernels():
    structure = molframe.read(DATA / "basic.cif")
    radii = molframe.chemistry.vdw_radii(structure)
    assert radii.shape == (structure.atom_count,)
    areas = molframe.surface.sasa(structure.coordinates, radii)
    assert areas.shape == radii.shape
    assert numpy.nansum(areas) > 0


def test_clashes_report_overlapping_pairs_only():
    structure = molframe.read(DATA / "basic.cif")
    loose = molframe.validation.clashes(structure, tolerance=10.0)
    tight = molframe.validation.clashes(structure, tolerance=0.0)
    assert len(tight) >= len(loose)
    assert len(loose) == len(loose.first) == len(loose.second) == len(loose.overlap)
    assert all(a < b for a, b in zip(tight.first, tight.second, strict=True))
    assert all(math.isfinite(value) and value > 0 for value in tight.overlap)
    with pytest.raises(ValueError):
        molframe.validation.clashes(structure, tolerance=-1.0)
    with pytest.raises(ValueError):
        molframe.validation.clashes(structure, backend="octree")


def test_neighbor_pairs_agree_across_backends_and_match_brute_force():
    rng = numpy.random.default_rng(7)
    xyz = rng.uniform(0, 12, size=(200, 3)).astype(numpy.float32)
    reference = None
    for backend in ("auto", "cell", "kd_tree", "brute_force"):
        first, second, distance = molframe.spatial.neighbor_pairs(xyz, 3.0, backend=backend)
        pairs = list(zip(first.tolist(), second.tolist(), strict=True))
        assert pairs == sorted(pairs) and all(a < b for a, b in pairs)
        if reference is None:
            reference = pairs
        assert pairs == reference
    delta = xyz[:, None, :] - xyz[None, :, :]
    full = numpy.sqrt((delta**2).sum(-1))
    expected = sorted(
        (int(i), int(j)) for i in range(200) for j in range(i + 1, 200) if full[i, j] <= 3.0
    )
    assert reference == expected
    assert numpy.all(distance <= 3.0 + 1e-4)


def test_neighbor_pairs_reject_a_bad_cutoff_and_backend():
    xyz = numpy.zeros((2, 3), dtype=numpy.float32)
    with pytest.raises(ValueError):
        molframe.spatial.neighbor_pairs(xyz, -1.0)
    with pytest.raises(ValueError):
        molframe.spatial.neighbor_pairs(xyz, 1.0, backend="octree")
