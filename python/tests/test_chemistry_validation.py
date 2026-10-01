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
