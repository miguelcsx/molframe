"""Element data, radii and validation through the Python contract."""

import math
from pathlib import Path

import numpy as np
import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-py" / "tests" / "data"
BACKEND_CHOICES = r"backend: octree .expected brute-force, cell-list, kd-tree, neighbor-list, auto"


def test_element_properties_are_case_insensitive_and_complete():
    carbon = molframe.chemistry.element("c")
    assert (carbon.symbol, carbon.atomic_number, carbon.period) == ("C", 6, 2)
    assert carbon.atomic_weight == pytest.approx(12.011)
    assert molframe.chemistry.element("Fe").group == 8
    with pytest.raises(ValueError, match='unknown element "Zz"'):
        molframe.chemistry.element("Zz")


def test_radius_sets_differ_and_unknown_sets_are_rejected():
    bondi = molframe.chemistry.vdw_radius("C", radii="bondi")
    assert bondi == pytest.approx(1.7)
    assert molframe.chemistry.vdw_radius("C", radii="alvarez") != bondi
    with pytest.raises(molframe.MolframeError, match="E7102"):
        molframe.chemistry.vdw_radius("C", radii="made_up")


def test_covalent_radii_come_from_their_own_table_and_reject_unknown_names():
    # The test profile names the Pyykkö single-bond table for covalent radii.
    assert molframe.chemistry.covalent_radius("C") == pytest.approx(0.75)
    assert molframe.chemistry.covalent_radius("c") == pytest.approx(0.75)
    assert molframe.chemistry.covalent_radius("C") != molframe.chemistry.vdw_radius("C")
    with pytest.raises(ValueError, match='unknown element "Zz"'):
        molframe.chemistry.covalent_radius("Zz")
    with pytest.raises(molframe.MolframeError, match="E7102"):
        molframe.chemistry.covalent_radius("C", radii="made_up")


def test_per_atom_covalent_radii_are_nan_where_a_named_table_is_silent(tmp_path):
    table = tmp_path / "carbon-only.covalent"
    table.write_text("name: carbon-only\nversion: 1\nsource: test\nlicence: CC0-1.0\n\nC 0.8\n")
    coordinates = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], dtype=np.float32)
    structure = molframe.chemistry.Molecule(["C", "O"], coordinates).to_structure()
    radii = molframe.chemistry.covalent_radii(structure, radii=str(table))
    assert radii.dtype == np.float32
    assert radii[0] == pytest.approx(0.8)
    assert math.isnan(radii[1])
    assert molframe.chemistry.covalent_radius("O", radii=str(table)) is None


def test_per_atom_radii_feed_the_surface_kernels():
    structure = molframe.read(DATA / "basic.cif")
    radii = molframe.chemistry.vdw_radii(structure)
    assert radii.shape == (structure.atom_count,)
    areas = molframe.surface.sasa(structure.coordinates, radii)
    assert areas.shape == radii.shape
    assert np.nansum(areas) > 0


def test_clashes_report_overlapping_pairs_only():
    structure = molframe.read(DATA / "basic.cif")
    loose = molframe.validation.clashes(structure, tolerance=10.0).value
    tight = molframe.validation.clashes(structure, tolerance=0.0).value
    assert len(tight) >= len(loose)
    assert len(loose) == len(loose.first) == len(loose.second) == len(loose.overlap)
    assert all(a < b for a, b in zip(tight.first, tight.second, strict=True))
    assert all(math.isfinite(value) and value > 0 for value in tight.overlap)
    with pytest.raises(ValueError, match="tolerance must be finite and non-negative"):
        molframe.validation.clashes(structure, tolerance=-1.0)
    with pytest.raises(ValueError, match=BACKEND_CHOICES):
        molframe.validation.clashes(structure, backend="octree")


def test_neighbor_pairs_agree_across_backends_and_match_brute_force():
    rng = np.random.default_rng(7)
    xyz = rng.uniform(0, 12, size=(200, 3)).astype(np.float32)
    reference = None
    for backend in ("auto", "cell", "kd_tree", "brute_force"):
        first, second, distance = molframe.spatial.neighbor_pairs(xyz, 3.0, backend=backend)
        pairs = list(zip(first.tolist(), second.tolist(), strict=True))
        assert pairs == sorted(pairs)
        assert all(a < b for a, b in pairs)
        if reference is None:
            reference = pairs
        assert pairs == reference
    delta = xyz[:, None, :] - xyz[None, :, :]
    full = np.sqrt((delta**2).sum(-1))
    expected = sorted(
        (int(i), int(j)) for i in range(200) for j in range(i + 1, 200) if full[i, j] <= 3.0
    )
    assert reference == expected
    assert np.all(distance <= 3.0 + 1e-4)


def test_neighbor_pairs_reject_a_bad_cutoff_and_backend():
    xyz = np.zeros((2, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="cutoff must be finite and non-negative"):
        molframe.spatial.neighbor_pairs(xyz, -1.0)
    with pytest.raises(ValueError, match=BACKEND_CHOICES):
        molframe.spatial.neighbor_pairs(xyz, 1.0, backend="octree")
