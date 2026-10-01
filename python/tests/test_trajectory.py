"""Whole-trajectory reading through the Python contract."""

import numpy
import pytest

import molframe

XYZ = """2
frame 0
C 0.0 0.0 0.0
O 1.0 0.0 0.0
2
frame 1
C 0.0 0.5 0.0
O 1.0 0.5 0.0
2
frame 2
C 0.0 1.0 0.0
O 1.0 1.0 0.0
"""


def test_a_multi_frame_file_becomes_a_read_only_frame_atom_xyz_array(tmp_path):
    path = tmp_path / "walk.xyz"
    path.write_text(XYZ)
    trajectory = molframe.trajectory.read(path)
    assert (len(trajectory), trajectory.n_frames, trajectory.n_atoms) == (3, 3, 2)
    assert trajectory.format == "xyz"
    assert trajectory.positions.shape == (3, 2, 3)
    assert trajectory.positions[2, 0, 1] == pytest.approx(1.0)
    assert not trajectory.positions.flags.writeable
    assert trajectory.times.shape == (3,)
    assert trajectory.positions is trajectory.positions


def test_an_unknown_extension_needs_an_explicit_format(tmp_path):
    path = tmp_path / "walk.dat"
    path.write_text(XYZ)
    with pytest.raises(ValueError):
        molframe.trajectory.read(path)
    assert molframe.trajectory.read(path, format="xyz").n_frames == 3
    with pytest.raises(ValueError):
        molframe.trajectory.read(path, format="mdcrd")
    with pytest.raises(ValueError):
        molframe.trajectory.read(tmp_path / "missing.xyz")


def test_frames_feed_the_comparison_kernels(tmp_path):
    path = tmp_path / "walk.xyz"
    path.write_text(XYZ)
    frames = molframe.trajectory.read(path).positions
    reference = numpy.ascontiguousarray(frames[0])
    assert molframe.geometry.rmsd(numpy.ascontiguousarray(frames[2]), reference) == pytest.approx(1.0)


def _tetrahedra(drifts):
    frames = []
    for shift in drifts:
        atoms = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        lines = ["4", "frame"] + [f"C {x + shift} {y} {z}" for x, y, z in atoms]
        frames.append("\n".join(lines))
    return "\n".join(frames) + "\n"


def test_rmsd_series_measures_shape_change_not_drift(tmp_path):
    path = tmp_path / "drift.xyz"
    path.write_text(_tetrahedra([0.0, 0.5, 1.0]))
    positions = molframe.trajectory.read(path).positions
    raw = molframe.trajectory.rmsd(positions, align=False)
    fitted = molframe.trajectory.rmsd(positions, align=True)
    assert raw.value[0] == pytest.approx(0.0)
    assert raw.value[2] == pytest.approx(1.0)  # the body drifted by one ångström
    assert fitted.value[2] == pytest.approx(0.0, abs=1e-5)  # but kept its shape
    assert fitted.status == "complete" and fitted.profile == "molframe-default-1.0"
    assert not fitted.value.flags.writeable
    with pytest.raises(ValueError):
        molframe.trajectory.rmsd(positions, reference=9)
    with pytest.raises(ValueError):
        molframe.trajectory.rmsd(numpy.zeros((2, 3, 2), dtype=numpy.float32))
