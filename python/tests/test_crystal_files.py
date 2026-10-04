"""Space groups, density maps and reflection tables, each checked against what it must equal."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import chemistry, crystal

BENCH = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"
MOL = """ethanol
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"""


def test_a_space_group_is_found_by_name_and_lists_its_operations():
    group = crystal.SpaceGroup.from_hermann_mauguin("P 21 21 21")
    assert group.international_number == 19
    assert len(group) == 4
    expressions = {operation.expression for operation in group.operations}
    assert expressions == {"x,y,z", "1/2-x,-y,1/2+z", "-x,1/2+y,1/2-z", "1/2+x,1/2-y,-z"}
    # The same group by its Hall symbol and by the catalogue's own spelling.
    by_hall = crystal.SpaceGroup.from_hall(group.hall_symbol)
    assert by_hall.hall_number == group.hall_number
    assert crystal.SpaceGroup.from_hermann_mauguin("p_2_1_2_1_2_1").hall_number == group.hall_number


def test_the_operations_act_on_fractional_coordinates_and_are_exact():
    group = crystal.SpaceGroup.from_hermann_mauguin("P 21 21 21")
    point = [0.1, 0.2, 0.3]
    images = {
        tuple(round(value, 12) for value in operation.apply(point))
        for operation in group.operations
    }
    assert images == {
        (0.1, 0.2, 0.3),
        (0.4, -0.2, 0.8),
        (-0.1, 0.7, 0.2),
        (0.6, 0.3, -0.3),
    }
    identity = [operation for operation in group.operations if operation.is_identity]
    assert len(identity) == 1
    screw = next(op for op in group.operations if op.expression == "1/2-x,-y,1/2+z")
    assert screw.rotation == ((-1, 0, 0), (0, -1, 0), (0, 0, 1))
    assert screw.translation == (0.5, 0.0, 0.5)


def test_every_catalogued_type_has_a_setting_and_the_operation_counts_follow_the_lattice():
    assert len(crystal.SpaceGroup.settings(1)) == 1
    assert len(crystal.SpaceGroup.settings(5)) > 1
    assert len(crystal.SpaceGroup.settings(230)) >= 1
    # P1 has one operation; the primitive cubic P m -3 m has 48.
    assert len(crystal.SpaceGroup.from_hermann_mauguin("P 1")) == 1
    assert len(crystal.SpaceGroup.from_hermann_mauguin("P m -3 m")) == 48
    with pytest.raises(molframe.MolframeError) as unknown:
        crystal.SpaceGroup.from_hermann_mauguin("P 99")
    assert unknown.value.code == "MOLFRAME-E6018"
    with pytest.raises(molframe.MolframeError):
        crystal.SpaceGroup.settings(231)


def test_a_structure_reports_the_space_group_it_carries():
    group = crystal.space_group(molframe.read(BENCH / "1crn.cif"))
    assert group is not None
    assert group.international_number == 4
    assert len(group) == 2
    ethanol = chemistry.read_mol(MOL).to_structure()
    assert crystal.space_group(ethanol) is None


def grid():
    rng = np.random.default_rng(4)
    return rng.normal(size=(3, 4, 5)).astype(np.float32)


def test_a_map_round_trips_through_mrc_and_keeps_its_geometry(tmp_path):
    values = grid()
    cell = crystal.UnitCell([10.0, 20.0, 30.0], [90.0, 90.0, 90.0])
    made = crystal.DensityMap(values, cell, space_group=1, labels=["test map"])
    assert made.dimensions == (5, 4, 3)
    assert made.sampling == (5, 4, 3)
    path = tmp_path / "map.mrc"
    path.write_bytes(made.to_mrc())
    read = crystal.read_mrc(path)
    assert np.array_equal(read.values, values)
    assert read.dimensions == (5, 4, 3)
    assert read.cell == cell
    assert read.space_group == 1
    assert read.labels == ["test map"]
    assert np.array_equal(crystal.read_mrc(path.read_bytes()).values, values)
    with pytest.raises(molframe.MolframeError) as truncated:
        crystal.read_mrc(path.read_bytes()[:500])
    assert truncated.value.code in {"MOLFRAME-E7101", "MOLFRAME-E1201"}


def test_map_statistics_and_histogram_equal_numpy():
    values = grid()
    made = crystal.DensityMap(values, crystal.UnitCell([5.0, 5.0, 5.0], [90.0, 90.0, 90.0]))
    found = made.statistics()
    flat = values.astype(np.float64).ravel()
    assert found.count == flat.size
    assert found.mean == pytest.approx(flat.mean())
    assert found.sigma == pytest.approx(flat.std())
    assert (found.minimum, found.maximum) == (float(values.min()), float(values.max()))
    counts = made.histogram(6, -2.0, 2.0)
    expected, _ = np.histogram(flat, bins=6, range=(-2.0, 2.0))
    assert np.array_equal(counts, expected)
    mask = values > 0
    masked = made.masked_statistics(mask)
    assert masked.count == int(mask.sum())
    assert masked.mean == pytest.approx(flat[mask.ravel()].mean())
    with pytest.raises(molframe.MolframeError):
        made.masked_statistics(np.zeros_like(mask))
    with pytest.raises(molframe.MolframeError):
        made.histogram(0, 0.0, 1.0)


def test_sampling_a_map_at_its_grid_points_returns_the_grid_values():
    values = grid()
    cell = crystal.UnitCell([10.0, 20.0, 30.0], [90.0, 90.0, 90.0])
    made = crystal.DensityMap(values, cell)
    nz, ny, nx = values.shape
    points, expected = [], []
    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                points.append([i * 10.0 / nx, j * 20.0 / ny, k * 30.0 / nz])
                expected.append(values[k, j, i])
    found = made.sample(np.array(points))
    assert np.allclose(found, expected, atol=1e-5)
    assert np.allclose(made.sample(np.array(points), method="cubic"), expected, atol=1e-4)
    midpoint = np.array([[0.5 * 10.0 / nx, 0.0, 0.0]])
    assert float(made.sample(midpoint)[0]) == pytest.approx(
        0.5 * (values[0, 0, 0] + values[0, 0, 1]), abs=1e-5
    )
    outside = np.array([[-5.0, 0.0, 0.0]])
    assert np.isnan(made.sample(outside)[0])
    wrapped = made.sample(np.array([[10.0 + points[0][0], 0.0, 0.0]]), boundary="periodic")
    assert float(wrapped[0]) == pytest.approx(float(values[0, 0, 0]), abs=1e-5)
    with pytest.raises(molframe.MolframeValueError, match="not a method"):
        made.sample(np.array(points), method="spline")  # type: ignore[arg-type]


def test_a_reflection_table_round_trips_through_mtz_with_its_missing_values():
    hkl = np.array([[1, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 3]], dtype=np.int32)
    table = crystal.ReflectionTable(
        hkl,
        {
            "FP": ("amplitude", [10.0, 20.0, float("nan"), 40.0]),
            "SIGFP": ("standard_deviation", [1.0, 2.0, 3.0, 4.0]),
            "FREE": ("flag", [0.0, 1.0, 0.0, 1.0]),
        },
        cell=crystal.UnitCell([30.0, 40.0, 50.0], [90.0, 90.0, 90.0]),
        space_group_number=19,
        space_group_name="P 21 21 21",
        title="synthetic",
    )
    assert table.row_count == 4
    assert table.types["FP"] == "amplitude"
    again = crystal.read_mtz(table.to_mtz())
    assert again.row_count == 4
    assert np.array_equal(again.miller_indices, hkl)
    found = np.asarray(again.column("FP"), dtype=np.float64)
    assert np.allclose(found[[0, 1, 3]], [10.0, 20.0, 40.0])
    assert np.isnan(found[2])
    assert np.allclose(again.column("SIGFP"), [1.0, 2.0, 3.0, 4.0])
    assert again.cell == table.cell
    assert again.space_group_number == 19
    assert {"H", "K", "L", "FP", "SIGFP", "FREE"} <= set(again.labels)
    with pytest.raises(molframe.MolframeKeyError):
        again.column("missing")


def test_a_reflection_table_refuses_columns_that_do_not_match_its_indices():
    hkl = np.array([[1, 0, 0], [0, 1, 0]], dtype=np.int32)
    with pytest.raises(molframe.MolframeValueError, match="values for 2 reflections"):
        crystal.ReflectionTable(hkl, {"FP": ("amplitude", [1.0])})
    with pytest.raises(molframe.MolframeError, match="not a reflection column type"):
        crystal.ReflectionTable(hkl, {"FP": ("weights", [1.0, 2.0])})
    with pytest.raises(molframe.MolframeError) as invalid:
        crystal.read_mtz(b"not an mtz file")
    assert invalid.value.code in {"MOLFRAME-E1201", "MOLFRAME-E7101"}
