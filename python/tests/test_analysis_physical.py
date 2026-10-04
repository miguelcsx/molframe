"""Physical analyses of one structure, each checked against an independent NumPy computation."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis

BENCH = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"

CIF_HEADER = (
    "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
    "_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
)


def structure_of(tmp_path: Path, rows: list[tuple[str, float, float, float]]) -> molframe.Structure:
    """Build a one-chain structure of single-atom residues from ``(element, x, y, z)`` rows."""
    lines = [
        f"ATOM {number} {element} {element}1 UNK A {number} {x} {y} {z}\n"
        for number, (element, x, y, z) in enumerate(rows, start=1)
    ]
    path = tmp_path / "t.cif"
    path.write_text(CIF_HEADER + "".join(lines))
    return molframe.read(path)


@pytest.fixture(scope="module")
def ubiquitin() -> molframe.Structure:
    return molframe.read(BENCH / "1ubq.cif")


def coordinates(structure: molframe.Structure) -> np.ndarray:
    return np.asarray(structure.coordinates, dtype=np.float64)


def indices(structure: molframe.Structure, query: str) -> np.ndarray:
    return np.asarray(structure.select(query).indices, dtype=np.int64)


def test_radial_distribution_of_one_group_matches_a_pair_histogram(ubiquitin):
    xyz = coordinates(ubiquitin)
    alpha = indices(ubiquitin, "name CA")
    low, high, bins, volume = 4.0, 24.0, 10, 5.0e4
    result = analysis.radial_distribution(
        ubiquitin,
        "name CA",
        "name CA",
        minimum_distance=low,
        maximum_distance=high,
        bins=bins,
        volume=volume,
    )
    gap = xyz[alpha][:, None, :] - xyz[alpha][None, :, :]
    distances = np.sqrt((gap**2).sum(axis=2))[np.triu_indices(len(alpha), k=1)]
    counts, edges = np.histogram(distances, bins=bins, range=(low, high))
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    pairs = len(alpha) * (len(alpha) - 1) / 2
    table = result.value
    assert result.status == "complete"
    assert table["count"].tolist() == counts.tolist()
    assert table["lower"] == pytest.approx(edges[:-1], abs=1e-4)
    assert table["upper"] == pytest.approx(edges[1:], abs=1e-4)
    assert table["distribution"] == pytest.approx(counts / (pairs * shell / volume), rel=1e-6)


def test_radial_distribution_between_two_groups_uses_every_cross_pair(ubiquitin):
    xyz = coordinates(ubiquitin)
    first, second = indices(ubiquitin, "name CA"), indices(ubiquitin, "resname HOH")
    table = analysis.radial_distribution(
        ubiquitin,
        "name CA",
        "resname HOH",
        minimum_distance=0.0,
        maximum_distance=12.0,
        bins=6,
        volume=1.0e4,
    ).value
    gap = xyz[first][:, None, :] - xyz[second][None, :, :]
    counts, _ = np.histogram(np.sqrt((gap**2).sum(axis=2)).ravel(), bins=6, range=(0.0, 12.0))
    assert table["count"].tolist() == counts.tolist()
    assert table["count"].sum() > 0


def test_coordination_numbers_count_the_neighbours_in_a_shell(ubiquitin):
    xyz = coordinates(ubiquitin)
    alpha, water = indices(ubiquitin, "name CA"), indices(ubiquitin, "resname HOH")
    low, high = 2.0, 8.0
    table = analysis.coordination_numbers(
        ubiquitin, "name CA", "resname HOH", minimum_distance=low, maximum_distance=high
    ).value
    gap = xyz[alpha][:, None, :] - xyz[water][None, :, :]
    distance = np.sqrt((gap**2).sum(axis=2))
    expected = ((distance >= low) & (distance <= high)).sum(axis=1)
    assert table["atom"].tolist() == alpha.tolist()
    assert table["count"].tolist() == expected.tolist()
    assert expected.sum() > 0


def test_coordination_within_one_group_counts_each_pair_at_both_ends(ubiquitin):
    xyz = coordinates(ubiquitin)
    alpha = indices(ubiquitin, "name CA")
    table = analysis.coordination_numbers(
        ubiquitin, "name CA", "name CA", minimum_distance=0.5, maximum_distance=6.0
    ).value
    gap = xyz[alpha][:, None, :] - xyz[alpha][None, :, :]
    distance = np.sqrt((gap**2).sum(axis=2))
    expected = ((distance >= 0.5) & (distance <= 6.0)).sum(axis=1)
    assert table["count"].tolist() == expected.tolist()


def test_linear_density_matches_a_weighted_histogram(ubiquitin):
    z = coordinates(ubiquitin)[:, 2]
    low, high, bins = float(z.min()) - 1.0, float(z.max()) + 1.0, 16
    count = analysis.linear_density(ubiquitin, axis="z", minimum=low, maximum=high, bins=bins).value
    expected, edges = np.histogram(z, bins=bins, range=(low, high))
    assert count["weight"] == pytest.approx(expected)
    assert count["density"] == pytest.approx(expected / np.diff(edges), rel=1e-5)

    symbol = {6: "C", 7: "N", 8: "O", 16: "S"}
    numbers = np.asarray(ubiquitin.elements)
    assert set(numbers.tolist()) <= symbol.keys()
    mass = np.array([molframe.chemistry.element(symbol[int(n)]).atomic_weight for n in numbers])
    by_mass = analysis.linear_density(
        ubiquitin, axis="z", minimum=low, maximum=high, bins=bins, weights="mass"
    ).value
    assert by_mass["weight"] == pytest.approx(
        np.histogram(z, bins=bins, range=(low, high), weights=mass)[0], rel=1e-6
    )

    charges = np.linspace(-1.0, 1.0, len(z))
    explicit = analysis.linear_density(
        ubiquitin, axis="z", minimum=low, maximum=high, bins=bins, weights=charges
    ).value
    assert explicit["weight"] == pytest.approx(
        np.histogram(z, bins=bins, range=(low, high), weights=charges)[0], abs=1e-9
    )


def test_linear_density_rejects_an_unknown_axis_and_unknown_weights(ubiquitin):
    with pytest.raises(molframe.MolframeError):
        analysis.linear_density(ubiquitin, axis="w", minimum=0.0, maximum=1.0, bins=2)
    with pytest.raises(ValueError, match="weights"):
        analysis.linear_density(
            ubiquitin, axis="x", minimum=0.0, maximum=1.0, bins=2, weights="charge"
        )


def test_density_map_matches_a_weighted_3d_histogram_and_reports_what_falls_outside(ubiquitin):
    xyz = coordinates(ubiquitin)
    origin = xyz.min(axis=0) + 2.0
    spacing, shape = (2.0, 2.5, 3.0), (10, 8, 6)
    grid = analysis.density_map(ubiquitin, origin=tuple(origin), spacing=spacing, shape=shape).value
    edges = [origin[axis] + spacing[axis] * np.arange(shape[axis] + 1) for axis in range(3)]
    counts, _ = np.histogramdd(xyz, bins=edges)
    voxel = float(np.prod(spacing))
    values = np.asarray(grid.values)
    assert values.shape == shape
    assert values == pytest.approx(counts / voxel, rel=1e-5, abs=1e-9)
    assert grid.excluded_weight == pytest.approx(len(xyz) - counts.sum())
    assert grid.excluded_weight > 0
    assert not values.flags.writeable


def test_leaflets_are_the_connected_components_of_the_site_graph(ubiquitin):
    xyz = coordinates(ubiquitin)
    alpha = indices(ubiquitin, "name CA")
    cutoff = 4.2
    table = analysis.leaflets(ubiquitin, "name CA", connection_distance=cutoff).value
    gap = xyz[alpha][:, None, :] - xyz[alpha][None, :, :]
    near = np.sqrt((gap**2).sum(axis=2)) <= cutoff
    label = -np.ones(len(alpha), dtype=int)
    for start in range(len(alpha)):
        if label[start] >= 0:
            continue
        label[start] = start
        stack = [start]
        while stack:
            current = stack.pop()
            for other in np.flatnonzero(near[current] & (label < 0)):
                label[other] = start
                stack.append(int(other))
    expected = {frozenset(alpha[label == value].tolist()) for value in set(label)}
    found: dict[int, set[int]] = {}
    for site, leaflet in zip(table["site"], table["leaflet"], strict=True):
        found.setdefault(int(leaflet), set()).add(int(site))
    assert {frozenset(group) for group in found.values()} == expected


def test_pore_profile_of_a_ring_is_the_ring_radius_less_atom_and_probe(tmp_path):
    ring = [
        ("C", 5.0 * np.cos(angle), 5.0 * np.sin(angle), z)
        for z in (0.0, 4.0)
        for angle in np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    ]
    structure = structure_of(tmp_path, [(e, round(x, 4), round(y, 4), z) for e, x, y, z in ring])
    probe = 0.5
    table = analysis.pore_profile(
        structure,
        axis_origin=(0.0, 0.0, 0.0),
        axis_direction=(0.0, 0.0, 1.0),
        start=0.0,
        end=4.0,
        samples=3,
        search_radius=1.0,
        grid_spacing=0.25,
        probe_radius=probe,
    ).value
    carbon = molframe.chemistry.vdw_radius("C")
    assert table["axial_coordinate"] == pytest.approx([0.0, 2.0, 4.0])
    # On a slice through a ring the axis is the widest point: ring radius - atom radius - probe.
    assert table["radius"][0] == pytest.approx(5.0 - carbon - probe, abs=0.05)
    assert table["radius"][2] == pytest.approx(5.0 - carbon - probe, abs=0.05)
    # Midway between the rings the nearest atom is farther, so the pore is wider.
    assert table["radius"][1] > table["radius"][0]
    assert np.hypot(table["centre_x"], table["centre_y"]).max() < 0.4


def test_surface_contacts_report_a_touching_pair_and_never_a_distant_atom(tmp_path):
    structure = structure_of(
        tmp_path,
        [("C", 0.0, 0.0, 0.0), ("C", 3.5, 0.0, 0.0), ("C", 50.0, 0.0, 0.0)],
    )
    carbon = molframe.chemistry.vdw_radius("C")
    assert carbon is not None
    options = {"probe": 1.4, "surface_density": 4.0, "minimum_area": 0.25}
    table = analysis.surface_contacts(structure, tolerance=0.2, **options).value
    assert (table["first"].tolist(), table["second"].tolist()) == ([0], [1])
    assert table["distance"][0] == pytest.approx(3.5)
    # 3.5 Å exceeds the 3.4 Å sum of the two radii, so a zero tolerance leaves no contact.
    assert 2 * carbon < 3.5
    assert len(analysis.surface_contacts(structure, tolerance=0.0, **options).value) == 0


def test_surface_contacts_of_a_protein_are_close_pairs_in_atom_order(ubiquitin):
    xyz = coordinates(ubiquitin)
    radii = np.asarray(molframe.chemistry.vdw_radii(ubiquitin), dtype=np.float64)
    tolerance = 0.5
    result = analysis.surface_contacts(
        ubiquitin, tolerance=tolerance, probe=1.4, surface_density=10.0, minimum_area=0.1
    )
    table = result.value
    first, second = table["first"].astype(int), table["second"].astype(int)
    assert result.status == "complete"
    assert len(table) > 0
    assert (first < second).all()
    separation = np.linalg.norm(xyz[first] - xyz[second], axis=1)
    assert table["distance"] == pytest.approx(separation, abs=1e-3)
    assert (separation <= radii[first] + radii[second] + tolerance + 1e-4).all()
    # A contact that is exposed is exposed at a tolerance that admits more pairs.
    looser = analysis.surface_contacts(
        ubiquitin, tolerance=1.5, probe=1.4, surface_density=10.0, minimum_area=0.1
    ).value
    assert len(looser) >= len(table)


def test_selections_accept_text_queries_and_selections_of_the_same_structure(ubiquitin):
    by_text = analysis.coordination_numbers(
        ubiquitin, "name CA", "name CA", minimum_distance=1.0, maximum_distance=5.0
    ).value
    by_selection = analysis.coordination_numbers(
        ubiquitin,
        ubiquitin.select("name CA"),
        ubiquitin.select("name CA"),
        minimum_distance=1.0,
        maximum_distance=5.0,
    ).value
    assert by_text["count"].tolist() == by_selection["count"].tolist()
    crambin = molframe.read(BENCH / "1crn.cif")
    with pytest.raises(ValueError, match="different structure"):
        analysis.coordination_numbers(
            crambin,
            ubiquitin.select("all"),
            "name CA",
            minimum_distance=1.0,
            maximum_distance=5.0,
        )
