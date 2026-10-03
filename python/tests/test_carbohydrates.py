import gzip
import importlib
import math
import shlex
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-chem" / "data"
CCD = DATA / "CCD-saccharides.cif"


def _category(source, name):
    """Read the simple single-line loops used by the deposited fixture oracle."""
    lines = source.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(f"_{name}."))
    headers = []
    while lines[start].startswith(f"_{name}."):
        headers.append(lines[start].strip().split(".", 1)[1])
        start += 1
    rows = []
    while start < len(lines) and not lines[start].startswith("#"):
        values = shlex.split(lines[start])
        assert len(values) == len(headers)
        rows.append(dict(zip(headers, values, strict=True)))
        start += 1
    return rows


@pytest.fixture
def deposited():
    source = gzip.decompress((DATA / "1HZH.cif.gz").read_bytes()).decode()
    structure = molframe.read(DATA / "1HZH.cif.gz")
    atoms = _category(source, "atom_site")
    assert len(atoms) == structure.atom_count
    assert all(structure.atoms[i].name == row["label_atom_id"] for i, row in enumerate(atoms))
    return source, structure, atoms


def test_all_carbohydrate_public_names_are_registered_and_namespace_identity_is_stable():
    chemistry = importlib.import_module("molframe.chemistry")
    assert molframe.chemistry is chemistry
    for name in (
        "SnfgSymbol",
        "RingGeometry",
        "Monosaccharide",
        "CarbohydrateLink",
        "CarbohydrateReport",
        "snfg_symbol",
        "carbohydrates",
    ):
        assert getattr(chemistry, name) is not None
    symbol = chemistry.snfg_symbol("NAG")
    assert isinstance(symbol, chemistry.SnfgSymbol)
    assert (symbol.abbreviation, symbol.name, symbol.shape) == (
        "GlcNAc",
        "N-Acetyl Glucosamine",
        "filled_cube",
    )
    assert symbol.color == 0x0090BC
    assert symbol.secondary_color is None
    assert chemistry.snfg_symbol("MAN").abbreviation == chemistry.snfg_symbol("BMA").abbreviation
    assert chemistry.snfg_symbol("GalN").secondary_color == 0xF1ECE1
    assert chemistry.snfg_symbol("ZZZ") is None
    with pytest.raises(AttributeError):
        symbol.color = 0


def _link_keys(links, atoms):
    return {
        (
            atoms[link.donor_atom]["label_asym_id"],
            int(atoms[link.donor_atom]["auth_seq_id"]),
            atoms[link.donor_atom]["label_atom_id"],
            int(atoms[link.acceptor_atom]["auth_seq_id"]),
            atoms[link.acceptor_atom]["label_atom_id"],
        )
        for link in links
    }


def _branch_keys(source):
    scheme = _category(source, "pdbx_branch_scheme")
    residues = {(row["entity_id"], row["num"]): row for row in scheme}
    expected = set()
    for row in _category(source, "pdbx_entity_branch_link"):
        donor = residues[row["entity_id"], row["entity_branch_list_num_1"]]
        acceptor = residues[row["entity_id"], row["entity_branch_list_num_2"]]
        assert donor["asym_id"] == acceptor["asym_id"]
        expected.add(
            (
                donor["asym_id"],
                int(donor["pdb_seq_num"]),
                row["atom_id_1"],
                int(acceptor["pdb_seq_num"]),
                row["atom_id_2"],
            )
        )
    return expected


def test_deposited_inventory_geometry_and_all_links_match_raw_cif_without_mutation(deposited):
    source, structure, atoms = deposited
    coordinates = np.array(structure.coordinates, copy=True)
    bond_count = structure.bond_count
    report = molframe.chemistry.carbohydrates(structure, CCD, version="wwPDB-2026-10-02")
    assert isinstance(report, molframe.chemistry.CarbohydrateReport)
    assert report.dictionary_version == "wwPDB-2026-10-02"
    assert len(report.monosaccharides) == 18
    assert report.incomplete_residues == report.terminal_links == []
    assert Counter(sugar.symbol.abbreviation for sugar in report.monosaccharides) == {
        "GlcNAc": 8,
        "Man": 6,
        "Gal": 3,
        "Fuc": 1,
    }
    observed = set()
    for sugar in report.monosaccharides:
        assert isinstance(sugar, molframe.chemistry.Monosaccharide)
        assert len(sugar.ring_atoms) == 6
        assert sugar.anomeric_atom in sugar.ring_atoms
        assert structure.residues[sugar.residue].name == atoms[sugar.ring_atoms[0]]["label_comp_id"]
        row = atoms[sugar.ring_atoms[0]]
        observed.add((row["label_asym_id"], int(row["auth_seq_id"]), row["label_comp_id"]))
        geometry = sugar.geometry
        assert isinstance(geometry, molframe.chemistry.RingGeometry)
        assert all(math.isfinite(value) for value in (*geometry.center, *geometry.normal))
        assert sum(value * value for value in geometry.normal) == pytest.approx(1.0, abs=1e-6)
        expected_center = coordinates[sugar.ring_atoms].mean(axis=0, dtype=np.float64)
        assert geometry.center == pytest.approx(expected_center, abs=1e-5)
        assert geometry.anomeric_direction is not None
        assert sum(value * value for value in geometry.anomeric_direction) == pytest.approx(1.0)
    expected = {
        (row["asym_id"], int(row["pdb_seq_num"]), row["mon_id"])
        for row in _category(source, "pdbx_branch_scheme")
    }
    assert observed == expected
    assert len(report.links) == 16
    assert _link_keys(report.links, atoms) == _branch_keys(source)
    for link in report.links:
        assert isinstance(link, molframe.chemistry.CarbohydrateLink)
        assert link.provenance == "file"
        assert link.acceptor is not None
        assert (
            structure.atoms[link.donor_atom].residue.index
            == report.monosaccharides[link.donor].residue
        )
        assert (
            structure.atoms[link.acceptor_atom].residue.index
            == report.monosaccharides[link.acceptor].residue
        )
    assert structure.bond_count == bond_count
    np.testing.assert_array_equal(structure.coordinates, coordinates)
    copied_rings = report.monosaccharides[0].ring_atoms
    copied_rings.clear()
    assert len(report.monosaccharides[0].ring_atoms) == 6
    with pytest.raises(AttributeError):
        report.dictionary_version = "changed"


def test_removed_declarations_recover_every_branch_link_with_inferred_provenance(
    deposited, tmp_path
):
    source, _, atoms = deposited
    # Remove only the deposited sugar/sugar struct_conn rows, retaining all other categories.
    conn = _category(source, "struct_conn")
    sugars = {"NAG", "MAN", "BMA", "GAL", "FUC"}
    assert not any(
        (row["ptnr1_label_comp_id"] in sugars) != (row["ptnr2_label_comp_id"] in sugars)
        for row in conn
    )
    lines = source.splitlines(keepends=True)
    start = next(i for i, line in enumerate(lines) if line.startswith("_struct_conn."))
    end = next(i for i in range(start, len(lines)) if lines[i].startswith("#"))
    headers = []
    while lines[start].startswith("_struct_conn."):
        headers.append(lines[start].strip().split(".", 1)[1])
        start += 1
    retained = [
        line
        for line in lines[start:end]
        if shlex.split(line)[headers.index("ptnr1_label_comp_id")] not in sugars
    ]
    path = tmp_path / "without-sugar-links.cif"
    path.write_text("".join(lines[:start] + retained + lines[end:]))
    structure = molframe.read(path)
    disabled = molframe.chemistry.carbohydrates(structure, CCD, spatial_fallback=False)
    assert len(disabled.monosaccharides) == 18
    # The public reader already inferred these valid edges; disabling new search
    # preserves input topology rather than erasing its provenance.
    assert _link_keys(disabled.links, atoms) == _branch_keys(source)
    assert all(link.provenance == "inferred_distance" for link in disabled.links)
    report = molframe.chemistry.carbohydrates(structure, str(CCD))
    assert len(report.links) == 16
    assert _link_keys(report.links, atoms) == _branch_keys(source)
    assert all(link.provenance == "inferred_distance" for link in report.links)
    assert report.terminal_links == []


def _synthetic(tmp_path, *, distance=1.4, occupied=False, missing=False, **options: object):
    name = options.get("name", "ND2")
    flat = options.get("flat", False)
    names = ["O5", "C1", "C2", "C3", "C4", "C5"]
    rows = []
    for i, atom_name in enumerate(names):
        if missing and atom_name == "C3":
            continue
        angle = (i - 1) * math.tau / 6
        x, y = (0, 0) if flat else (1.5 * math.cos(angle), 1.5 * math.sin(angle))
        rows.append(
            f"HETATM {len(rows) + 1} {'O' if i == 0 else 'C'} {atom_name} . NAG A 1 {x} {y} 0"
        )
    if occupied:
        rows.append(f"HETATM {len(rows) + 1} O O1 . NAG A 1 2.8 0 0")
    rows.extend(
        (
            f"ATOM {len(rows) + 1} N {name} . ASN B 1 {1.5 + distance} 0 0",
            f"ATOM {len(rows) + 2} C CG . ASN B 1 {2.8 + distance} 0 0",
        )
    )
    path = tmp_path / "synthetic.cif"
    fields = (
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_seq_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
    )
    path.write_text(
        "data_test\nloop_\n"
        + "\n".join(f"_atom_site.{item}" for item in fields)
        + "\n"
        + "\n".join(rows)
        + "\n"
    )
    return molframe.read(path)


@pytest.mark.parametrize(("distance", "expected"), [(1.4, 1), (2.0, 1), (2.001, 0), (0.8, 0)])
def test_terminal_distance_inference_has_conservative_boundary(tmp_path, distance, expected):
    structure = _synthetic(tmp_path, distance=distance)
    report = molframe.chemistry.carbohydrates(structure, CCD)
    assert len(report.terminal_links) == expected
    assert len(report.monosaccharides) == 1
    if expected:
        link = report.terminal_links[0]
        assert link.acceptor is None
        assert link.provenance == "inferred_distance"
        assert structure.atoms[link.acceptor_atom].name == "ND2"
    disabled = molframe.chemistry.carbohydrates(structure, CCD, spatial_fallback=False)
    # At 1.4 Å the public reader already supplied a valid distance-inferred edge.
    assert len(disabled.terminal_links) == (1 if distance == 1.4 else 0)


@pytest.mark.parametrize("options", [{"occupied": True}, {"name": "NOT_A_SITE"}])
def test_occupied_donors_and_unknown_acceptor_chemistry_are_not_linked(tmp_path, options):
    report = molframe.chemistry.carbohydrates(_synthetic(tmp_path, **options), CCD)
    assert report.terminal_links == []


def test_missing_ring_topology_and_degenerate_geometry_remain_explicit(tmp_path):
    missing = molframe.chemistry.carbohydrates(_synthetic(tmp_path, missing=True), CCD)
    assert missing.monosaccharides == []
    assert missing.incomplete_residues == [0]
    report = molframe.chemistry.carbohydrates(_synthetic(tmp_path, flat=True), CCD)
    assert len(report.monosaccharides) == 1
    assert report.monosaccharides[0].geometry is None
    bare = molframe.chemistry.carbohydrates(_synthetic(tmp_path))
    assert bare.dictionary_version is None
    assert len(bare.monosaccharides) == 1
    assert bare.incomplete_residues == []


def test_dictionary_and_structure_arguments_fail_at_the_public_boundary(tmp_path):
    structure = _synthetic(tmp_path)
    with pytest.raises(OSError, match="MOLFRAME-E7101"):
        molframe.chemistry.carbohydrates(structure, tmp_path / "missing.cif")
    with pytest.raises(TypeError):
        molframe.chemistry.carbohydrates("not a structure")
