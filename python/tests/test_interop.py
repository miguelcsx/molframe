import ctypes
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.ipc
import pyarrow.parquet
import pytest

import molframe
from molframe import interop

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "crates" / "molframe-py" / "tests" / "data"
FOUR_HHB = ROOT / "crates" / "molframe-bench" / "data" / "4hhb.cif"
CCD = ROOT / "crates" / "molframe-chem" / "data" / "CCD-amino-acids.cif"


@pytest.fixture(scope="module")
def hemoglobin():
    return molframe.read(FOUR_HHB)


def test_the_hierarchy_streams_whole_tables_to_arrow(hemoglobin):
    atoms = pa.table(hemoglobin.atoms)
    assert atoms.num_rows == hemoglobin.atom_count
    assert atoms.column("atom_index").to_pylist()[:3] == [0, 1, 2]
    assert pa.table(hemoglobin.residues).num_rows == hemoglobin.residue_count
    assert pa.table(hemoglobin.chains).num_rows == hemoglobin.chain_count
    assert pa.table(hemoglobin.bonds).num_rows == hemoglobin.bond_count == len(hemoglobin.bonds)
    assert pa.table(hemoglobin.bonds).schema.names == [
        "bond_index",
        "atom_a",
        "atom_b",
        "order",
        "provenance",
    ]


def test_a_view_streams_only_its_own_rows(hemoglobin):
    residue = hemoglobin.residues[10]
    rows = pa.table(residue.atoms)
    assert rows.num_rows == len(residue.atoms)
    first = residue.atoms[0].index
    assert rows.column("atom_index").to_pylist() == list(range(first, first + len(residue.atoms)))
    chain = hemoglobin.chains[1]
    chain_residues = pa.table(chain.residues)
    assert chain_residues.num_rows == len(chain.residues)
    assert chain_residues.column("chain_index").to_pylist() == [chain.index] * len(chain.residues)


def test_coordinates_arrive_as_a_fixed_size_float_list(hemoglobin):
    column = pa.table(hemoglobin.atoms).column("coordinates")
    assert pa.types.is_fixed_size_list(column.type)
    assert column.type.list_size == 3
    assert column.type.value_type == pa.float32()
    flat = np.asarray(column.combine_chunks().flatten()).reshape(-1, 3)
    assert np.array_equal(flat, np.asarray(hemoglobin.coordinates))


def test_a_result_table_streams_with_its_own_column_types(hemoglobin):
    rules = [
        molframe.chemistry.PolymerRoleRule(name, 1 << bit, component_kind=1)
        for bit, name in enumerate(("N", "CA", "C", "O"))
    ]
    annotated = molframe.chemistry.apply_polymer_role_profile(
        hemoglobin, CCD, rules, profile_id="backbone", version="wwPDB-2026-10-03"
    ).structure
    table = molframe.analysis.dssp(annotated)
    arrow = pa.table(table)
    assert arrow.column_names == table.names
    assert arrow.num_rows == len(table)
    for name in table.names:
        assert arrow.column(name).to_pylist() == list(table[name])


def test_a_contact_table_still_streams(hemoglobin):
    contacts = molframe.analysis.atom_contacts(hemoglobin, 3.0)
    assert pa.table(contacts).num_rows == len(contacts)


def test_coordinates_import_through_dlpack_as_an_independent_copy(hemoglobin):
    tensor = interop.coordinates(hemoglobin)
    assert tensor.shape == (hemoglobin.atom_count, 3)
    assert tensor.__dlpack_device__() == (1, 0)
    imported = np.from_dlpack(tensor)
    assert imported.dtype == np.float32
    assert np.array_equal(imported, np.asarray(hemoglobin.coordinates))
    editable = imported.copy()
    editable[0, 0] = 1.0e6
    assert np.asarray(hemoglobin.coordinates)[0, 0] != 1.0e6
    second = np.from_dlpack(tensor)
    assert second.ctypes.data != imported.ctypes.data


def test_a_dlpack_capsule_nobody_takes_is_freed_not_leaked(hemoglobin):
    tensor = interop.coordinates(hemoglobin)
    capsule = tensor.__dlpack__()
    name = ctypes.pythonapi.PyCapsule_GetName
    name.restype = ctypes.c_char_p
    name.argtypes = [ctypes.py_object]
    assert name(capsule) == b"dltensor"
    del capsule  # the destructor frees it; a double free would crash here
    taken = np.from_dlpack(tensor)
    assert taken.shape[1] == 3


def test_dlpack_refuses_a_device_it_cannot_serve(hemoglobin):
    with pytest.raises(BufferError, match="CPU"):
        interop.coordinates(hemoglobin).__dlpack__(dl_device=(2, 0))


def test_a_graph_carries_pyg_conventions(hemoglobin):
    graph = interop.graph(
        hemoglobin,
        nodes="residues",
        edges="contacts",
        cutoff=8.0,
        direction="symmetric",
        node_features=["atom_count", "b_factor"],
        edge_features=["distance"],
    )
    assert graph.node_count == hemoglobin.residue_count
    assert graph.edge_index.shape == (2, graph.edge_count)
    assert str(graph.edge_index.dtype) == "int64"
    assert graph.node_feature_names == ["atom_count", "b_factor"]
    assert graph.node_features.shape == (graph.node_count, 2)
    assert graph.edge_features.shape == (graph.edge_count, 1)
    sources, targets = graph.edge_index
    assert sources.min() >= 0
    assert targets.max() < graph.node_count
    pairs = set(zip(sources.tolist(), targets.tolist(), strict=True))
    assert all((b, a) in pairs for a, b in pairs), "symmetric graphs store both orientations"
    assert float(graph.edge_features.max()) <= 8.0 + 1e-3
    assert np.array_equal(
        graph.node_features[:, 0], [len(residue.atoms) for residue in hemoglobin.residues]
    )


def test_a_graph_states_every_decision_or_is_refused(hemoglobin):
    with pytest.raises(TypeError):
        interop.graph(hemoglobin, nodes="atoms", edges="bonds")  # type: ignore[call-arg]
    with pytest.raises(
        molframe.PolicyError, match=r"edges: contacts .expected contacts with a cutoff"
    ):
        interop.graph(hemoglobin, nodes="atoms", edges="contacts", direction="undirected")
    with pytest.raises(molframe.PolicyError, match="node_feature: charm"):
        interop.graph(
            hemoglobin,
            nodes="atoms",
            edges="bonds",
            direction="undirected",
            node_features=["charm"],
        )
    with pytest.raises(molframe.MolframeError) as unsupported:
        interop.graph(
            hemoglobin,
            nodes="residues",
            edges="contacts",
            cutoff=4.0,
            direction="undirected",
            node_features=["element"],
        )
    assert unsupported.value.code == "MOLFRAME-E6103"


def test_a_missing_feature_is_refused_unless_a_fill_value_is_stated(hemoglobin):
    request = {
        "nodes": "atoms",
        "edges": "bonds",
        "direction": "undirected",
        "node_features": ["partial_charge"],
    }
    with pytest.raises(molframe.MolframeError, match="partial_charge"):
        interop.graph(hemoglobin, **request)
    filled = interop.graph(hemoglobin, missing=0.0, **request)
    assert not np.any(filled.node_features)


def test_atom_tables_round_trip_through_parquet_and_ipc(tmp_path, hemoglobin):
    parquet = tmp_path / "atoms.parquet"
    ipc = tmp_path / "atoms.arrow"
    interop.write_atoms_parquet(hemoglobin, parquet, metadata={"source": "4hhb"})
    interop.write_atoms_ipc(hemoglobin, ipc, metadata={"source": "4hhb"})
    expected = pa.table(hemoglobin.atoms)
    from_parquet = pyarrow.parquet.read_table(parquet)
    assert from_parquet.num_rows == expected.num_rows
    assert (
        from_parquet.column("coordinates").to_pylist() == expected.column("coordinates").to_pylist()
    )
    assert from_parquet.schema.metadata[b"source"] == b"4hhb"
    with pyarrow.ipc.open_file(ipc) as reader:
        from_ipc = reader.read_all()
    assert from_ipc.num_rows == expected.num_rows
    assert from_ipc.schema.metadata[b"source"] == b"4hhb"


def test_a_failed_write_names_the_output_code(tmp_path, hemoglobin):
    with pytest.raises(OSError, match="table file I/O failed") as raised:
        interop.write_atoms_parquet(hemoglobin, tmp_path / "missing" / "atoms.parquet")
    assert raised.value.code == "MOLFRAME-E7901"


def test_a_manifest_filters_splits_and_batches_without_loading_structures():
    dataset = interop.Dataset.from_manifest(DATA / "dataset.json")
    assert len(dataset) == 2
    first = dataset[0]
    assert first.id == "basic-a"
    assert first.resolution == pytest.approx(1.8)
    sharp = dataset.filter(resolution_below=2.0)
    assert 0 < len(sharp) < len(dataset)
    assert all(entry.resolution < 2.0 for entry in sharp.entries)
    split = dataset.split("temporal", ratios=(0.5, 0.25, 0.25))
    assert len(split.train) + len(split.validation) + len(split.test) == len(dataset)
    assert split.warnings == []
    batches = dataset.batches(2)
    assert sum(len(batch) for batch in batches) == len(dataset)
    assert max(len(batch) for batch in batches) <= 2
    loaded = sharp.load(0)
    assert loaded.atom_count == sharp[0].atom_count


def test_a_random_split_needs_a_seed_and_says_what_it_cannot_promise():
    dataset = interop.Dataset.from_manifest(DATA / "dataset.json")
    with pytest.raises(molframe.PolicyError, match="random with a seed"):
        dataset.split("random", ratios=(0.5, 0.25, 0.25))
    split = dataset.split("random", ratios=(0.5, 0.25, 0.25), seed=3)
    again = dataset.split("random", ratios=(0.5, 0.25, 0.25), seed=3)
    assert [e.id for e in split.train.entries] == [e.id for e in again.train.entries]
    assert len(split.warnings) == 1
    assert "redundancy" in split.warnings[0]


def test_datasets_are_built_from_entries_and_refuse_duplicates(tmp_path):
    entries = [
        interop.ManifestEntry("a", DATA / "basic.pdb", 2, tags=["x"], statistics={"n": 1.0}),
        interop.ManifestEntry("b", DATA / "basic.pdb", 2),
    ]
    dataset = interop.Dataset(entries)
    assert [entry.id for entry in dataset.entries] == ["a", "b"]
    assert dataset[-1].id == "b"
    assert dataset[0].tags == ["x"]
    assert dataset[0].statistics == {"n": 1.0}
    with pytest.raises(IndexError):
        dataset[2]
    with pytest.raises(molframe.MolframeError) as duplicate:
        interop.Dataset([entries[0], entries[0]])
    assert duplicate.value.code == "MOLFRAME-E7101"
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"entries": [{"id": "", "path": "x", "atom_count": 1}]}))
    with pytest.raises(molframe.MolframeError):
        interop.Dataset.from_manifest(manifest)
