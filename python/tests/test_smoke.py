"""Curated Python facade, ownership, and Workflow smoke tests."""

import subprocess
import sys
from pathlib import Path

import numpy
import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-py" / "tests" / "data"
IMPORT_BUDGET_SECONDS = 0.20


def test_import_stays_under_the_budget():
    source = (
        "import time; start = time.perf_counter(); "
        "import molframe; print(time.perf_counter() - start)"
    )
    elapsed = subprocess.run(
        [sys.executable, "-c", source],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    assert float(elapsed) < IMPORT_BUDGET_SECONDS


def test_root_is_curated_and_domains_use_final_names():
    expected = {
        "Structure",
        "Selection",
        "Query",
        "Reader",
        "Workflow",
        "CompiledWorkflow",
        "read",
        "geometry",
        "analysis",
        "trajectory",
        "sequence",
        "crystal",
        "validation",
        "motif",
        "chemistry",
        "formats",
    }
    assert expected <= set(molframe.__all__)
    for removed in ["Plan", "Batch", "AtomIndex", "StructureData", "geom", "traj", "seq"]:
        assert not hasattr(molframe, removed)


def test_read_selection_and_coordinate_ownership():
    structure = molframe.read(DATA / "basic.pdb")
    assert structure.atom_count == 2
    assert structure.residue_count == 1
    assert structure.chain_count == 1
    assert structure.model_count == 1
    assert len(structure.atoms) == 2
    assert len(structure.models) == 1
    assert len(structure.chains) == 1
    assert len(structure.residues) == 1
    assert structure.chains[0].label == "A"
    assert structure.chains["A"].residues[0].atoms[1].name == "CA"
    first = structure.coordinates
    second = structure.coordinates
    assert first.dtype == numpy.float32
    assert not first.flags.writeable
    assert numpy.shares_memory(first, second)

    selection = structure.select("name CA")
    compiled_query = molframe.Query("name CA")
    compiled = compiled_query.select(structure)
    through_structure = structure.select(compiled_query)
    assert len(selection) == 1
    assert selection.indices.tolist() == compiled.indices.tolist()
    assert compiled.indices.tolist() == through_structure.indices.tolist()
    assert selection.to_coordinates().shape == (1, 3)


def test_reader_reuses_owner_backed_input():
    payload = (DATA / "basic.cif").read_bytes()
    reader = molframe.Reader(payload, name="basic.cif")
    assert reader.byte_length == len(payload)
    assert len(reader.read().atoms) == 2
    assert len(reader.read().atoms) == 2


def test_non_contiguous_arrays_are_rejected_with_actionable_guidance():
    coordinates = numpy.zeros((4, 6), dtype=numpy.float32)[:, ::2]
    assert not coordinates.flags.c_contiguous
    with pytest.raises(ValueError, match="numpy.ascontiguousarray"):
        molframe.geometry.centroid(coordinates)


def test_contact_table_columns_are_aligned_and_owner_backed():
    structure = molframe.read(DATA / "basic.pdb")
    table = molframe.analysis.atom_contacts(structure, 3.0)
    assert len(table.first) == len(table.second) == len(table.distance) == len(table)
    assert table.first.dtype == numpy.uint32
    assert table.distance.dtype == numpy.float32
    first = table.first
    assert not first.flags.writeable
    del table
    assert first.shape[0] >= 0


def test_structure_edit_is_transactional():
    structure = molframe.read(DATA / "basic.pdb")
    editor = structure.edit()
    editor.rename_chain(0, "B")
    edited = editor.finish()
    assert len(edited.atoms) == len(structure.atoms)
    with pytest.raises(RuntimeError, match="already been finished"):
        editor.finish()


def test_workflow_compiles_explains_reuses_and_matches_eager():
    workflow = molframe.Workflow()
    mobile = workflow.input("mobile", kind="coordinates")
    reference = workflow.input("reference", kind="coordinates")
    score = workflow.rmsd(mobile, reference)
    workflow.output("score", score)
    compiled = workflow.compile()

    coordinates = numpy.ascontiguousarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=numpy.float32,
    )
    with pytest.raises(ValueError, match="copy=True"):
        compiled.run({"mobile": coordinates, "reference": coordinates})
    result = compiled.run(
        {"mobile": coordinates, "reference": coordinates},
        copy=True,
    )
    assert result["score"] == molframe.geometry.rmsd(coordinates, coordinates)
    explanation = compiled.explain()
    assert explanation["physical_node_count"] == 3
    assert [node["cost"] for node in explanation["nodes"]] == [
        "borrow",
        "borrow",
        "borrow",
    ]


def test_workflow_materializes_owner_backed_contact_tables():
    workflow = molframe.Workflow()
    structure = workflow.input("structure", kind="structure")
    contacts = workflow.atom_contacts(structure, 3.0)
    workflow.output("contacts", contacts)
    compiled = workflow.compile()
    result = compiled.run({"structure": molframe.read(DATA / "basic.pdb")})
    table = result["contacts"]
    assert len(table.first) == len(table.second) == len(table.distance) == len(table)
    assert compiled.explain()["spatial_consumers"] == 1
