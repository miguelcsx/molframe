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


def test_version_is_the_installed_distribution_version():
    from importlib.metadata import version

    assert molframe.__version__ == version("molframe")


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
    atom = structure.atoms[1]
    assert atom.residue is not None
    assert atom.residue.index == 0
    assert atom.residue.name == "GLY"
    residue = structure.residues[0]
    assert residue.atom("CA").index == atom.index
    assert residue.atom("missing") is None
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
    residues = selection.residues()
    assert len(residues) == 1
    assert residues[0].index == 0

    editor = structure.edit()
    editor.rename_chain(0, "B")
    mismatched = editor.finish()
    with pytest.raises(ValueError, match="stale"):
        selection.residues(mismatched)


def test_atom_component_name_preserves_native_chemistry_identity():
    structure = molframe.read(DATA / "basic.pdb")

    assert structure.atoms[0].component_name == "GLY"
    assert structure.atoms[1].component_name == "GLY"
def test_reader_reuses_owner_backed_input():
    payload = (DATA / "basic.cif").read_bytes()
    reader = molframe.Reader(payload, name="basic.cif")
    assert reader.byte_length == len(payload)
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


def test_named_queries_resolve_to_closed_queries():
    structure = molframe.read(DATA / "basic.pdb")
    aliases = molframe.QueryAliases()
    aliases.define("first", molframe.Query("index 0"))
    aliases.define("both", molframe.Query("$first or index 1"))
    assert aliases.names == ["both", "first"]
    assert "first" in aliases and len(aliases) == 2
    reference = molframe.Query("$both and not $first")
    assert reference.references == ["both", "first"]
    closed = aliases.resolve(reference)
    assert closed.references == []
    assert structure.select(closed).indices.tolist() == [1]
    assert closed.fingerprint == molframe.Query(closed.source).fingerprint


def test_named_query_failures_are_value_errors():
    aliases = molframe.QueryAliases()
    aliases.define("loop", molframe.Query("$loop"))
    with pytest.raises(ValueError, match="E4006"):
        aliases.resolve(molframe.Query("$loop"))
    with pytest.raises(ValueError, match="E4005"):
        aliases.resolve(molframe.Query("$missing"))
    with pytest.raises(ValueError):
        aliases.define("1bad", molframe.Query("all"))


def test_select_accepts_text_and_compiled_queries_alike():
    structure = molframe.read(DATA / "basic.pdb")
    assert (
        structure.select("index 1").indices.tolist()
        == structure.select(molframe.Query("index 1")).indices.tolist()
        == structure.select(molframe.sel.all() & molframe.Query("index 1")).indices.tolist()
    )


def test_file_reads_infer_bonds_and_leave_their_input_unchanged():
    # A file read applies the default perception pass, so the GLY backbone
    # carries its standard bond without the caller asking for it.
    structure = molframe.read(DATA / "basic.pdb")
    assert structure.bond_count == 1
    bonded = structure.infer_bonds()
    assert bonded.bond_count == 1
    assert structure.bond_count == 1
    with pytest.raises(ValueError):
        structure.infer_bonds(scale=0.0)


def test_a_failed_query_quotes_itself_and_names_the_fix():
    structure = molframe.read(DATA / "basic.pdb")
    with pytest.raises(molframe.QueryError) as raised:
        structure.select("name CA and bogus")
    message = str(raised.value)
    assert "MOLFRAME-E4004" in message
    assert "bogus" in message
    assert "^" in message
    assert "help:" in message
    assert isinstance(raised.value, ValueError)


def test_a_valid_but_suspicious_query_warns():
    structure = molframe.read(DATA / "basic.pdb")
    with pytest.warns(molframe.QueryWarning, match="W4001"):
        structure.select("name CA and name N or name C")


def test_atom_names_and_elements_match_as_written():
    structure = molframe.read(DATA / "basic.pdb")
    assert len(structure.select("name CA")) == 1
    assert len(structure.select("element c")) == len(structure.select("element C"))
