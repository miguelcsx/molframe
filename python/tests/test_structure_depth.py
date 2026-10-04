from pathlib import Path

import numpy as np
import pytest

import molframe

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "crates" / "molframe-py" / "tests" / "data"
BENCH = ROOT / "crates" / "molframe-bench" / "data"


def anisou_pdb() -> bytes:
    atom = "ATOM      1  CA  ALA A   1      11.104   6.134  -6.504  1.00 10.00           C  \n"
    other = "ATOM      2  CB  ALA A   1      12.560   6.195  -6.504  0.50 20.00           C  \n"
    anisou = "ANISOU    1  CA  ALA A   1     1000   2000   3000    100    200    300       C  \n"
    return (atom + anisou + other + "END\n").encode()


@pytest.fixture(scope="module")
def crambin():
    return molframe.read(BENCH / "1crn.cif")


def test_entry_metadata_states_only_what_the_file_says(crambin):
    metadata = crambin.metadata
    assert metadata.id == "1CRN"
    assert metadata.method == "X-RAY DIFFRACTION"
    assert metadata.resolution == pytest.approx(1.5)
    assert "1CRN" in repr(metadata)
    assert molframe.read(DATA / "basic.pdb").metadata.id is None
    assert molframe.read(BENCH / "2m7c.pdb").metadata.resolution is None


def test_a_cell_is_a_unit_cell_and_a_missing_one_is_none(crambin):
    cell = crambin.cell
    assert cell is not None
    assert cell.lengths == pytest.approx((40.96, 18.65, 22.52))
    assert cell.angles == pytest.approx((90.0, 90.77, 90.0))
    fractional = cell.to_fractional(cell.to_cartesian([0.25, 0.5, 0.75]))
    assert fractional == pytest.approx([0.25, 0.5, 0.75])
    assert cell == molframe.crystal.UnitCell(cell.lengths, cell.angles)
    assert molframe.read(BENCH / "2m7c.pdb").cell is None


def test_bulk_columns_agree_with_the_atom_handles(crambin):
    assert crambin.elements.shape == (crambin.atom_count,)
    assert crambin.occupancies.shape == (crambin.atom_count,)
    assert crambin.b_factors.shape == (crambin.atom_count,)
    for index in (0, 100, crambin.atom_count - 1):
        atom = crambin.atoms[index]
        assert crambin.elements[index] == atom.atomic_number
        assert crambin.occupancies[index] == pytest.approx(atom.occupancy)
        assert crambin.b_factors[index] == pytest.approx(atom.b_factor)
    assert np.array_equal(crambin.elements[:2], [7, 6])


def test_atoms_residues_and_chains_carry_their_attributes(crambin):
    atom = crambin.atoms[1]
    assert (atom.name, atom.element, atom.atomic_number) == ("CA", "C", 6)
    assert atom.serial == 2
    assert atom.altloc is None
    assert atom.anisotropy is None
    residue = atom.residue
    assert residue is not None
    assert residue.number == 1
    assert residue.auth_number == 1
    assert residue.name == "THR"
    assert residue.insertion_code is None
    assert not residue.is_hetero
    chain = residue.chain
    assert chain is not None
    assert chain.label == "A"
    assert chain.auth_label == "A"
    assert chain.polymer_kind == "protein"
    assert chain.entity == 0
    assert residue.secondary_structure == molframe.SecondaryStructure.Strand
    assert residue.secondary_source == molframe.SecondarySource.File


def test_secondary_source_names_where_each_state_came_from(crambin):
    sources = crambin.secondary_source
    assert len(sources) == crambin.residue_count
    assert set(sources) <= {
        molframe.SecondarySource.File,
        molframe.SecondarySource.Dssp,
        molframe.SecondarySource.CaOnly,
        molframe.SecondarySource.Unassigned,
    }
    assert molframe.SecondarySource.File in sources


def test_anisotropy_is_none_when_unstated_and_a_table_when_present(crambin):
    assert crambin.anisotropy is None
    structure = molframe.read(anisou_pdb(), name="a.pdb")
    table = structure.anisotropy
    assert table is not None
    assert table.names == ["atom", "u11", "u22", "u33", "u12", "u13", "u23"]
    assert list(table["atom"]) == [0]
    assert [float(table[name][0]) for name in table.names[1:]] == pytest.approx(
        [0.1, 0.2, 0.3, 0.01, 0.02, 0.03]
    )
    assert structure.atoms[0].anisotropy == pytest.approx((0.1, 0.2, 0.3, 0.01, 0.02, 0.03))
    assert structure.atoms[1].anisotropy is None
    assert structure.occupancies[1] == pytest.approx(0.5)


def test_entities_are_apart_from_the_chains_that_copy_them():
    structure = molframe.read(BENCH / "4hhb.cif")
    entities = structure.entities
    assert len(entities) == 5
    polymers = [
        entity for entity in (entities[i] for i in range(len(entities))) if entity.kind == "polymer"
    ]
    assert len(polymers) == 2
    alpha = polymers[0]
    assert len(alpha.sequence) == 141
    assert alpha.sequence[0] == "VAL"
    assert len(alpha.chains) == 2
    for chain_index in alpha.chains:
        assert structure.chains[chain_index].entity == alpha.index
    assert entities[-1].index == len(entities) - 1
    with pytest.raises(IndexError):
        entities[len(entities)]


def test_annotations_carry_values_and_which_atoms_record_them():
    structure = molframe.read(DATA / "modelcif_plddt.cif")
    annotations = structure.annotations
    assert "plddt" in annotations
    assert "plddt" in annotations.names
    assert len(annotations) == len(annotations.names)
    column = annotations["plddt"]
    assert column.kind == "real"
    assert column.values.shape == (structure.atom_count,)
    assert column.present.dtype == np.bool_
    assert column.present.shape == column.values.shape
    assert column.present.any()
    assert not column.present.all(), "some atoms record no confidence, and say so"
    recorded = column.values[column.present]
    assert np.all((recorded >= 0.0) & (recorded <= 100.0))
    kinds = annotations["component_kind"]
    assert kinds.kind in {"integer", "symbol", "real", "boolean"}
    with pytest.raises(KeyError):
        annotations["nope"]


def test_a_topology_edit_publishes_a_new_structure_and_leaves_the_original(crambin):
    plain = molframe.read(DATA / "basic.pdb")
    editor = plain.edit()
    editor.rename_chain(0, "Z")
    renamed = editor.finish()
    assert renamed.chains[0].label == "Z"
    assert plain.chains[0].label == "A"
    with pytest.raises(molframe.MolframeError, match="already been finished"):
        editor.rename_chain(0, "Y")
    with pytest.raises(molframe.ConsistencyError, match="extensions"):
        crambin.edit().rename_chain(0, "Z")


def test_deleting_a_selection_removes_exactly_those_atoms(crambin):
    selection = crambin.select("element S")
    assert 0 < len(selection) < crambin.atom_count
    refused = crambin.edit()
    with pytest.raises(molframe.ConsistencyError) as raised:
        refused.delete(selection)
    assert raised.value.code == "MOLFRAME-E3014"
    editor = crambin.edit()
    editor.clear_extensions()
    editor.delete(selection)
    smaller = editor.finish()
    assert smaller.atom_count == crambin.atom_count - len(selection)
    with pytest.warns(molframe.QueryWarning, match="not present"):
        assert len(smaller.select("element S")) == 0
    assert crambin.atom_count == len(crambin.coordinates)


def test_a_coordinate_edit_is_private_until_committed_and_published_as_a_copy(crambin):
    before = np.asarray(crambin.coordinates).copy()
    scoped = crambin.edit().coordinates()
    positions = scoped.positions()
    assert positions.shape == (crambin.atom_count, 3)
    assert positions.flags.writeable
    positions += np.float32(1.0)
    assert np.array_equal(np.asarray(crambin.coordinates), before), "the original is untouched"
    moved = scoped.commit()
    assert np.allclose(np.asarray(moved.coordinates), before + 1.0)
    positions += np.float32(1.0)
    assert np.allclose(np.asarray(moved.coordinates), before + 1.0), "a published copy is final"
    assert np.allclose(np.asarray(scoped.commit().coordinates), before + 2.0)


def test_a_coordinate_edit_is_charged_to_the_context_and_refuses_a_missing_model(crambin):
    with pytest.raises(molframe.MemoryBudgetError):
        crambin.edit().coordinates(context=molframe.ExecutionContext(memory_budget=16))
    scoped = crambin.edit().coordinates()
    with pytest.raises(molframe.MolframeError) as raised:
        scoped.positions(model=5)
    assert raised.value.code is not None
