from pathlib import Path

import numpy as np
import pytest

import molframe

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "crates" / "molframe-py" / "tests" / "data"
FOUR_HHB_CIF = ROOT / "crates" / "molframe-bench" / "data" / "4hhb.cif"
FOUR_HHB_PDB = ROOT / "crates" / "molframe-bench" / "data" / "4hhb.pdb"
FOUR_HHB_BCIF = ROOT / "crates" / "molframe-bench" / "data" / "4hhb.bcif"

WATER_SDF = (
    "water\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0\n  1  3  1  0\nM  END\n$$$$\n"
)
WATER_MOL2 = (
    "@<TRIPOS>MOLECULE\nwat\n 3 2 1 0 0\nSMALL\nNO_CHARGES\n\n\n"
    "@<TRIPOS>ATOM\n"
    " 1 OW 0.0 0.0 0.0 O.3 1 WAT 0.0\n"
    " 2 HW1 0.757 0.586 0.0 H 1 WAT 0.0\n"
    " 3 HW2 -0.757 0.586 0.0 H 1 WAT 0.0\n"
    "@<TRIPOS>BOND\n 1 1 2 1\n 2 1 3 1\n"
)
SALT = (
    "data_salt\n_cell_length_a 5.64\n_cell_length_b 5.64\n_cell_length_c 5.64\n"
    "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
    "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n"
    "_atom_site_fract_y\n_atom_site_fract_z\nNa1 Na 0.0 0.0 0.0\nCl1 Cl 0.5 0.5 0.5\n"
)


@pytest.mark.parametrize(
    ("text", "name", "atoms", "bonds"),
    [(WATER_SDF, "sdf", 3, 2), (WATER_MOL2, "mol2", 3, 2), (SALT, "smallcif", 2, 0)],
)
def test_small_molecule_formats_read_through_format(text, name, atoms, bonds):
    structure = molframe.read(text.encode(), format=name)
    assert structure.atom_count == atoms
    assert structure.bond_count == bonds
    options = molframe.ReadOptions(format=name)
    assert molframe.read(text.encode(), options=options).atom_count == atoms
    assert molframe.Reader(text.encode()).read(format=name).atom_count == atoms


def test_a_suffix_is_a_format_name_and_unknown_names_are_refused():
    assert molframe.read(WATER_SDF.encode(), format="mol").atom_count == 3
    with pytest.raises(molframe.PolicyError, match=r"format: xyz .expected auto, mmcif") as raised:
        molframe.read(WATER_SDF.encode(), format="xyz")
    assert raised.value.code == "MOLFRAME-E6104"


def test_stating_the_format_twice_is_refused():
    with pytest.raises(molframe.MolframeValueError, match="stated twice"):
        molframe.read(WATER_SDF.encode(), format="sdf", options=molframe.ReadOptions(format="sdf"))


def test_read_options_carry_every_decision_by_name():
    options = molframe.ReadOptions(
        format="pdb",
        mode="strict",
        first_model_only=True,
        discard_hydrogens=True,
        missing_element="infer_from_atom_name",
        skip_categories=["REMARK"],
        max_nesting_depth=8,
    )
    assert options.format == "pdb"
    assert options.mode == "strict"
    assert options.first_model_only
    assert not options.coordinates_only
    assert options.discard_hydrogens
    assert options.missing_element == "infer_from_atom_name"
    assert options.ambiguous_residue_boundary == "reject"
    assert options.skip_categories == ["REMARK"]
    assert options.only_categories is None
    with pytest.raises(molframe.PolicyError, match=r"mode: lenient .expected strict"):
        molframe.ReadOptions(mode="lenient")
    with pytest.raises(molframe.MolframeValueError, match="not both"):
        molframe.ReadOptions(only_categories=["a"], skip_categories=["b"])


def test_a_category_filter_changes_what_is_read_and_never_the_atoms():
    full = molframe.read(FOUR_HHB_CIF)
    lean = molframe.read(
        FOUR_HHB_CIF,
        options=molframe.ReadOptions(skip_categories=["struct_conf", "struct_sheet_range"]),
    )
    assert lean.atom_count == full.atom_count
    assert np.array_equal(np.asarray(lean.coordinates), np.asarray(full.coordinates))
    helix = {str(each) for each in full.secondary_structure}
    assert len(helix) > 1


def test_model_and_hydrogen_options_apply():
    everything = molframe.read(FOUR_HHB_PDB)
    first = molframe.read(FOUR_HHB_PDB, options=molframe.ReadOptions(first_model_only=True))
    assert first.model_count == 1
    assert first.atom_count <= everything.atom_count
    coordinates_only = molframe.read(
        FOUR_HHB_PDB, options=molframe.ReadOptions(coordinates_only=True)
    )
    assert coordinates_only.bond_count < everything.bond_count, "perception is skipped"
    assert everything.bond_count > 0


def test_read_with_diagnostics_returns_the_findings_as_diagnostics():
    structure, findings = molframe.read_with_diagnostics(FOUR_HHB_PDB)
    assert structure.atom_count > 0
    assert isinstance(findings, list)
    assert all(isinstance(each, molframe.Diagnostic) for each in findings)
    assert all(each.code.startswith("MOLFRAME-") for each in findings)


def test_a_read_that_fails_carries_the_findings():
    with pytest.raises(molframe.MolframeError) as raised:
        molframe.read(b"this is not a structure", name="x.txt")
    assert raised.value.code == "MOLFRAME-E1001"
    assert raised.value.findings


def test_a_read_runs_under_a_context_and_obeys_its_cancellation():
    context = molframe.ExecutionContext(workers=1)
    assert molframe.read(FOUR_HHB_PDB, context=context).atom_count > 0
    cancelled = molframe.ExecutionContext()
    cancelled.cancel()
    with pytest.raises(molframe.Cancelled):
        molframe.read(FOUR_HHB_PDB, context=cancelled)
    starved = molframe.ExecutionContext(memory_budget=1024)
    with pytest.raises(molframe.MemoryBudgetError):
        molframe.read(FOUR_HHB_PDB, context=starved)


@pytest.mark.parametrize("path", [FOUR_HHB_CIF, FOUR_HHB_PDB, FOUR_HHB_BCIF])
def test_batches_cover_every_atom_in_file_order(path):
    whole = molframe.read(path, options=molframe.ReadOptions(first_model_only=True))
    batches = list(molframe.open_structure_batches(path, rows=1000))
    assert len(batches) > 1
    assert all(len(batch) <= 1000 for batch in batches)
    positions = np.concatenate([batch.positions for batch in batches])
    assert positions.shape == (whole.atom_count, 3)
    assert np.array_equal(positions, np.asarray(whole.coordinates))
    elements = np.concatenate([batch.elements for batch in batches])
    assert elements.shape == (whole.atom_count,)
    names = [name for batch in batches for name in batch.atom_names]
    assert names[:4] == [whole.atoms[index].name for index in range(4)]
    assert {chain for batch in batches for chain in batch.chains} >= {"A", "B"}


def test_a_batch_marks_unrecorded_values_with_nan_not_zero():
    batch = next(iter(molframe.open_structure_batches(FOUR_HHB_CIF)))
    assert batch.occupancies.shape == (len(batch),)
    assert not np.any(batch.occupancies == 0.0) or np.isnan(batch.occupancies).sum() >= 0
    assert batch.positions.shape == (len(batch), 3)
    assert isinstance(batch.diagnostics, list)


def test_batches_are_bounded_by_the_context_budget():
    starved = molframe.ExecutionContext(memory_budget=2048)
    with pytest.raises(molframe.MolframeError):
        list(molframe.open_structure_batches(FOUR_HHB_CIF, context=starved))


def test_a_format_without_a_bounded_reader_is_refused_not_read_whole(tmp_path):
    path = tmp_path / "water.sdf"
    path.write_text(WATER_SDF)
    with pytest.raises(molframe.MolframeError):
        list(molframe.open_structure_batches(path))


def test_write_states_its_decisions(tmp_path):
    structure = molframe.read(DATA / "basic.pdb")
    pdb = tmp_path / "renamed.pdb"
    molframe.formats.write(structure, pdb, chain_map={"A": "B"})
    assert molframe.read(pdb).chains[0].label == "B"
    odd = tmp_path / "x.dat"
    from_cif = molframe.read(DATA / "basic.cif")
    molframe.formats.write(from_cif, odd, format="mmcif")
    assert molframe.read(odd, format="mmcif").atom_count == from_cif.atom_count
    sdf = tmp_path / "water.sdf"
    molframe.formats.write(molframe.read(WATER_SDF.encode(), format="sdf"), sdf)
    assert molframe.read(sdf).atom_count == 3
    with pytest.raises(molframe.MolframeError):
        molframe.formats.write(structure, tmp_path / "x.cif", memory_limit=1)
    with pytest.raises(molframe.PolicyError):
        molframe.formats.write(structure, tmp_path / "y", format="nonsense")
