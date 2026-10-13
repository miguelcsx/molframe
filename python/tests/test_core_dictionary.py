"""Component chemistry without a user-supplied dictionary."""

import warnings
from pathlib import Path

import pytest

import molframe

ROOT = Path(__file__).resolve().parents[2]
CCD = ROOT / "crates/molframe-chem/data/CCD-amino-acids.cif"
ROWS = (
    "data_s\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n"
    "_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n"
    "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
)


def _read(tmp_path: Path, rows: str) -> molframe.Structure:
    path = tmp_path / "s.cif"
    path.write_text(ROWS + rows)
    return molframe.read(path)


def test_annotate_defaults_to_the_bundled_core_set(tmp_path: Path) -> None:
    structure = _read(
        tmp_path, "1 O O HOH A . 0 0 0\n2 H H1 HOH A . 0.96 0 0\n3 H H2 HOH A . -0.24 0.93 0\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        annotated = molframe.chemistry.annotate(structure)
    assert annotated.bond_count == 2


def test_unknown_component_names_itself_and_how_to_supply_a_dictionary(tmp_path: Path) -> None:
    structure = _read(tmp_path, "1 C C1 ZZZ A . 0 0 0\n")
    with pytest.warns(UserWarning, match="ZZZ.*components"):
        molframe.chemistry.annotate(structure)


def test_user_dictionary_layers_over_the_core_set(tmp_path: Path) -> None:
    structure = _read(
        tmp_path, "1 O O HOH A . 0 0 0\n2 H H1 HOH A . 0.96 0 0\n3 H H2 HOH A . -0.24 0.93 0\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        molframe.chemistry.annotate(structure, CCD)


def test_in_file_definitions_resolve_without_a_dictionary(tmp_path: Path) -> None:
    path = tmp_path / "s.cif"
    path.write_text(
        ROWS
        + "1 C C1 ZZZ A . 0 0 0\n2 O O1 ZZZ A . 1.2 0 0\n"
        + "#\nloop_\n_chem_comp.id\n_chem_comp.type\n_chem_comp.name\n"
        + "ZZZ non-polymer 'TEST'\n#\n"
        + "loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n"
        + "_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n"
        + "_chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_leaving_atom_flag\n"
        + "_chem_comp_atom.pdbx_stereo_config\n"
        + "ZZZ C1 C 0 N N N\nZZZ O1 O 0 N N N\n#\n"
        + "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n"
        + "_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n"
        + "_chem_comp_bond.pdbx_aromatic_flag\n_chem_comp_bond.pdbx_stereo_config\n"
        + "ZZZ C1 O1 DOUB N N\n#\n"
    )
    structure = molframe.read(path)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        annotated = molframe.chemistry.annotate(structure)
    assert annotated.bond_count == 1
