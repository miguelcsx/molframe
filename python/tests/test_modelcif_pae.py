"""Dense PAE: ModelCIF round trip, zero-copy views and the JSON importer."""

from pathlib import Path

import numpy as np
import pytest

import molframe

GRID = [
    [0.25, 1.5, 31.75],
    [2.5, 0.25, 4.25],
    [9.5, 3.75, 0.5],
]
NAMES = ("MET", "ALA", "GLY")


def _model_cif() -> str:
    lines = [
        "data_pae",
        "loop_",
        "_ma_qa_metric.id",
        "_ma_qa_metric.name",
        "_ma_qa_metric.type",
        "_ma_qa_metric.mode",
        "2 PAE PAE local-pairwise",
        "#",
        "loop_",
        "_ma_qa_metric_local_pairwise.ordinal_id",
        "_ma_qa_metric_local_pairwise.model_id",
        "_ma_qa_metric_local_pairwise.label_asym_id_1",
        "_ma_qa_metric_local_pairwise.label_seq_id_1",
        "_ma_qa_metric_local_pairwise.label_comp_id_1",
        "_ma_qa_metric_local_pairwise.label_asym_id_2",
        "_ma_qa_metric_local_pairwise.label_seq_id_2",
        "_ma_qa_metric_local_pairwise.label_comp_id_2",
        "_ma_qa_metric_local_pairwise.metric_id",
        "_ma_qa_metric_local_pairwise.metric_value",
    ]
    for index in range(9):
        a, b = divmod(index, 3)
        lines.append(f"{index + 1} 1 A {a + 1} {NAMES[a]} A {b + 1} {NAMES[b]} 2 {GRID[a][b]}")
    lines += [
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    ]
    for residue, name in enumerate(NAMES):
        lines.append(f"ATOM {residue + 1} C CA {name} A {residue + 1} {residue}.0 0.0 0.0")
    lines.append("#")
    return "\n".join(lines) + "\n"


@pytest.fixture
def structure(tmp_path: Path) -> molframe.Structure:
    path = tmp_path / "model.cif"
    path.write_text(_model_cif())
    return molframe.read(path)


def test_pae_is_stored_in_the_narrowest_exact_layout(structure: molframe.Structure) -> None:
    matrix = molframe.modelcif.pae(structure)
    assert matrix is not None
    assert matrix.size == 3
    assert matrix.dtype == "uint16"
    assert matrix.denominator == 100
    assert matrix.is_exact
    assert matrix.nbytes == 18
    values = matrix.array().astype(np.float64) / matrix.denominator
    assert np.array_equal(values, np.array(GRID))
    assert matrix.get(0, 2) == 31.75


def test_array_is_a_read_only_zero_copy_view(structure: molframe.Structure) -> None:
    matrix = molframe.modelcif.pae(structure)
    assert matrix is not None
    first = matrix.array()
    second = matrix.array()
    assert not first.flags.owndata
    assert not first.flags.writeable
    assert np.shares_memory(first, second)


def test_structure_without_pae_has_none(tmp_path: Path) -> None:
    path = tmp_path / "atoms.cif"
    text = _model_cif()
    path.write_text("data_x\n" + text[text.index("loop_\n_atom_site") :])
    assert molframe.modelcif.pae(molframe.read(path)) is None


def test_written_mmcif_reads_back_to_the_same_matrix(
    structure: molframe.Structure, tmp_path: Path
) -> None:
    text = molframe.formats.to_mmcif(structure)
    path = tmp_path / "again.cif"
    path.write_text(text)
    again = molframe.modelcif.pae(molframe.read(path))
    original = molframe.modelcif.pae(structure)
    assert again is not None
    assert original is not None
    assert np.array_equal(again.array(), original.array())
    assert again.denominator == original.denominator


def test_pae_json_importer_keeps_integers_in_one_byte(tmp_path: Path) -> None:
    path = tmp_path / "pae.json"
    path.write_text(
        '[{"predicted_aligned_error":[[0,1],[2,0]],"max_predicted_aligned_error":31.75}]'
    )
    matrix = molframe.modelcif.read_pae_json(path)
    assert matrix.dtype == "uint8"
    assert matrix.declared_maximum == 31.75
    assert np.array_equal(matrix.array(), np.array([[0, 1], [2, 0]]))


def test_forced_quantisation_reports_its_error(tmp_path: Path) -> None:
    path = tmp_path / "pae.json"
    path.write_text('[{"predicted_aligned_error":[[0.1,300.5],[0.3,0.4]]}]')
    matrix = molframe.modelcif.read_pae_json(path, quantise="uint8")
    assert matrix.dtype == "uint8"
    assert not matrix.is_exact
    assert matrix.max_quantisation_error > 0


def test_budget_and_bad_quantise_are_refused(tmp_path: Path) -> None:
    path = tmp_path / "pae.json"
    path.write_text('[{"predicted_aligned_error":[[0,1],[2,0]]}]')
    with pytest.raises(ValueError, match=r"budget|exceed|bytes"):
        molframe.modelcif.read_pae_json(path, max_bytes=1)
    with pytest.raises(ValueError, match="quantise"):
        molframe.modelcif.read_pae_json(path, quantise="int4")  # pyright: ignore[reportArgumentType]


PLDDT_FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "crates"
    / "molframe-py"
    / "tests"
    / "data"
    / "modelcif_plddt.cif"
)


def test_confidence_bands_agree_with_selection_and_carry_colours() -> None:
    structure = molframe.read(PLDDT_FIXTURE)
    bands = molframe.modelcif.plddt_bands(structure)
    assert bands is not None
    names = molframe.modelcif.plddt_band_names()
    assert names == ["very_low", "low", "confident", "very_high"]
    colours = molframe.modelcif.plddt_band_colors()
    assert colours.shape == (4, 3)
    assert colours.dtype == np.uint8
    scored = bands != molframe.modelcif.PLDDT_UNSCORED
    assert scored.any()
    assert not scored.all()
    for code, name in enumerate(names):
        selected = structure.select(f"plddt_{name}")
        assert len(selected) == int((bands == code).sum()), name
