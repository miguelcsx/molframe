"""A ragged ensemble keeps its atoms reachable instead of reporting itself empty.

A multi-model file whose models hold different atoms has no single atom axis, so
the outer structure carries the models rather than a flat topology. Its counts are
therefore zero, which a caller cannot tell apart from an empty molecule. These
tests pin the flag and the accessor that make the difference observable, and check
that the atoms really are there rather than merely advertised.
"""

from pathlib import Path

import molframe

# Model 1 holds two atoms, model 2 holds one. The models differ, which is what
# makes the reader materialise them as independent structures.
RAGGED = """MODEL        1
ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      1.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1      5.000   0.000   0.000  1.00 20.00           N
ENDMDL
END
"""

# The same atom set in both models: one topology, so the flat accessors apply.
DENSE = """MODEL        1
ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      1.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1      5.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      6.000   0.000   0.000  1.00 20.00           C
ENDMDL
END
"""


def written(tmp_path: Path, text: str) -> molframe.Structure:
    path = tmp_path / "ensemble.pdb"
    path.write_text(text)
    return molframe.read(path)


def test_a_ragged_ensemble_says_so_and_exposes_each_model(tmp_path: Path) -> None:
    structure = written(tmp_path, RAGGED)

    assert structure.model_count == 2
    assert structure.is_ragged_ensemble
    # The shell has no topology of its own: these counts are zero, which is why
    # the caller must not read them as "this file has no atoms".
    assert structure.atom_count == 0

    models = structure.ragged_models
    assert len(models) == 2
    assert [model.atom_count for model in models] == [2, 1]
    assert [model.chain_count for model in models] == [1, 1]
    # The atoms are really present, not merely counted. `elements` holds atomic
    # numbers: 7 is nitrogen and 6 is carbon.
    assert models[0].coordinates.shape == (2, 3)
    assert models[1].coordinates.shape == (1, 3)
    assert models[0].elements.tolist() == [7, 6]
    assert [atom.element for atom in models[0].atoms] == ["N", "C"]
    assert not any(model.is_ragged_ensemble for model in models)


def test_a_ragged_models_atoms_are_the_rows_the_file_holds(tmp_path: Path) -> None:
    structure = written(tmp_path, RAGGED)

    first, second = structure.ragged_models
    # Each model resolves its own rows, in deposition order.
    assert first.coordinates.tolist()[1] == [1.0, 0.0, 0.0]
    assert second.coordinates.tolist()[0] == [5.0, 0.0, 0.0]
    # And each answers a selection over its own atoms.
    assert len(first.select("all").indices) == 2
    assert len(second.select("all").indices) == 1


def test_a_dense_ensemble_is_not_ragged_and_needs_no_model_iteration(tmp_path: Path) -> None:
    structure = written(tmp_path, DENSE)

    assert structure.model_count == 2
    assert not structure.is_ragged_ensemble
    assert structure.ragged_models == []
    # One shared topology, so the flat accessors answer directly.
    assert structure.atom_count == 2
    assert structure.coordinates.shape == (2, 3)
