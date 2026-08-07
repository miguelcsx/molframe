#!/usr/bin/env python3
"""Compare pdbiox's unambiguous parse surface with gemmi and Biotite."""

from pathlib import Path
from sys import argv

import gemmi
import numpy as np
import pdbiox
from biotite.structure.io import load_structure


def pdbiox_values(path: Path):
    structure = pdbiox.read(path)
    return {
        "models": structure.model_count,
        "chains": structure.chain_count,
        "residues": structure.residue_count,
        "names": [structure.atoms[i].name for i in range(structure.atom_count)],
        "elements": [structure.atoms[i].element for i in range(structure.atom_count)],
        "xyz": np.asarray(structure.xyz),
    }


def gemmi_values(path: Path):
    structure = gemmi.read_structure(str(path))
    atoms = list(structure[0].all())
    return {
        "models": len(structure),
        "chains": len(structure[0]),
        "residues": sum(len(chain) for chain in structure[0]),
        "names": [item.atom.name for item in atoms],
        "elements": [item.atom.element.name for item in atoms],
        "xyz": np.array(
            [[item.atom.pos.x, item.atom.pos.y, item.atom.pos.z] for item in atoms]
        ),
    }


def biotite_values(path: Path):
    structure = load_structure(path, model=None)
    if structure.coord.ndim == 3:
        model_count = structure.coord.shape[0]
        structure = structure[0]
    else:
        model_count = 1
    return {
        "models": model_count,
        "chains": len(np.unique(structure.chain_id)),
        "residues": len(
            np.unique(np.stack((structure.chain_id, structure.res_id), axis=1), axis=0)
        ),
        "names": structure.atom_name.tolist(),
        "elements": structure.element.tolist(),
        "xyz": structure.coord,
    }


def compare(path: Path):
    implementations = {
        "pdbiox": pdbiox_values(path),
        "gemmi": gemmi_values(path),
        "biotite": biotite_values(path),
    }
    expected = implementations["pdbiox"]
    for name, actual in implementations.items():
        for field in ("models", "chains", "residues", "names", "elements"):
            if actual[field] != expected[field]:
                raise AssertionError(
                    f"{path}: {name} differs for {field}: "
                    f"{actual[field]!r} != {expected[field]!r}"
                )
        if not np.allclose(actual["xyz"], expected["xyz"], rtol=0.0, atol=1e-6):
            raise AssertionError(f"{path}: {name} differs for coordinates")
    print(f"PASS {path}")


if __name__ == "__main__":
    for argument in argv[1:]:
        compare(Path(argument))
