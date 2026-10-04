"""Match ligand SASA and burial with FreeSASA's Lee-Richards implementation.

Both programs use the same selected atoms, supplied Bondi radii, 20 slices and
float32-representable probe. No chemical bond inference is involved. Per-atom
areas and the fixed <1 square angstrom ligand-burial conclusion are retained.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from matched_interactions import AMINO_ACIDS, entity_kind, rdkit_transfer
from rdkit.Chem import rdFreeSASA

import molframe
from molframe import chemistry, surface


def compare(path: Path) -> list[dict]:
    structure = molframe.read(path)
    selected = list(map(int, structure.select("not element H").indices))
    protein = [
        i
        for i in selected
        if structure.atoms[i].residue.name in AMINO_ACIDS
        and entity_kind(structure, structure.atoms[i].residue) == "polymer"
    ]
    ligands = {}
    for i in selected:
        residue = structure.atoms[i].residue
        if entity_kind(structure, residue) == "non_polymer":
            ligands.setdefault(residue.index, []).append(i)
    radii = np.asarray(chemistry.vdw_radii(structure, radii="bondi"))
    probe = float(np.float32(1.4))
    options = rdFreeSASA.SASAOpts()
    options.algorithm = rdFreeSASA.SASAAlgorithm.LeeRichards
    options.probeRadius = probe
    rows = []
    for residue, atoms in sorted(ligands.items()):
        if len(atoms) < 5 or not protein:
            continue
        indices = protein + atoms
        values = np.ascontiguousarray(radii[indices], dtype=np.float32)
        if not np.isfinite(values).all():
            rows.append({"ligand_residue": residue, "refused": "missing Bondi radius"})
            continue
        xyz = np.ascontiguousarray(np.asarray(structure.coordinates)[indices])
        native = np.asarray(surface.lee_richards(xyz, values, probe=probe, slices=20))
        molecule = rdkit_transfer(structure, indices)
        total = rdFreeSASA.CalcSASA(molecule, list(map(float, values)), opts=options)
        reference = np.array([atom.GetDoubleProp("SASA") for atom in molecule.GetAtoms()])
        native_ligand, reference_ligand = native[len(protein) :], reference[len(protein) :]
        rows.append(
            {
                "ligand_residue": residue,
                "atoms": len(indices),
                "within_roundoff_tolerance": bool(np.max(np.abs(native - reference)) <= 1e-8),
                "max_atom_error_A2": float(np.max(np.abs(native - reference))),
                "total_error_A2": float(abs(native.sum() - total)),
                "native_ligand_area_A2": float(native_ligand.sum()),
                "reference_ligand_area_A2": float(reference_ligand.sum()),
                "native_buried_fraction": float(np.mean(native_ligand < 1.0)),
                "reference_buried_fraction": float(np.mean(reference_ligand < 1.0)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--ids", nargs="+", default=["1hvr", "1hsg", "3ptb", "1a4w"])
    arguments = parser.parse_args()
    manifest = json.loads(arguments.corpus.read_text())
    rows = []
    for record in manifest["files"]:
        path = arguments.corpus.parent / record["path"]
        if path.suffix != ".cif" or path.stem not in arguments.ids:
            continue
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            msg = f"frozen input changed: {path}"
            raise ValueError(msg)
        rows.extend({"id": path.stem, **row} for row in compare(path))
    print(
        json.dumps(
            {
                "algorithm": "Lee-Richards, 20 slices",
                "absolute_tolerance_A2": 1e-8,
                "molframe": molframe.__version__,
                "results": rows,
            },
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
