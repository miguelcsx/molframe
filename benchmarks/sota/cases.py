#!/usr/bin/env python3
"""Semantically checked competitor workloads for process-isolated benchmarks."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import pathlib
from typing import Any, Callable


Result = dict[str, int | float | str | None]
Workload = Callable[[pathlib.Path, pathlib.Path | None], Result]


def package_version(distribution: str) -> str:
    return importlib.metadata.version(distribution)


def gemmi_structure(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    import gemmi

    structure = gemmi.read_structure(str(path))
    atom_rows = 0
    chains = 0
    residues = 0
    for model in structure:
        for chain in model:
            chains += 1
            for residue in chain:
                residues += 1
                atom_rows += len(residue)
    return {
        "result_digest": atom_rows,
        "model_count": len(structure),
        "chain_count": chains,
        "residue_count": residues,
        "version": package_version("gemmi"),
    }


def biopython_structure(
    path: pathlib.Path,
    parser: Any,
    distribution: str = "biopython",
) -> Result:
    import warnings

    from Bio import BiopythonWarning

    warnings.simplefilter("ignore", BiopythonWarning)
    structure = parser.get_structure("benchmark", str(path))
    atom_rows = 0
    models = 0
    chains = 0
    residues = 0
    for model in structure:
        models += 1
        for chain in model:
            chains += 1
            for residue in chain:
                residues += 1
                atom_rows += len(residue.get_unpacked_list())
    return {
        "result_digest": atom_rows,
        "model_count": models,
        "chain_count": chains,
        "residue_count": residues,
        "version": package_version(distribution),
    }


def biopython_mmcif(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    from Bio.PDB import MMCIFParser

    return biopython_structure(path, MMCIFParser(QUIET=True))


def biopython_mmcif_fast(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    from Bio.PDB import FastMMCIFParser

    return biopython_structure(path, FastMMCIFParser(QUIET=True))


def biopython_pdb(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    from Bio.PDB import PDBParser

    return biopython_structure(path, PDBParser(QUIET=True))


def biotite_array(array: Any, distribution: str = "biotite") -> Result:
    from biotite.structure import AtomArrayStack

    if isinstance(array, AtomArrayStack):
        model_count = int(array.stack_depth())
        atom_rows = int(array.stack_depth() * array.array_length())
    else:
        model_count = 1
        atom_rows = int(array.array_length())
    return {
        "result_digest": atom_rows,
        "model_count": model_count,
        "chain_count": None,
        "residue_count": None,
        "version": package_version(distribution),
    }


def biotite_cif(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    import biotite.structure.io.pdbx as pdbx

    document = pdbx.CIFFile.read(path)
    return biotite_array(pdbx.get_structure(document, altloc="all"))


def biotite_bcif(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    import biotite.structure.io.pdbx as pdbx

    document = pdbx.BinaryCIFFile.read(path)
    return biotite_array(pdbx.get_structure(document, altloc="all"))


def biotite_pdb(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    import biotite.structure.io.pdb as pdb

    document = pdb.PDBFile.read(path)
    return biotite_array(document.get_structure(altloc="all"))


def mmtf_python(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    import mmtf

    decoded = mmtf.parse(str(path))
    return {
        "result_digest": int(decoded.num_atoms),
        "model_count": int(decoded.num_models),
        "chain_count": int(decoded.num_chains),
        "residue_count": int(decoded.num_groups),
        "version": package_version("mmtf-python"),
    }


def mrcfile_block(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    import mrcfile
    import numpy as np

    shape = (96, 96, 96)
    output = np.empty(shape, dtype=np.float32)
    with mrcfile.mmap(path, mode="r", permissive=False) as density:
        np.copyto(output, density.data[509:605, 271:367, 113:209])
        pointer = output.__array_interface__["data"][0]
        np.copyto(output, density.data[17:113, 401:497, 701:797])
        if output.__array_interface__["data"][0] != pointer:
            raise RuntimeError("mrcfile output buffer changed between equal blocks")
    return {
        "result_digest": int(output.size),
        "model_count": None,
        "checksum": float(output.sum(dtype=np.float64)),
        "version": package_version("mrcfile"),
    }


def mdanalysis_xtc(path: pathlib.Path, topology: pathlib.Path | None) -> Result:
    import MDAnalysis as mda

    if topology is None:
        raise ValueError("MDAnalysis XTC requires --topology")
    universe = mda.Universe(str(topology), str(path))
    frames = 0
    checksum = 0.0
    for timestep in universe.trajectory:
        frames += 1
        if timestep.n_atoms:
            checksum += float(timestep.positions[0].sum())
    return {
        "result_digest": frames * 1_000_000_000 + int(universe.atoms.n_atoms),
        "model_count": frames,
        "atom_count": int(universe.atoms.n_atoms),
        "checksum": checksum,
        "version": package_version("MDAnalysis"),
    }


def mdtraj_xtc(path: pathlib.Path, topology: pathlib.Path | None) -> Result:
    import mdtraj

    if topology is None:
        raise ValueError("MDTraj XTC requires --topology")
    frames = 0
    atoms = 0
    checksum = 0.0
    for chunk in mdtraj.iterload(str(path), top=str(topology), chunk=2):
        frames += int(chunk.n_frames)
        atoms = int(chunk.n_atoms)
        if chunk.n_frames and chunk.n_atoms:
            checksum += 10.0 * float(chunk.xyz[:, 0, :].sum())
    return {
        "result_digest": frames * 1_000_000_000 + atoms,
        "model_count": frames,
        "atom_count": atoms,
        "checksum": checksum,
        "version": package_version("mdtraj"),
    }


def biotite_xtc(path: pathlib.Path, _: pathlib.Path | None) -> Result:
    from biotite.structure.io.xtc import XTCFile

    coordinates = XTCFile.read(path).get_coord()
    frames, atoms, _ = coordinates.shape
    return {
        "result_digest": int(frames) * 1_000_000_000 + int(atoms),
        "model_count": int(frames),
        "atom_count": int(atoms),
        "checksum": float(coordinates[:, 0, :].sum(dtype="float64")),
        "version": package_version("biotite"),
    }


WORKLOADS: dict[str, Workload] = {
    "gemmi_structure": gemmi_structure,
    "biopython_mmcif": biopython_mmcif,
    "biopython_mmcif_fast": biopython_mmcif_fast,
    "biopython_pdb": biopython_pdb,
    "biotite_cif": biotite_cif,
    "biotite_bcif": biotite_bcif,
    "biotite_pdb": biotite_pdb,
    "mmtf_python": mmtf_python,
    "mrcfile_block": mrcfile_block,
    "mdanalysis_xtc": mdanalysis_xtc,
    "mdtraj_xtc": mdtraj_xtc,
    "biotite_xtc": biotite_xtc,
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("workload", choices=sorted(WORKLOADS))
    parser.add_argument("path", type=pathlib.Path)
    parser.add_argument("--topology", type=pathlib.Path)
    arguments = parser.parse_args()
    result = WORKLOADS[arguments.workload](arguments.path, arguments.topology)
    result.update({"schema_version": 1, "workload": arguments.workload})
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
