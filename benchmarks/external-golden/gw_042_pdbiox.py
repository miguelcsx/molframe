#!/usr/bin/env python3
"""pdbiox side of the matched-policy structure-interface workflow."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import tempfile

import numpy as np
import pdbiox


CHAINS = ("A", "B")
CUTOFF = 4.0
PROBE = 1.4
POINTS = 100


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atom_key(chain: object, residue: object, atom: object) -> str:
    chain_id = chain.auth_label or chain.label
    sequence = residue.auth_seq_id
    insertion = residue.ins_code or ""
    name = atom.auth_name or atom.name
    altloc = atom.altloc or ""
    return f"{chain_id}:{sequence}:{insertion}:{residue.auth_name or residue.name}:{name}:{altloc}"


def selected_atoms(structure: object) -> tuple[list[int], list[str], np.ndarray, np.ndarray, list[bool]]:
    indices: list[int] = []
    keys: list[str] = []
    coordinates: list[list[float]] = []
    radii: list[float] = []
    first: list[bool] = []
    for chain in structure.chains:
        chain_id = chain.auth_label or chain.label
        if chain_id not in CHAINS:
            continue
        for residue in chain.residues:
            for atom in residue.atoms:
                if atom.element == "H" or atom.coord is None:
                    continue
                element = pdbiox.Element.from_symbol(atom.element or "")
                if element is None:
                    raise ValueError(f"unknown element for {atom_key(chain, residue, atom)}")
                radius = element.vdw_radius(pdbiox.RadiusSet.Bondi)
                if radius is None:
                    raise ValueError(f"missing Bondi radius for {atom.element}")
                indices.append(atom.index)
                keys.append(atom_key(chain, residue, atom))
                coordinates.append(list(atom.coord))
                radii.append(radius)
                first.append(chain_id == CHAINS[0])
    if len(keys) != len(set(keys)):
        raise ValueError("selected author atom identifiers are not unique")
    return (
        indices,
        keys,
        np.asarray(coordinates, dtype=np.float32),
        np.asarray(radii, dtype=np.float32),
        first,
    )


def contacts(structure: object, indices: list[int], keys: list[str], first: list[bool]) -> dict[str, object]:
    left = pdbiox.Selection([index for index, member in zip(indices, first) if member])
    right = pdbiox.Selection([index for index, member in zip(indices, first) if not member])
    table = pdbiox.analysis.atom_contacts_between(
        structure, left, right, CUTOFF, pdbiox.SpatialBackend.CellList
    )
    key_by_index = dict(zip(indices, keys))
    pairs = sorted(
        f"{key_by_index[int(a)]}|{key_by_index[int(b)]}"
        for a, b in zip(table.first, table.second)
    )
    residues = sorted({":".join(key.split(":", 4)[:4]) for pair in pairs for key in pair.split("|")})
    return {
        "count": len(pairs),
        "pair_sha256": hashlib.sha256("\n".join(pairs).encode()).hexdigest(),
        "interface_residues": residues,
        "minimum_distance": float(np.min(table.distance)) if len(table) else None,
        "maximum_distance": float(np.max(table.distance)) if len(table) else None,
    }


def sasa(coordinates: np.ndarray, radii: np.ndarray, first: list[bool]) -> dict[str, float]:
    membership = np.asarray(first, dtype=np.bool_)
    areas = np.asarray(
        pdbiox.surface.shrake_rupley(coordinates, radii, PROBE, POINTS), dtype=np.float64
    )
    buried = pdbiox.surface.buried_surface(
        coordinates, radii, membership, PROBE, POINTS
    )
    return {
        "first_alone": buried.first_alone,
        "second_alone": buried.second_alone,
        "together": buried.together,
        "buried": buried.buried,
        "together_per_atom_sum": float(np.sum(areas)),
    }


def alpha_carbons(structure: object) -> tuple[np.ndarray, np.ndarray, list[str]]:
    by_chain: dict[str, dict[str, list[float]]] = {chain: {} for chain in CHAINS}
    for chain in structure.chains:
        chain_id = chain.auth_label or chain.label
        if chain_id not in CHAINS:
            continue
        for residue in chain.residues:
            try:
                atom = residue.atoms["CA"]
            except KeyError:
                continue
            if atom.coord is None:
                continue
            residue_key = f"{residue.auth_seq_id}:{residue.ins_code or ''}"
            by_chain[chain_id][residue_key] = list(atom.coord)
    shared = sorted(set(by_chain[CHAINS[0]]) & set(by_chain[CHAINS[1]]))
    fixed = np.asarray([by_chain[CHAINS[0]][key] for key in shared], dtype=np.float32)
    mobile = np.asarray([by_chain[CHAINS[1]][key] for key in shared], dtype=np.float32)
    return fixed, mobile, shared


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    arguments = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="pdbiox-gw042-") as directory:
        source = Path(directory) / "input.cif"
        source.write_bytes(gzip.decompress(arguments.fixture.read_bytes()))
        structure = pdbiox.read(source)
        indices, keys, coordinates, radii, first = selected_atoms(structure)
        assembly = structure.assembly("1", limit=100)
        fixed, mobile, matched = alpha_carbons(structure)
        superposition = pdbiox.superpose(mobile, fixed)
        export = Path(directory) / "roundtrip.cif"
        export.write_text(
            pdbiox.write_mmcif(
                structure, generate_connection_ids=True, connection_type="covale"
            )
        )
        roundtrip = pdbiox.read(export)
        report = {
            "schema_version": 1,
            "workflow": "GW-042",
            "implementation": "pdbiox",
            "version": importlib.metadata.version("pdbiox"),
            "fixture_sha256": sha256(arguments.fixture),
            "policy": {
                "model": 1,
                "identifier_namespace": "author",
                "chains": list(CHAINS),
                "atoms": "non-hydrogen atoms with coordinates",
                "contact": f"inclusive Euclidean distance <= {CUTOFF} angstrom",
                "sasa": f"Bondi radii, probe {PROBE} angstrom, {POINTS} Fibonacci points",
                "superposition": "unweighted rigid fit of author-residue-matched CA atoms",
            },
            "source": {
                "atoms": structure.atom_count,
                "models": structure.model_count,
                "author_chains": len({chain.auth_label for chain in structure.chains}),
            },
            "assembly": {
                "id": "1",
                "instances": assembly.instance_count,
                "chains": len(assembly.chains()),
                "atoms": len(assembly.atoms()),
            },
            "selection": {"atoms": len(keys), "first_atoms": sum(first)},
            "contacts": contacts(structure, indices, keys, first),
            "sasa": sasa(coordinates, radii, first),
            "superposition": {"matched_atoms": len(matched), "rmsd": superposition.rmsd},
            "export": {
                "format": "mmCIF",
                "roundtrip_atoms": roundtrip.atom_count,
                "roundtrip_models": roundtrip.model_count,
            },
        }
        print(json.dumps(report, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
