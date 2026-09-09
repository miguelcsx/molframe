#!/usr/bin/env python3
"""Gemmi, Biotite and Biopython side of the structure-interface workflow."""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import tempfile

import gemmi
import numpy as np
import biotite.structure as struc
from biotite.structure.io import pdbx
from Bio.PDB import MMCIFIO, MMCIFParser, ShrakeRupley, Superimposer


CHAINS = ("A", "B")
CUTOFF = 4.0
PROBE = 1.4
POINTS = 100
BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
         "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "MG": 1.73}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def key(chain: str, seq: int, insertion: str, component: str, atom: str, altloc: str) -> str:
    normalized_altloc = altloc.replace("\x00", "").strip()
    return f"{chain}:{seq}:{insertion.strip()}:{component}:{atom.strip()}:{normalized_altloc}"


def gemmi_atoms(model: object) -> tuple[list[str], np.ndarray, np.ndarray, list[bool]]:
    keys: list[str] = []
    coordinates: list[list[float]] = []
    radii: list[float] = []
    first: list[bool] = []
    for chain in model:
        if chain.name not in CHAINS:
            continue
        for residue in chain:
            for atom in residue:
                element = atom.element.name.upper()
                if element == "H":
                    continue
                radius = BONDI.get(element)
                if radius is None:
                    raise ValueError(f"missing Bondi radius for {element}")
                keys.append(key(chain.name, residue.seqid.num, residue.seqid.icode, residue.name,
                                atom.name, atom.altloc))
                coordinates.append([atom.pos.x, atom.pos.y, atom.pos.z])
                radii.append(radius)
                first.append(chain.name == CHAINS[0])
    if len(keys) != len(set(keys)):
        raise ValueError("Gemmi author atom identifiers are not unique")
    return keys, np.asarray(coordinates), np.asarray(radii), first


def biotite_atoms(array: object) -> tuple[list[str], np.ndarray, list[bool]]:
    selection = np.isin(array.chain_id, CHAINS) & (array.element != "H")
    selected = array[selection]
    keys = [
        key(chain, int(sequence), insertion, component, atom, "")
        for chain, sequence, insertion, component, atom in zip(
            selected.chain_id, selected.res_id, selected.ins_code,
            selected.res_name, selected.atom_name
        )
    ]
    if len(keys) != len(set(keys)):
        raise ValueError("Biotite author atom identifiers are not unique")
    return keys, np.asarray(selected.coord), (selected.chain_id == CHAINS[0]).tolist()


def biopython_atoms(model: object) -> tuple[list[str], np.ndarray, list[bool]]:
    keys: list[str] = []
    coordinates: list[np.ndarray] = []
    first: list[bool] = []
    for chain in model:
        if chain.id not in CHAINS:
            continue
        for residue in chain:
            _, sequence, insertion = residue.id
            for atom in residue:
                if atom.element.upper() == "H":
                    continue
                keys.append(key(chain.id, sequence, insertion, residue.resname,
                                atom.name, atom.get_altloc()))
                coordinates.append(atom.coord)
                first.append(chain.id == CHAINS[0])
    if len(keys) != len(set(keys)):
        raise ValueError("Biopython author atom identifiers are not unique")
    return keys, np.asarray(coordinates), first


def contact_summary(keys: list[str], coordinates: np.ndarray, first: list[bool]) -> dict[str, object]:
    membership = np.asarray(first)
    left = np.flatnonzero(membership)
    right = np.flatnonzero(~membership)
    pairs: list[str] = []
    distances: list[float] = []
    cutoff_squared = CUTOFF * CUTOFF
    for start in range(0, len(left), 256):
        local = left[start:start + 256]
        delta = coordinates[local, None, :] - coordinates[right, :]
        squared = np.einsum("ijk,ijk->ij", delta, delta)
        rows, columns = np.nonzero(squared <= cutoff_squared)
        for row, column in zip(rows.tolist(), columns.tolist()):
            pairs.append(f"{keys[int(local[row])]}|{keys[int(right[column])]}" )
            distances.append(float(np.sqrt(squared[row, column])))
    pairs.sort()
    residues = sorted({":".join(item.split(":", 4)[:4]) for pair in pairs for item in pair.split("|")})
    return {
        "count": len(pairs),
        "pair_sha256": hashlib.sha256("\n".join(pairs).encode()).hexdigest(),
        "interface_residues": residues,
        "minimum_distance": min(distances) if distances else None,
        "maximum_distance": max(distances) if distances else None,
    }


def biotite_sasa(array: object, radii: np.ndarray, first: list[bool]) -> dict[str, float]:
    membership = np.asarray(first)

    def total(selection: np.ndarray) -> float:
        values = struc.sasa(array[selection], probe_radius=PROBE, ignore_ions=False,
                            point_number=POINTS, point_distr="Fibonacci",
                            vdw_radii=radii[selection])
        return float(np.nansum(values))

    first_alone = total(membership)
    second_alone = total(~membership)
    together = total(np.ones(len(array), dtype=bool))
    return {"first_alone": first_alone, "second_alone": second_alone,
            "together": together, "buried": first_alone + second_alone - together}


def biopython_sasa(model: object) -> dict[str, float]:
    selected = copy.deepcopy(model)
    for chain in list(selected):
        if chain.id not in CHAINS:
            selected.detach_child(chain.id)
            continue
        for residue in list(chain):
            for atom in list(residue):
                if atom.element.upper() == "H":
                    residue.detach_child(atom.id)
    calculator = ShrakeRupley(probe_radius=PROBE, n_points=POINTS, radii_dict=BONDI)

    def area(chains: tuple[str, ...]) -> float:
        candidate = copy.deepcopy(selected)
        for chain in list(candidate):
            if chain.id not in chains:
                candidate.detach_child(chain.id)
        calculator.compute(candidate, level="A")
        return float(sum(atom.sasa for atom in candidate.get_atoms()))

    first_alone = area((CHAINS[0],))
    second_alone = area((CHAINS[1],))
    together = area(CHAINS)
    return {"first_alone": first_alone, "second_alone": second_alone,
            "together": together, "buried": first_alone + second_alone - together}


def matched_ca_gemmi(model: object) -> tuple[list[object], list[object]]:
    atoms: dict[str, dict[str, object]] = {chain: {} for chain in CHAINS}
    for chain in model:
        if chain.name not in CHAINS:
            continue
        for residue in chain:
            atom = residue.find_atom("CA", "*")
            if atom is not None:
                atoms[chain.name][f"{residue.seqid.num}:{residue.seqid.icode.strip()}"] = atom
    shared = sorted(set(atoms[CHAINS[0]]) & set(atoms[CHAINS[1]]))
    return ([atoms[CHAINS[0]][item] for item in shared],
            [atoms[CHAINS[1]][item] for item in shared])


def superpositions(gemmi_model: object, biotite_array: object, bio_model: object) -> dict[str, object]:
    fixed, mobile = matched_ca_gemmi(gemmi_model)
    gemmi_fit = gemmi.superpose_positions([atom.pos for atom in mobile],
                                          [atom.pos for atom in fixed])
    mask_a = (biotite_array.chain_id == CHAINS[0]) & (biotite_array.atom_name == "CA")
    mask_b = (biotite_array.chain_id == CHAINS[1]) & (biotite_array.atom_name == "CA")
    bio_fixed = [residue["CA"] for residue in bio_model[CHAINS[0]] if "CA" in residue]
    bio_mobile = [residue["CA"] for residue in bio_model[CHAINS[1]] if "CA" in residue]
    bio_fit = Superimposer()
    bio_fit.set_atoms(bio_fixed, bio_mobile)
    _, biotite_transform = struc.superimpose(biotite_array.coord[mask_a],
                                              biotite_array.coord[mask_b])
    fitted = biotite_transform.apply(biotite_array.coord[mask_b])
    biotite_rmsd = float(struc.rmsd(biotite_array.coord[mask_a], fitted))
    return {
        "gemmi": {"matched_atoms": len(fixed), "rmsd": float(gemmi_fit.rmsd)},
        "biotite": {"matched_atoms": int(np.count_nonzero(mask_a)), "rmsd": biotite_rmsd},
        "biopython": {"matched_atoms": len(bio_fixed), "rmsd": float(bio_fit.rms)},
    }


def exports(directory: Path, gemmi_structure: object, biotite_array: object,
            bio_structure: object) -> dict[str, object]:
    gemmi_path = directory / "gemmi.cif"
    gemmi_structure.make_mmcif_document().write_file(str(gemmi_path))
    gemmi_roundtrip = gemmi.read_structure(str(gemmi_path))[0]

    biotite_path = directory / "biotite.cif"
    biotite_file = pdbx.CIFFile()
    pdbx.set_structure(biotite_file, biotite_array)
    biotite_file.write(biotite_path)
    biotite_roundtrip = pdbx.get_structure(
        pdbx.CIFFile.read(biotite_path), model=1, use_author_fields=True, altloc="first"
    )

    biopython_path = directory / "biopython.cif"
    biopython_writer = MMCIFIO()
    biopython_writer.set_structure(bio_structure)
    biopython_writer.save(str(biopython_path))
    biopython_roundtrip = next(
        MMCIFParser(QUIET=True).get_structure("roundtrip", biopython_path).get_models()
    )
    return {
        "gemmi": {"roundtrip_atoms": sum(1 for _ in gemmi_roundtrip.all())},
        "biotite": {"roundtrip_atoms": len(biotite_roundtrip)},
        "biopython": {"roundtrip_atoms": sum(1 for _ in biopython_roundtrip.get_atoms())},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    arguments = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="pdbiox-gw042-reference-") as directory:
        source = Path(directory) / "input.cif"
        source.write_bytes(gzip.decompress(arguments.fixture.read_bytes()))
        gemmi_structure = gemmi.read_structure(str(source))
        gemmi_model = gemmi_structure[0]
        cif = pdbx.CIFFile.read(source)
        biotite_array = pdbx.get_structure(cif, model=1, use_author_fields=True,
                                           altloc="first", extra_fields=["occupancy"])
        bio_structure = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(
            "1aon", source
        )
        bio_model = next(bio_structure.get_models())
        keys, coordinates, radii, first = gemmi_atoms(gemmi_model)
        biotite_keys, biotite_coordinates, biotite_first = biotite_atoms(biotite_array)
        biopython_keys, biopython_coordinates, biopython_first = biopython_atoms(bio_model)
        selected_biotite = biotite_array[np.isin(biotite_array.chain_id, CHAINS)
                                         & (biotite_array.element != "H")]
        biotite_radii = np.asarray([BONDI[element.upper()] for element in selected_biotite.element])
        gemmi_assembly = gemmi.make_assembly(gemmi_structure.assemblies[0], gemmi_model,
                                             gemmi.HowToNameCopiedChain.AddNumber)
        biotite_assembly = pdbx.get_assembly(cif, assembly_id="1", model=1,
                                             use_author_fields=True, altloc="first")
        report = {
            "schema_version": 1,
            "workflow": "GW-042",
            "implementation": "Gemmi+Biotite+Biopython",
            "versions": {name: importlib.metadata.version(name) for name in
                         ("gemmi", "biotite", "biopython")},
            "fixture_sha256": sha256(arguments.fixture),
            "source": {
                "gemmi": {"atoms": sum(1 for _ in gemmi_model.all()), "author_chains": len(gemmi_model)},
                "biotite": {"atoms": len(biotite_array), "author_chains": len(set(biotite_array.chain_id))},
                "biopython": {"atoms": sum(1 for _ in bio_model.get_atoms()), "author_chains": len(bio_model)},
            },
            "assembly": {
                "gemmi": {"atoms": sum(1 for _ in gemmi_assembly.all()), "author_chains": len(gemmi_assembly)},
                "biotite": {"atoms": len(biotite_assembly), "author_chains": len(set(biotite_assembly.chain_id))},
                "biopython": {"comparable": False, "reason": "MMCIFParser does not construct biological assemblies"},
            },
            "selection": {"atoms": len(keys), "first_atoms": sum(first)},
            "contacts": {
                "gemmi": contact_summary(keys, coordinates, first),
                "biotite": contact_summary(
                    biotite_keys, biotite_coordinates, biotite_first
                ),
                "biopython": contact_summary(
                    biopython_keys, biopython_coordinates, biopython_first
                ),
            },
            "sasa": {
                "biotite": biotite_sasa(selected_biotite, biotite_radii, biotite_first),
                "biopython": biopython_sasa(bio_model),
                "gemmi": {"comparable": False, "reason": "Gemmi 0.7.5 exposes no SASA calculation"},
            },
            "superposition": superpositions(gemmi_model, biotite_array, bio_model),
            "export": exports(Path(directory), gemmi_structure, biotite_array, bio_structure),
        }
        print(json.dumps(report, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
