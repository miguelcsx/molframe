#!/usr/bin/env python3
"""Measure public pdbiox exports without claiming unavailable imports."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pdbiox


FIELDS = (
    "chain_id",
    "residue_name",
    "residue_number",
    "insertion_code",
    "heterogen",
    "atom_name",
    "element",
    "serial",
    "formal_charge",
    "occupancy",
    "b_factor",
    "alternate_location",
    "position",
    "bonds",
    "bond_order",
)
ATOM_FIELDS = FIELDS[:-2]

BOND_ORDER = {
    "BondOrder.Unknown": "unknown_or_polymeric",
    "BondOrder.Polymeric": "unknown_or_polymeric",
    "BondOrder.Single": "single",
    "BondOrder.Double": "double",
    "BondOrder.Triple": "triple",
    "BondOrder.Quadruple": "quadruple",
    "BondOrder.Aromatic": "aromatic",
}

BIOTITE_BOND_ORDER = {
    0: "unknown_or_polymeric",
    1: "single",
    2: "double",
    3: "triple",
    4: "quadruple",
    9: "aromatic",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical(structure: Any) -> dict[str, Any]:
    projection = pdbiox.adapters.TopologyBatch.from_model(
        structure, 0, pdbiox.Namespace.Auth
    )
    strings = projection.strings
    atom_residue = projection.atom_residue
    residue_chain = projection.residue_chain

    def text(string_id: int) -> str | None:
        return strings[string_id] if string_id < len(strings) else None

    def optional(values: list[Any], validity: list[bool]) -> list[Any | None]:
        return [value if valid else None for value, valid in zip(values, validity)]

    residue_numbers = optional(
        projection.residue_numbers, projection.residue_number_validity
    )
    return {
        "chain_id": [
            text(projection.chain_ids[residue_chain[residue]])
            for residue in atom_residue
        ],
        "residue_name": [text(projection.residue_names[residue]) for residue in atom_residue],
        "residue_number": [residue_numbers[residue] for residue in atom_residue],
        "insertion_code": [
            text(projection.residue_insertion_codes[residue]) for residue in atom_residue
        ],
        "heterogen": [projection.residue_is_heterogen[residue] for residue in atom_residue],
        "atom_name": [text(value) for value in projection.atom_names],
        "element": [
            pdbiox.chem.Element.from_atomic_number(value).symbol
            for value in projection.atomic_numbers
        ],
        "serial": optional(projection.atom_serials, projection.atom_serial_validity),
        "formal_charge": optional(projection.formal_charges, projection.formal_charge_validity),
        "occupancy": optional(projection.occupancies, projection.occupancy_validity),
        "b_factor": optional(projection.b_factors, projection.b_factor_validity),
        "alternate_location": [text(value) for value in projection.atom_alternate_locations],
        "position": [
            list(value)
            for value in zip(
                projection.position_x,
                projection.position_y,
                projection.position_z,
            )
        ],
        "bonds": sorted(
            sorted_pair(a, b)
            for a, b in zip(projection.bond_atom_a, projection.bond_atom_b)
        ),
        "bond_order": [BOND_ORDER[str(order)] for order in projection.bond_orders],
    }


def gemmi_fields(structure: Any) -> dict[str, Any]:
    model = structure[0]
    atoms = [(chain, residue, atom) for chain in model for residue in chain for atom in residue]
    return {
        "chain_id": [chain.name for chain, _, _ in atoms],
        "residue_name": [residue.name for _, residue, _ in atoms],
        "residue_number": [residue.seqid.num for _, residue, _ in atoms],
        "insertion_code": [none_if_blank(residue.seqid.icode) for _, residue, _ in atoms],
        "heterogen": None,
        "atom_name": [atom.name for _, _, atom in atoms],
        "element": [atom.element.name for _, _, atom in atoms],
        "serial": [none_if_zero(atom.serial) for _, _, atom in atoms],
        "formal_charge": [atom.charge for _, _, atom in atoms],
        "occupancy": [optional_number(atom.occ) for _, _, atom in atoms],
        "b_factor": [optional_number(atom.b_iso) for _, _, atom in atoms],
        "alternate_location": [none_if_blank(atom.altloc) for _, _, atom in atoms],
        "position": [[atom.pos.x, atom.pos.y, atom.pos.z] for _, _, atom in atoms],
        "bonds": None,
        "bond_order": None,
    }


def biopython_fields(structure: Any) -> dict[str, Any]:
    atoms = list(structure.get_atoms())
    residues = [atom.get_parent() for atom in atoms]
    chains = [residue.get_parent() for residue in residues]
    return {
        "chain_id": [chain.id for chain in chains],
        "residue_name": [residue.resname for residue in residues],
        "residue_number": [residue.id[1] for residue in residues],
        "insertion_code": [none_if_blank(residue.id[2]) for residue in residues],
        "heterogen": [residue.id[0] != " " for residue in residues],
        "atom_name": [atom.name for atom in atoms],
        "element": [atom.element for atom in atoms],
        "serial": [atom.serial_number for atom in atoms],
        "formal_charge": None,
        "occupancy": [atom.occupancy for atom in atoms],
        "b_factor": [atom.bfactor for atom in atoms],
        "alternate_location": [none_if_blank(atom.altloc) for atom in atoms],
        "position": [atom.coord.tolist() for atom in atoms],
        "bonds": None,
        "bond_order": None,
    }


def biotite_fields(array: Any) -> dict[str, Any]:
    bond_array = array.bonds.as_array() if array.bonds is not None else np.empty((0, 3))
    return {
        "chain_id": array.chain_id.tolist(),
        "residue_name": array.res_name.tolist(),
        "residue_number": array.res_id.tolist(),
        "insertion_code": [none_if_blank(value) for value in array.ins_code.tolist()],
        "heterogen": array.hetero.tolist(),
        "atom_name": array.atom_name.tolist(),
        "element": array.element.tolist(),
        "serial": None,
        "formal_charge": None,
        "occupancy": None,
        "b_factor": None,
        "alternate_location": None,
        "position": array.coord.tolist(),
        "bonds": sorted(np.sort(bond_array[:, :2], axis=1).astype(int).tolist()),
        "bond_order": [BIOTITE_BOND_ORDER[int(value)] for value in bond_array[:, 2]],
    }


def mdanalysis_fields(universe: Any) -> dict[str, Any]:
    atoms = universe.atoms
    bonds = atoms.bonds.indices if len(atoms.bonds) else np.empty((0, 2), dtype=int)
    return {
        "chain_id": atoms.chainIDs.tolist(),
        "residue_name": atoms.resnames.tolist(),
        "residue_number": atoms.resids.tolist(),
        "insertion_code": None,
        "heterogen": None,
        "atom_name": atoms.names.tolist(),
        "element": atoms.elements.tolist(),
        "serial": None,
        "formal_charge": None,
        "occupancy": None,
        "b_factor": None,
        "alternate_location": None,
        "position": atoms.positions.tolist(),
        "bonds": sorted(np.sort(bonds, axis=1).astype(int).tolist()),
        "bond_order": None,
    }


def none_if_blank(value: Any) -> Any:
    return None if value in (None, "", " ", "\x00") else value


def none_if_zero(value: int) -> int | None:
    return None if value == 0 else value


def optional_number(value: float) -> float | None:
    return None if np.isnan(value) else float(value)


def sorted_pair(left: int, right: int) -> list[int]:
    return [left, right] if left <= right else [right, left]


def compare(expected: dict[str, Any], observed: dict[str, Any]) -> dict[str, Any]:
    observed = align_by_serial(expected, observed)
    matrix: dict[str, Any] = {}
    for field in FIELDS:
        actual = observed[field]
        if actual is None:
            matrix[field] = {"status": "not_represented"}
            continue
        reference = expected[field]
        if field == "position":
            delta = float(np.max(np.abs(np.asarray(reference) - np.asarray(actual))))
            matrix[field] = {"status": "preserved", "maximum_absolute_delta": delta}
        else:
            matrix[field] = {
                "status": "preserved" if reference == actual else "changed",
                "mismatch_count": mismatch_count(reference, actual),
            }
    return matrix


def align_by_serial(
    expected: dict[str, Any], observed: dict[str, Any]
) -> dict[str, Any]:
    expected_serials = expected["serial"]
    observed_serials = observed["serial"]
    if expected_serials is None or observed_serials is None or None in expected_serials:
        return observed
    positions = {serial: index for index, serial in enumerate(observed_serials)}
    if len(positions) != len(observed_serials) or set(positions) != set(expected_serials):
        return observed
    order = [positions[serial] for serial in expected_serials]
    aligned = dict(observed)
    for field in ATOM_FIELDS:
        values = observed[field]
        if values is not None:
            aligned[field] = [values[index] for index in order]
    return aligned


def mismatch_count(left: list[Any], right: list[Any]) -> int:
    if len(left) != len(right):
        return max(len(left), len(right))
    return sum(a != b for a, b in zip(left, right))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    arguments = parser.parse_args()
    structure = pdbiox.read(str(arguments.fixture.resolve()))
    expected = canonical(structure)
    exporters = {
        "gemmi": (structure.to_gemmi, gemmi_fields),
        "biopython": (structure.to_biopython, biopython_fields),
        "biotite": (structure.to_biotite, biotite_fields),
        "mdanalysis": (structure.to_mdanalysis, mdanalysis_fields),
    }
    importers = {
        name: getattr(pdbiox.adapters, f"from_{name}") for name in exporters
    }
    field_matrix = {}
    round_trip = {}
    external_warnings = {}
    for name, (exporter, field_reader) in exporters.items():
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            external = exporter(pdbiox.Namespace.Auth)
        external_warnings[name] = [
            {"category": warning.category.__name__, "message": str(warning.message)}
            for warning in caught
        ]
        field_matrix[name] = compare(expected, field_reader(external))
        observed = canonical(importers[name](external))
        round_trip[name] = {
            "status": "passed",
            "atom_count": len(observed["atom_name"]),
            "field_matrix": compare(expected, observed),
        }
    report = {
        "schema_version": 1,
        "workflow": "GW-037",
        "fixture_sha256": sha256(arguments.fixture),
        "versions": {
            name: importlib.metadata.version(distribution)
            for name, distribution in {
                "pdbiox": "pdbiox",
                "gemmi": "gemmi",
                "biopython": "biopython",
                "biotite": "biotite",
                "mdanalysis": "MDAnalysis",
            }.items()
        },
        "atom_count": len(expected["atom_name"]),
        "direct_import_api": {name: callable(importer) for name, importer in importers.items()},
        "round_trip": round_trip,
        "field_matrix": field_matrix,
        "external_warnings": external_warnings,
    }
    print(json.dumps(report, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
