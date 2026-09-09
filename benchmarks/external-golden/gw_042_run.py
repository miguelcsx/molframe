#!/usr/bin/env python3
"""Run GW-042 in one pinned environment and persist evidence only on success."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(driver: Path, fixture: Path) -> dict[str, object]:
    output = subprocess.run(
        [sys.executable, str(driver), str(fixture)], check=True, capture_output=True, text=True
    )
    return json.loads(output.stdout)


def relative_delta(first: float, second: float) -> float:
    return abs(first - second) / max(abs(first), abs(second), 1.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    arguments = parser.parse_args()
    directory = Path(__file__).resolve().parent
    pdbiox_driver = directory / "gw_042_pdbiox.py"
    reference_driver = directory / "gw_042_reference.py"
    candidate = run(pdbiox_driver, arguments.fixture)
    reference = run(reference_driver, arguments.fixture)
    failures: list[str] = []
    for field in ("atoms", "author_chains"):
        expected = candidate["source"][field]
        for library, values in reference["source"].items():
            if values[field] != expected:
                failures.append(f"source {field}: pdbiox={expected}, {library}={values[field]}")
    for library in ("gemmi", "biotite"):
        if reference["assembly"][library]["atoms"] != candidate["assembly"]["atoms"]:
            failures.append(f"assembly atoms differ for {library}")
    for library, values in reference["contacts"].items():
        for field in ("count", "pair_sha256", "interface_residues"):
            if candidate["contacts"][field] != values[field]:
                failures.append(f"contact {field} differs for {library}")
    rmsd_deltas: dict[str, float] = {}
    for library, values in reference["superposition"].items():
        if values["matched_atoms"] != candidate["superposition"]["matched_atoms"]:
            failures.append(f"superposition atom count differs for {library}")
        delta = abs(values["rmsd"] - candidate["superposition"]["rmsd"])
        rmsd_deltas[library] = delta
        if delta > 1.0e-3:
            failures.append(f"superposition RMSD delta {library}={delta}")
    sasa_deltas: dict[str, dict[str, float]] = {}
    for library, tolerance in (("biotite", 0.02), ("biopython", 0.05)):
        sasa_deltas[library] = {}
        for field in ("first_alone", "second_alone", "together", "buried"):
            delta = relative_delta(candidate["sasa"][field], reference["sasa"][library][field])
            sasa_deltas[library][field] = delta
            if delta > tolerance:
                failures.append(f"SASA relative delta {library}.{field}={delta}")
    if candidate["export"]["roundtrip_atoms"] != candidate["source"]["atoms"]:
        failures.append("pdbiox mmCIF export changed atom count")
    for library, values in reference["export"].items():
        if values["roundtrip_atoms"] != candidate["source"]["atoms"]:
            failures.append(f"{library} mmCIF export changed atom count")
    if failures:
        raise ValueError("; ".join(failures))
    result = {
        "schema_version": 1,
        "workflow": "GW-042",
        "status": "passed",
        "fixture": {"path": str(arguments.fixture), "sha256": candidate["fixture_sha256"]},
        "drivers": {
            "pdbiox": {"path": str(pdbiox_driver.relative_to(directory.parent.parent)),
                       "sha256": digest(pdbiox_driver), "version": candidate["version"]},
            "reference": {"path": str(reference_driver.relative_to(directory.parent.parent)),
                          "sha256": digest(reference_driver), "versions": reference["versions"]},
        },
        "policy": candidate["policy"],
        "shape": candidate["source"],
        "assembly": candidate["assembly"],
        "selection": candidate["selection"],
        "contacts": candidate["contacts"],
        "sasa": {"pdbiox": candidate["sasa"], "maximum_relative_deltas": sasa_deltas,
                 "non_comparable": {"gemmi": reference["sasa"]["gemmi"]["reason"]}},
        "superposition": {"pdbiox": candidate["superposition"],
                          "maximum_absolute_rmsd_deltas": rmsd_deltas},
        "export": {"pdbiox": candidate["export"], "references": reference["export"]},
        "contract_notes": {
            "contacts": "Each library parsed the same atom identity contract; each parse independently passes the same inclusive Euclidean predicate.",
            "assembly_chain_namespace": "Assembly atom cardinality is comparable; pdbiox reports label-asym instances while the references report author-chain cardinality.",
            "assembly_biopython": reference["assembly"]["biopython"]["reason"],
            "sasa_gemmi": reference["sasa"]["gemmi"]["reason"],
        },
    }
    target = directory / "gw_042_result.json"
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
