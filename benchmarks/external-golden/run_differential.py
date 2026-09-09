#!/usr/bin/env python3
"""Run a registered external golden pair and compare its machine output."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DRIVER_ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_driver(python: Path, driver: Path, fixture: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(python), str(driver), str(fixture)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    if completed.stderr:
        raise RuntimeError(f"{driver.name} wrote to stderr: {completed.stderr.rstrip()}")
    value = json.loads(completed.stdout)
    if not isinstance(value, dict):
        raise ValueError(f"{driver.name} did not produce a JSON object")
    return value


def maximum_delta(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        raise ValueError(f"array lengths differ: {len(left)} != {len(right)}")
    return max((abs(float(a) - float(b)) for a, b in zip(left, right)), default=0.0)


def compare_gw_043(pdbiox: dict[str, Any], reference: dict[str, Any]) -> dict[str, Any]:
    for field in ("fixture_sha256", "units", "policy", "frame_count", "atom_count"):
        if pdbiox.get(field) != reference.get(field):
            raise ValueError(f"GW-043 policy field differs: {field}")
    tolerances = {"rmsd": 1.0e-4, "rmsf": 1.0e-4}
    deltas: dict[str, dict[str, float]] = {}
    for quantity, tolerance in tolerances.items():
        deltas[quantity] = {}
        for library in ("mdanalysis", "mdtraj"):
            delta = maximum_delta(pdbiox[quantity], reference[library][quantity])
            deltas[quantity][library] = delta
            if delta > tolerance:
                raise ValueError(
                    f"GW-043 {quantity} pdbiox/{library} delta {delta} > {tolerance}"
                )
    return {"tolerances": tolerances, "maximum_absolute_deltas": deltas}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow", choices=("GW-043",))
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    fixture = ROOT / "target/stress-data/sota/lattice-100k-100.xtc"
    pdbiox_driver = DRIVER_ROOT / "gw_043_pdbiox.py"
    reference_driver = DRIVER_ROOT / "gw_043_reference.py"
    try:
        pdbiox = run_driver(arguments.python, pdbiox_driver, fixture)
        reference = run_driver(arguments.python, reference_driver, fixture)
        comparison = compare_gw_043(pdbiox, reference)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError, json.JSONDecodeError) as error:
        print(f"FAIL  {arguments.workflow}: {error}", file=sys.stderr)
        return 1
    report = {
        "schema_version": 1,
        "workflow": arguments.workflow,
        "status": "passed",
        "fixture": {
            "path": str(fixture.relative_to(ROOT)),
            "sha256": sha256(fixture),
        },
        "drivers": {
            "pdbiox": {
                "path": str(pdbiox_driver.relative_to(ROOT)),
                "sha256": sha256(pdbiox_driver),
                "version": pdbiox["version"],
            },
            "reference": {
                "path": str(reference_driver.relative_to(ROOT)),
                "sha256": sha256(reference_driver),
                "versions": reference["versions"],
            },
        },
        "policy": pdbiox["policy"],
        "shape": [pdbiox["frame_count"], pdbiox["atom_count"], 3],
        "comparison": comparison,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS  {arguments.workflow}: {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
