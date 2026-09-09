#!/usr/bin/env python3
"""Run GW-037 and persist evidence only after all contract checks pass."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DRIVER = Path(__file__).with_name("gw_037_pdbiox.py")
EXPECTED_FIXTURE_SHA256 = "b44610fb0a3d2c003dff40ca0d251e33232146de950439df5bb659a6d58ca3f8"
FIXTURE_SOURCE = "https://files.rcsb.org/download/1CJB.cif"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(report: dict[str, Any]) -> None:
    if report.get("workflow") != "GW-037":
        raise ValueError("driver reported the wrong workflow")
    if report.get("fixture_sha256") != EXPECTED_FIXTURE_SHA256:
        raise ValueError("fixture hash differs from the pinned 1CJB corpus")
    imports = report.get("direct_import_api")
    if imports != {name: True for name in ("gemmi", "biopython", "biotite", "mdanalysis")}:
        raise ValueError("direct-import availability changed; update the round-trip policy")
    round_trips = report.get("round_trip")
    if not isinstance(round_trips, dict) or set(round_trips) != set(imports):
        raise ValueError("round-trip evidence is incomplete")
    for library, result in round_trips.items():
        if result.get("status") != "passed":
            raise ValueError(f"{library} round-trip did not pass")
        matrix = result.get("field_matrix")
        if not isinstance(matrix, dict) or "position" not in matrix:
            raise ValueError(f"{library} round-trip field matrix is incomplete")
        position = matrix["position"]
        if position.get("status") != "preserved":
            raise ValueError(f"{library} round-trip did not preserve coordinates")
        if float(position.get("maximum_absolute_delta", float("inf"))) > 1.0e-5:
            raise ValueError(f"{library} round-trip coordinate delta exceeds 1e-5 angstrom")
    matrices = report.get("field_matrix")
    if not isinstance(matrices, dict) or set(matrices) != set(imports):
        raise ValueError("field-loss matrix is incomplete")
    for library, matrix in matrices.items():
        if not isinstance(matrix, dict) or "position" not in matrix:
            raise ValueError(f"{library} field-loss matrix is incomplete")
        position = matrix["position"]
        if position.get("status") != "preserved":
            raise ValueError(f"{library} did not preserve coordinates")
        if float(position.get("maximum_absolute_delta", float("inf"))) > 1.0e-5:
            raise ValueError(f"{library} coordinate delta exceeds 1e-5 angstrom")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        completed = subprocess.run(
            [str(arguments.python), str(DRIVER), str(arguments.fixture)],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        if completed.stderr:
            raise RuntimeError(f"driver wrote to stderr: {completed.stderr.rstrip()}")
        report = json.loads(completed.stdout)
        validate(report)
    except subprocess.CalledProcessError as error:
        detail = error.stderr.rstrip() if error.stderr else str(error)
        print(f"FAIL  GW-037: {detail}", file=sys.stderr)
        return 1
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        print(f"FAIL  GW-037: {error}", file=sys.stderr)
        return 1
    result = {
        "schema_version": 1,
        "workflow": "GW-037",
        "status": "passed",
        "fixture": {
            "id": "mmcif_1cjb",
            "source": FIXTURE_SOURCE,
            "sha256": sha256(arguments.fixture),
        },
        "driver": {
            "path": str(DRIVER.relative_to(ROOT)),
            "sha256": sha256(DRIVER),
        },
        **{key: value for key, value in report.items() if key not in {"schema_version", "workflow", "fixture_sha256"}},
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"PASS  GW-037: {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
