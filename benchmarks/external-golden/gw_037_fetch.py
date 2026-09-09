#!/usr/bin/env python3
"""Fetch the pinned public 1CJB corpus used by GW-037."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import urllib.request


SOURCE = "https://files.rcsb.org/download/1CJB.cif"
EXPECTED_SHA256 = "b44610fb0a3d2c003dff40ca0d251e33232146de950439df5bb659a6d58ca3f8"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    with urllib.request.urlopen(SOURCE, timeout=30) as response:
        payload = response.read()
    observed = hashlib.sha256(payload).hexdigest()
    if observed != EXPECTED_SHA256:
        raise ValueError(f"1CJB SHA-256 differs: {observed}")
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_bytes(payload)
    print(f"FETCH  GW-037: {SOURCE}")
    print(f"SHA256 {observed}")
    print(f"OUTPUT {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
