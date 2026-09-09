#!/usr/bin/env python3
"""pdbiox side of GW-043 using only public trajectory/geometry APIs."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pdbiox


SELECTION_SIZE = 4096


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def selected_frames(path: Path) -> np.ndarray:
    reader = pdbiox.traj.read_trajectory(path)
    if reader.n_atoms < SELECTION_SIZE:
        raise ValueError(f"fixture has {reader.n_atoms} atoms, expected at least {SELECTION_SIZE}")
    frames: list[np.ndarray] = []
    while (frame := reader.read_next()) is not None:
        positions = np.asarray(frame.positions, dtype=np.float32)
        frames.append(positions[:SELECTION_SIZE])
    if not frames:
        raise ValueError("fixture contains no trajectory frames")
    return np.stack(frames)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    arguments = parser.parse_args()
    frames = selected_frames(arguments.fixture)
    rmsd = np.asarray(
        pdbiox.traj.rmsd_to_reference(frames, 0, pdbiox.traj.FrameAlignment.Rigid),
        dtype=np.float64,
    )
    reference = frames[0]
    aligned = np.empty_like(frames)
    for index, frame in enumerate(frames):
        aligned[index] = pdbiox.superpose(frame, reference).transform.apply_all(frame)
    rmsf = np.asarray(pdbiox.rmsf(aligned), dtype=np.float64)
    report = {
        "schema_version": 1,
        "workflow": "GW-043",
        "implementation": "pdbiox",
        "version": importlib.metadata.version("pdbiox"),
        "fixture_sha256": sha256(arguments.fixture),
        "units": "angstrom",
        "policy": {
            "selection": f"atom_index < {SELECTION_SIZE}",
            "alignment": "unweighted rigid Kabsch to frame 0",
            "rmsf": "per-atom population RMSF after frame alignment",
        },
        "frame_count": int(frames.shape[0]),
        "atom_count": int(frames.shape[1]),
        "rmsd": rmsd.tolist(),
        "rmsf": rmsf.tolist(),
    }
    print(json.dumps(report, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
