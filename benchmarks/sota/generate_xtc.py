#!/usr/bin/env python3
"""Generate a deterministic large XTC stream and its minimal GRO topology."""

from __future__ import annotations

import argparse
import math
import os
import pathlib
import tempfile

import numpy as np
from mdtraj.formats import XTCTrajectoryFile


def positions(atoms: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    side = math.ceil(atoms ** (1.0 / 3.0))
    index = np.arange(atoms, dtype=np.int64)
    x = (index % side).astype(np.float32) * np.float32(0.01)
    y = ((index // side) % side).astype(np.float32) * np.float32(0.01)
    z = (index // (side * side)).astype(np.float32) * np.float32(0.01)
    return x, y, z


def write_topology(path: pathlib.Path, xyz: tuple[np.ndarray, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="ascii", dir=path.parent, delete=False
    ) as stream:
        temporary = pathlib.Path(stream.name)
        try:
            stream.write("pdbiox deterministic XTC benchmark\n")
            stream.write(f"{xyz[0].size}\n")
            for index, (x, y, z) in enumerate(zip(*xyz, strict=True)):
                serial = index % 99_999 + 1
                stream.write(
                    f"{serial:5d}{'RES':<5}{'C':>5}{serial:5d}"
                    f"{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}\n"
                )
            stream.write("  10.00000  10.00000  10.00000\n")
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)
    path.chmod(0o644)


def write_trajectory(
    path: pathlib.Path, xyz: tuple[np.ndarray, ...], frames: int
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = pathlib.Path(stream.name)
    try:
        coordinates = np.empty((1, xyz[0].size, 3), dtype=np.float32)
        coordinates[0, :, 0] = xyz[0]
        coordinates[0, :, 1] = xyz[1]
        coordinates[0, :, 2] = xyz[2]
        box = np.eye(3, dtype=np.float32)[None, :, :] * np.float32(10.0)
        with XTCTrajectoryFile(str(temporary), "w", force_overwrite=True) as writer:
            for frame in range(frames):
                shift = np.float32(frame * 0.001)
                coordinates[0, :, 0] = xyz[0] + shift
                writer.write(
                    coordinates,
                    time=np.asarray([frame * 2.0], dtype=np.float32),
                    step=np.asarray([frame], dtype=np.int32),
                    box=box,
                )
        os.replace(temporary, path)
        path.chmod(0o644)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("xtc", type=pathlib.Path)
    parser.add_argument("topology", type=pathlib.Path)
    parser.add_argument("--atoms", type=int, default=100_000)
    parser.add_argument("--frames", type=int, default=100)
    arguments = parser.parse_args()
    if arguments.atoms < 1 or arguments.frames < 1:
        parser.error("--atoms and --frames must be positive")
    if arguments.xtc.exists() or arguments.topology.exists():
        parser.error("refusing to overwrite an existing fixture")
    xyz = positions(arguments.atoms)
    write_topology(arguments.topology, xyz)
    try:
        write_trajectory(arguments.xtc, xyz, arguments.frames)
    except BaseException:
        arguments.topology.unlink(missing_ok=True)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
