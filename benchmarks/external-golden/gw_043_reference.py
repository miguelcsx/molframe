#!/usr/bin/env python3
"""MDAnalysis + MDTraj side of GW-043 under the policy used by pdbiox."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.coordinates.XTC import XTCReader
import mdtraj as md
from mdtraj.formats import XTCTrajectoryFile
import numpy as np


SELECTION_SIZE = 4096


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_mdanalysis(path: Path) -> np.ndarray:
    reader = XTCReader(str(path), convert_units=True)
    if reader.n_atoms < SELECTION_SIZE:
        raise ValueError(f"fixture has {reader.n_atoms} atoms, expected at least {SELECTION_SIZE}")
    frames = [timestep.positions[:SELECTION_SIZE].copy() for timestep in reader]
    reader.close()
    return np.stack(frames).astype(np.float32, copy=False)


def load_mdtraj(path: Path) -> np.ndarray:
    batches: list[np.ndarray] = []
    with XTCTrajectoryFile(str(path), "r") as source:
        while True:
            coordinates, _, _, _ = source.read(n_frames=16)
            if coordinates.shape[0] == 0:
                break
            batches.append(coordinates[:, :SELECTION_SIZE, :])
    return np.concatenate(batches).astype(np.float32, copy=False) * np.float32(10.0)


def mdanalysis_results(frames: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    universe = mda.Universe.empty(SELECTION_SIZE, trajectory=True)
    universe.load_new(frames, format=mda.coordinates.memory.MemoryReader)
    align.AlignTraj(
        universe,
        universe,
        select="all",
        ref_frame=0,
        in_memory=True,
        match_atoms=False,
    ).run()
    aligned = np.stack([timestep.positions.copy() for timestep in universe.trajectory])
    rmsd = np.asarray(
        [rms.rmsd(frame, aligned[0], center=False, superposition=False) for frame in aligned]
    )
    universe.trajectory.rewind()
    rmsf = rms.RMSF(universe.atoms).run().results.rmsf
    return rmsd, np.asarray(rmsf)


def mdtraj_results(frames_angstrom: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    trajectory = md.Trajectory(frames_angstrom * np.float32(0.1), topology=None)
    rmsd = md.rmsd(trajectory, trajectory, 0, precentered=False).astype(np.float64) * 10.0
    trajectory.superpose(trajectory, 0)
    rmsf = md.rmsf(trajectory, trajectory, 0, precentered=False).astype(np.float64) * 10.0
    return rmsd, rmsf


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    arguments = parser.parse_args()
    mda_frames = load_mdanalysis(arguments.fixture)
    mdtraj_frames = load_mdtraj(arguments.fixture)
    coordinate_delta = float(np.max(np.abs(mda_frames - mdtraj_frames)))
    if coordinate_delta > 2.0e-5:
        raise ValueError(f"MDAnalysis/MDTraj coordinate mismatch: {coordinate_delta}")
    mda_rmsd, mda_rmsf = mdanalysis_results(mda_frames)
    mdtraj_rmsd, mdtraj_rmsf = mdtraj_results(mdtraj_frames)
    report = {
        "schema_version": 1,
        "workflow": "GW-043",
        "implementation": "MDAnalysis+MDTraj",
        "versions": {
            "MDAnalysis": importlib.metadata.version("MDAnalysis"),
            "mdtraj": importlib.metadata.version("mdtraj"),
        },
        "fixture_sha256": sha256(arguments.fixture),
        "units": "angstrom",
        "policy": {
            "selection": f"atom_index < {SELECTION_SIZE}",
            "alignment": "unweighted rigid Kabsch to frame 0",
            "rmsf": "per-atom population RMSF after frame alignment",
        },
        "frame_count": int(mda_frames.shape[0]),
        "atom_count": int(mda_frames.shape[1]),
        "reader_max_abs_delta": coordinate_delta,
        "mdanalysis": {"rmsd": mda_rmsd.tolist(), "rmsf": mda_rmsf.tolist()},
        "mdtraj": {"rmsd": mdtraj_rmsd.tolist(), "rmsf": mdtraj_rmsf.tolist()},
    }
    print(json.dumps(report, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
