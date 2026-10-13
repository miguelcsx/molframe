"""Builds the Top8000 Ramachandran and rotamer tables molframe reads.

Reads a checkout of the Richardson Laboratory reference data and writes two
files in the versioned table format. Pure standard library; see README.md for
provenance and licences.

    python3 -I build_reference.py \
        --top8000 DIR   checkout of rlabduke/reference_data
        --commit SHA    the commit that checkout is at
        --out DIR       destination (default: this directory)
"""  # noqa: INP001

import argparse
import math
import struct
from itertools import groupby
from pathlib import Path

RAMA_GRIDS = [
    ("general", "rama8000-general-noGPIVpreP.data"),
    ("gly", "rama8000-gly-sym.data"),
    ("transpro", "rama8000-transpro.data"),
    ("cispro", "rama8000-cispro.data"),
    ("prepro", "rama8000-prepro-noGP.data"),
    ("ileval", "rama8000-ileval-nopreP.data"),
]
ROTA_TABLES = [
    "arg",
    "asn",
    "asp",
    "cys",
    "gln",
    "glu",
    "his",
    "ile",
    "leu",
    "lys",
    "met",
    "phetyr",
    "pro",
    "ser",
    "thr",
    "trp",
    "val",
]
ROTA_FAVOURED = 0.02
ROTA_ALLOWED = 0.003
LOG_FLOOR = -6.0
RAMA_BINS = 180
GRID_HEADER_FIELDS = 6
VARINT_MORE = 0x80
FULL_CIRCLE = 359.0
SOURCE = (
    "Top8000 {} contour grids, Richardson Laboratory, Duke University "
    "(https://github.com/rlabduke/reference_data commit {}); Williams et al., "
    "Protein Sci 27:293 (2018); Hintze et al., Proteins 84:1177 (2016); "
    "log-quantised and re-encoded by molframe"
)
RAMA_CUTOFFS = "favoured: 0.02\noutlier: 0.0005 0.001 0.001 0.002 0.001 0.001\n"


def with_header(magic: bytes, name: str, kind: str, commit: str, extra: str) -> bytes:
    """Magic, header length and the identity header of a binary table."""
    header = (
        "# SPDX-License-Identifier: CC-BY-4.0\n"
        "# SPDX-FileCopyrightText: 2026 Richardson Laboratory, Duke University; "
        "reformatted by Miguel Cardenas\n"
        f"format: 1\nname: {name}\nversion: {commit[:8]}\n"
        f"source: {SOURCE.format(kind, commit)}\nlicence: CC-BY-4.0\n{extra}"
    ).encode()
    return magic + struct.pack("<H", len(header)) + header


Grid = tuple[list[tuple[float, float, int]], list[tuple[list[float], float]]]


def read_grid(path: Path) -> Grid:
    """Parse one contour grid into its axis (low, high, bins) and its cells."""
    dims: list[tuple[float, float, int]] = []
    cells: list[tuple[list[float], float]] = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if line.startswith("#"):
            if (
                len(parts) == GRID_HEADER_FIELDS
                and parts[1].startswith("x")
                and parts[1].endswith(":")
            ):
                dims.append((float(parts[2]), float(parts[3]), int(parts[4])))
        elif parts:
            cells.append(([float(x) for x in parts[:-1]], float(parts[-1])))
    return dims, cells


def cell_index(dims: list[tuple[float, float, int]], coords: list[float]) -> int:
    """Row-major index of the grid cell holding `coords`, clamped to the grid."""
    index = 0
    for (lo, hi, bins), x in zip(dims, coords, strict=True):
        width = (hi - lo) / bins
        i = math.floor((x - lo) / width)
        index = index * bins + min(max(i, 0), bins - 1)
    return index


def quantise(value: float) -> int:
    """Map a density to one byte on a log scale; zero stays zero."""
    if value <= 0.0:
        return 0
    q = 1 + round(254 * (math.log10(value) - LOG_FLOOR) / -LOG_FLOOR)
    return min(max(q, 1), 255)


def build_rama(top: Path, out: Path, commit: str) -> None:
    """Write `top8000-rama.rama`: one quantised 180 x 180 grid per category."""
    blob = bytearray(with_header(b"MFR2", "top8000-rama", "Ramachandran", commit, RAMA_CUTOFFS))
    blob.append(len(RAMA_GRIDS))
    folder = top / "Top8000" / "Top8000_ramachandran_pct_contour_grids"
    for name, filename in RAMA_GRIDS:
        dims, cells = read_grid(folder / filename)
        if [d[2] for d in dims] != [RAMA_BINS, RAMA_BINS]:
            message = f"unexpected Ramachandran grid shape {dims}"
            raise ValueError(message)
        grid = bytearray(RAMA_BINS * RAMA_BINS)
        for coords, value in cells:
            grid[cell_index(dims, coords)] = quantise(value)
        blob += name.encode().ljust(8, b"\0") + grid
    (out / "top8000-rama.rama").write_bytes(blob)


def varint(n: int) -> bytearray:
    """LEB128 encoding of a non-negative integer."""
    out = bytearray()
    while n >= VARINT_MORE:
        out.append(n & 0x7F | VARINT_MORE)
        n >>= 7
    out.append(n)
    return out


def build_rota(top: Path, out: Path, commit: str) -> None:
    """Write `top8000-rotamer.rotamer`: run-length coded class grids."""
    folder = top / "Top8000" / "Top8000_rotamer_pct_contour_grids"
    blob = bytearray(with_header(b"MFT2", "top8000-rotamer", "rotamer", commit, ""))
    blob.append(len(ROTA_TABLES))
    for name in ROTA_TABLES:
        dims, cells = read_grid(folder / f"rota8000-{name}.data")
        total = 1
        for _, _, bins in dims:
            total *= bins
        classes = bytearray(total)
        for coords, value in cells:
            classes[cell_index(dims, coords)] = (
                2 if value >= ROTA_FAVOURED else 1 if value >= ROTA_ALLOWED else 0
            )
        payload = bytearray()
        for cls, run in groupby(classes):
            payload += varint(len(list(run)) << 2 | cls)
        label = name.encode()
        blob += bytes([len(label)]) + label + bytes([len(dims)])
        for lo, hi, bins in dims:
            blob += struct.pack("<HB", bins, 1 if hi - lo < FULL_CIRCLE else 0)
        blob += struct.pack("<I", len(payload)) + payload
    (out / "top8000-rotamer.rotamer").write_bytes(blob)


def main() -> None:
    """Parse the command line and write both tables."""
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--top8000", required=True, type=Path)
    parser.add_argument("--commit", required=True, help="reference_data commit")
    parser.add_argument("--out", default=here, type=Path)
    args = parser.parse_args()
    build_rama(args.top8000, args.out, args.commit)
    build_rota(args.top8000, args.out, args.commit)


if __name__ == "__main__":
    main()
