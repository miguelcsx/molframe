#!/usr/bin/env python3
"""Builds the reference tables embedded by molframe-validate.

Reads two downloaded sources and writes four small files next to this
script. Pure standard library; see README.md for provenance and licences.

    python3 -I build_reference.py \
        --top8000 DIR   checkout of rlabduke/reference_data
        --ost FILE      OpenStructure stereo_chemical_props.txt
        --out DIR       destination (default: this directory)
"""
import argparse
import math
import os
import struct
from itertools import groupby

RAMA_GRIDS = [
    ("general", "rama8000-general-noGPIVpreP.data"),
    ("gly", "rama8000-gly-sym.data"),
    ("transpro", "rama8000-transpro.data"),
    ("cispro", "rama8000-cispro.data"),
    ("prepro", "rama8000-prepro-noGP.data"),
    ("ileval", "rama8000-ileval-nopreP.data"),
]
ROTA_TABLES = ["arg", "asn", "asp", "cys", "gln", "glu", "his", "ile", "leu",
               "lys", "met", "phetyr", "pro", "ser", "thr", "trp", "val"]
ROTA_FAVOURED = 0.02
ROTA_ALLOWED = 0.003
LOG_FLOOR = -6.0
CHI_ATOMS = """\
ARG ARG N CA CB CG | CA CB CG CD | CB CG CD NE | CG CD NE CZ
ASN ASN N CA CB CG | CA CB CG OD1
ASP ASP N CA CB CG | CA CB CG OD1
CYS CYS N CA CB SG
GLN GLN N CA CB CG | CA CB CG CD | CB CG CD OE1
GLU GLU N CA CB CG | CA CB CG CD | CB CG CD OE1
HIS HIS N CA CB CG | CA CB CG ND1
ILE ILE N CA CB CG1 | CA CB CG1 CD1
LEU LEU N CA CB CG | CA CB CG CD1
LYS LYS N CA CB CG | CA CB CG CD | CB CG CD CE | CG CD CE NZ
MET MET N CA CB CG | CA CB CG SD | CB CG SD CE
PHE PHE N CA CB CG | CA CB CG CD1
PRO PRO N CA CB CG
SER SER N CA CB OG
THR THR N CA CB OG1
TRP TRP N CA CB CG | CA CB CG CD1
TYR PHE N CA CB CG | CA CB CG CD1
VAL VAL N CA CB CG1
"""


def read_grid(path):
    dims, cells = [], []
    for line in open(path):
        if line.startswith("#"):
            parts = line.split()
            if len(parts) == 6 and parts[1].startswith("x") and parts[1].endswith(":"):
                dims.append((float(parts[2]), float(parts[3]), int(parts[4])))
            continue
        parts = line.split()
        if parts:
            cells.append(([float(x) for x in parts[:-1]], float(parts[-1])))
    return dims, cells


def cell_index(dims, coords):
    index = 0
    for (lo, hi, bins), x in zip(dims, coords):
        width = (hi - lo) / bins
        i = int(math.floor((x - lo) / width))
        index = index * bins + min(max(i, 0), bins - 1)
    return index


def quantise(value):
    if value <= 0.0:
        return 0
    q = 1 + round(254 * (math.log10(value) - LOG_FLOOR) / -LOG_FLOOR)
    return min(max(q, 1), 255)


def build_rama(top, out):
    blob = bytearray(b"MFRM" + bytes([len(RAMA_GRIDS)]))
    folder = os.path.join(top, "Top8000", "Top8000_ramachandran_pct_contour_grids")
    for name, filename in RAMA_GRIDS:
        dims, cells = read_grid(os.path.join(folder, filename))
        assert [d[2] for d in dims] == [180, 180], dims
        grid = bytearray(180 * 180)
        for coords, value in cells:
            grid[cell_index(dims, coords)] = quantise(value)
        blob += name.encode().ljust(8, b"\0") + grid
    open(os.path.join(out, "rama.bin"), "wb").write(blob)


def varint(n):
    out = bytearray()
    while n >= 0x80:
        out.append(n & 0x7F | 0x80)
        n >>= 7
    out.append(n)
    return out


def build_rota(top, out):
    folder = os.path.join(top, "Top8000", "Top8000_rotamer_pct_contour_grids")
    blob = bytearray(b"MFRT" + bytes([len(ROTA_TABLES)]))
    for name in ROTA_TABLES:
        dims, cells = read_grid(os.path.join(folder, f"rota8000-{name}.data"))
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
            blob += struct.pack("<HB", bins, 1 if hi - lo < 359.0 else 0)
        blob += struct.pack("<I", len(payload)) + payload
    open(os.path.join(out, "rotamer.bin"), "wb").write(blob)


def build_chi(out):
    rows = ["# residue\ttable\tchi1 | chi2 ... (atom names)"]
    for line in CHI_ATOMS.strip().splitlines():
        head, _, tail = line.partition(" ")
        table, _, chis = tail.partition(" ")
        rows.append(f"{head}\t{table.lower().replace('phe', 'phetyr')}\t{chis}")
    open(os.path.join(out, "chi_atoms.tsv"), "w").write("\n".join(rows) + "\n")


# Peptide-link terms, which the OpenStructure table does not carry. C-N and the
# two angles at C and N come from the constants AlphaFold ships (derived there
# from Engh and Huber); O-C-N comes from the CCP4 Monomer Library link_TRANS
# and link_PTRANS entries. Residue PRO is the residue that follows the carbonyl.
LINK_ROWS = [
    ("bond", "*", "C-N", "1.329", "0.014", "AF"),
    ("bond", "PRO", "C-N", "1.341", "0.016", "AF"),
    ("angle", "*", "CA-C-N", "116.568", "1.995", "AF"),
    ("angle", "*", "C-N-CA", "121.352", "2.315", "AF"),
    ("angle", "*", "O-C-N", "123.469", "1.5", "CCP4"),
    ("angle", "PRO", "O-C-N", "121.016", "1.53", "CCP4"),
]


def build_geometry(ost, out):
    lines = [
        "# kind\tresidue\tatoms\tmean\tsd\tsource",
        "# EH = Engh and Huber 1991 via OpenStructure; AF = AlphaFold",
        "# residue_constants.py; CCP4 = CCP4 Monomer Library (see README.md).",
        "# Residue * = any residue. Peptide-link terms are keyed by the residue",
        "# that follows the carbonyl; the angle rows without a PRO entry apply to",
        "# every residue.",
    ]
    section = None
    for raw in open(ost):
        raw = raw.rstrip("\n")
        if raw.startswith("Bond"):
            section = "bond"
        elif raw.startswith("Angle"):
            section = "angle"
        elif raw.startswith("Non-bonded"):
            section = None
        elif section and raw.strip() and not raw.startswith("-"):
            atoms, residue, mean, sd = raw.split()
            lines.append(f"{section}\t{residue}\t{atoms}\t{mean}\t{sd}\tEH")
    for kind, residue, atoms, mean, sd, source in LINK_ROWS:
        lines.append(f"{kind}\t{residue}\t{atoms}\t{mean}\t{sd}\t{source}")
    open(os.path.join(out, "geometry.tsv"), "w").write("\n".join(lines) + "\n")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument("--top8000", required=True)
    parser.add_argument("--ost", required=True)
    parser.add_argument("--out", default=here)
    args = parser.parse_args()
    build_rama(args.top8000, args.out)
    build_rota(args.top8000, args.out)
    build_chi(args.out)
    build_geometry(args.ost, args.out)


if __name__ == "__main__":
    main()
