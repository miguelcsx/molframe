#!/usr/bin/env python3
"""Generate pdbiox's compact offline space-group catalogue from spglib."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import spglib


MAGIC = b"PDBXSG01"
SETTING_COUNT = 530
DENOMINATOR = 12


def text(value: str) -> bytes:
    encoded = value.encode("utf-8")
    if len(encoded) > 255:
        raise ValueError(f"catalogue string is too long: {value!r}")
    return bytes([len(encoded)]) + encoded


def generate() -> bytes:
    output = bytearray(MAGIC)
    output.extend(struct.pack("<H", SETTING_COUNT))
    for hall_number in range(1, SETTING_COUNT + 1):
        kind = spglib.get_spacegroup_type(hall_number)
        symmetry = spglib.get_symmetry_from_database(hall_number)
        if kind is None or symmetry is None:
            raise RuntimeError(f"spglib has no Hall setting {hall_number}")
        rotations = symmetry["rotations"]
        translations = symmetry["translations"]
        if len(rotations) != len(translations):
            raise RuntimeError(f"operation columns differ for Hall setting {hall_number}")
        output.extend(struct.pack("<HH", hall_number, kind.number))
        output.extend(text(kind.international_short))
        output.extend(text(kind.international_full))
        output.extend(text(kind.hall_symbol))
        output.extend(text(kind.choice))
        output.extend(struct.pack("<H", len(rotations)))
        for rotation, translation in zip(rotations, translations, strict=True):
            for value in rotation.flat:
                integer = int(value)
                if integer not in (-1, 0, 1):
                    raise ValueError(f"unexpected rotation coefficient {integer}")
                output.extend(struct.pack("<b", integer))
            for value in translation:
                numerator = round(float(value) * DENOMINATOR) % DENOMINATOR
                if abs(float(value) - numerator / DENOMINATOR) > 1e-12:
                    raise ValueError(f"inexact translation {value}")
                output.append(numerator)
    return bytes(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--require-version", default="2.7.0")
    args = parser.parse_args()
    if spglib.__version__ != args.require_version:
        raise RuntimeError(
            f"expected spglib {args.require_version}, found {spglib.__version__}"
        )
    payload = generate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(payload)
    print(f"wrote {len(payload)} bytes to {args.output}")


if __name__ == "__main__":
    main()
