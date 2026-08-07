#!/usr/bin/env python3
"""Consumer proof for the Arrow PyCapsule boundary.

Run against an installed pdbiox wheel with pyarrow, polars and pandas present.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("POLARS_UNKNOWN_EXTENSION_TYPE_BEHAVIOR", "load_as_storage")

import pandas as pd
import polars as pl
import pyarrow as pa

import pdbiox


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: check-arrow-python.py STRUCTURE")
        return 2

    structure = pdbiox.read(Path(sys.argv[1]))
    atom_table = pa.table(structure.atoms)
    coordinates = atom_table.column("coordinates").chunk(0)
    arrow_address = coordinates.values.buffers()[1].address
    numpy_address = structure.xyz.__array_interface__["data"][0]

    assert arrow_address == numpy_address
    assert atom_table.num_rows == structure.atom_count
    assert pa.table(structure.residues).num_rows == structure.residue_count
    assert pa.table(structure.chains).num_rows == structure.chain_count
    assert pa.table(structure.bonds).num_rows == len(structure.bonds)

    polars_frame = pl.DataFrame(structure.atoms)
    pandas_frame = pa.table(structure.atoms).to_pandas()
    assert polars_frame.height == structure.atom_count
    assert isinstance(pandas_frame, pd.DataFrame)
    assert len(pandas_frame) == structure.atom_count

    print(
        "Arrow consumers passed: "
        f"rows={structure.atom_count}, coordinates_zero_copy=true"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
