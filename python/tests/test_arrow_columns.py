"""Atom tables export categorical text as shared Arrow dictionaries.

Values are checked against gemmi on 4HHB and, when a copy is available
(``MOLFRAME_3J3Q``), on the 3J3Q capsid.
"""

import os
from pathlib import Path

import pyarrow as pa
import pyarrow.ipc
import pyarrow.parquet
import pytest

import molframe
from molframe import interop

gemmi = pytest.importorskip("gemmi")

ROOT = Path(__file__).resolve().parents[2]
FOUR_HHB = ROOT / "crates" / "molframe-bench" / "data" / "4hhb.cif"
CATEGORICAL = ("chain_id", "label_asym_id", "resname", "element_symbol", "atom_name", "altloc")


def _structures():
    cases = [("4hhb", FOUR_HHB)]
    big = os.environ.get("MOLFRAME_3J3Q")
    if big:
        cases.append(("3j3q", Path(big)))
    return cases


@pytest.fixture(params=_structures(), ids=lambda case: case[0], scope="module")
def case(request):
    _, path = request.param
    if not path.exists():
        pytest.skip(f"{path} is not available")
    return molframe.read(path), gemmi.read_structure(str(path))


def _gemmi_atoms(model):
    for chain in model:
        for residue in chain:
            for atom in residue:
                yield chain, residue, atom


def test_categorical_columns_are_dictionaries_and_values_match_gemmi(case):
    structure, reference = case
    table = pa.table(structure.atoms)
    for name in CATEGORICAL:
        assert pa.types.is_dictionary(table.schema.field(name).type), name
    assert "b_factor" in table.column_names
    chain = table.column("chain_id").to_pylist()
    resname = table.column("resname").to_pylist()
    name = table.column("atom_name").to_pylist()
    element = table.column("element_symbol").to_pylist()
    altloc = table.column("altloc").to_pylist()
    b_factor = table.column("b_factor").to_numpy()
    # gemmi groups atoms by chain; molframe keeps file order, so compare the
    # per-atom records as multisets.
    ours = sorted(
        (
            chain[i],
            resname[i],
            name[i],
            element[i].upper(),
            altloc[i] or "\0",
            round(b_factor[i], 2),
        )
        for i in range(table.num_rows)
    )
    theirs = sorted(
        (
            c.name,
            residue.name,
            atom.name,
            atom.element.name.upper(),
            atom.altloc,
            round(atom.b_iso, 2),
        )
        for c, residue, atom in _gemmi_atoms(reference[0])
    )
    assert len(ours) == len(theirs)
    assert ours == theirs


def test_dictionaries_are_small_and_shared_across_chunks(case):
    structure, _ = case
    batches = pa.table(structure.atoms).to_batches()
    for name in CATEGORICAL:
        dictionaries = [batch.column(name).dictionary for batch in batches]
        assert len(dictionaries[0]) < 5000, name
        assert all(d.equals(dictionaries[0]) for d in dictionaries)
        if len(batches) > 1:
            first = dictionaries[0].buffers()[1].address
            assert all(d.buffers()[1].address == first for d in dictionaries), name


def test_files_carry_schema_metadata_and_round_trip_through_every_engine(tmp_path, case):
    structure, _ = case
    parquet = tmp_path / "atoms.parquet"
    ipc = tmp_path / "atoms.arrow"
    interop.write_atoms_parquet(structure, parquet, metadata={"source": "fixture"})
    interop.write_atoms_ipc(structure, ipc, metadata={"source": "fixture"})
    expected = pa.table(structure.atoms)

    from_parquet = pyarrow.parquet.read_table(parquet)
    with pyarrow.ipc.open_file(ipc) as reader:
        from_ipc = reader.read_all()
    for table in (from_parquet, from_ipc):
        metadata = table.schema.metadata
        assert metadata[b"source"] == b"fixture"
        assert metadata[b"molframe:schema_version"] == b"2"
        assert b"angstrom" in metadata[b"molframe:units"]
        assert int(metadata[b"molframe:atom_count"]) == structure.atom_count
        for name in CATEGORICAL:
            assert table.column(name).to_pylist() == expected.column(name).to_pylist()

    polars = pytest.importorskip("polars")
    frame = polars.read_parquet(parquet)
    assert frame.height == expected.num_rows
    assert frame["chain_id"].to_list() == expected.column("chain_id").to_pylist()
    assert frame["resname"].n_unique() == len(set(expected.column("resname").to_pylist()))

    duckdb = pytest.importorskip("duckdb")
    rows = duckdb.sql(
        f"select chain_id, count(*) from read_parquet('{parquet}') "  # noqa: S608
        "where atom_name = 'CA' group by chain_id order by chain_id"
    ).fetchall()
    ca = list(_count(expected, "chain_id", "atom_name", "CA"))
    assert rows == ca
    on_arrow = duckdb.sql("select count(*) from from_ipc").fetchone()
    assert on_arrow == (expected.num_rows,)

    pandas = pytest.importorskip("pandas")
    categorical = from_parquet.to_pandas()["chain_id"]
    assert isinstance(categorical.dtype, pandas.CategoricalDtype)


def _count(table, group, column, value):
    counts = {}
    for g, v in zip(table.column(group).to_pylist(), table.column(column).to_pylist(), strict=True):
        if v == value:
            counts[g] = counts.get(g, 0) + 1
    return sorted(counts.items())
