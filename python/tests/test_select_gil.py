"""Selection releases the GIL, honours cancellation, and name tests stay cheap."""

import threading
import time

import numpy as np
import pytest

import molframe

RESIDUES = 20_000
NAMES = ("N", "CA", "C", "O", "CB")


def _large_structure() -> molframe.Structure:
    """Build a synthetic 100k-atom protein-like chain, built without any data file."""
    header = (
        "data_big\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n"
    )
    rows = []
    serial = 0
    for residue in range(1, RESIDUES + 1):
        for slot, name in enumerate(NAMES):
            serial += 1
            x = (residue % 100) * 3.8 + slot * 0.3
            y = ((residue // 100) % 100) * 3.8
            z = (residue // 10_000) * 3.8
            coords = f"{x:.2f} {y:.2f} {z:.2f}"
            rows.append(f"ATOM {serial} {name[0]} {name} ALA A {residue} {coords} A {residue}")
    return molframe.read("\n".join([header, *rows, ""]).encode(), name="big.cif")


@pytest.fixture(scope="module")
def big() -> molframe.Structure:
    return _large_structure()


def _ticks_during(work) -> tuple[float, int]:
    """Run `work` while a thread counts 1 ms ticks; return (seconds, ticks)."""
    ticks = 0
    stop = threading.Event()

    def tick() -> None:
        nonlocal ticks
        while not stop.is_set():
            ticks += 1
            time.sleep(0.001)

    ticker = threading.Thread(target=tick)
    ticker.start()
    try:
        start = time.perf_counter()
        work()
        elapsed = time.perf_counter() - start
    finally:
        stop.set()
        ticker.join()
    return elapsed, ticks


def test_a_long_selection_lets_other_python_threads_run(big):
    query = "within 6 of (resid 1:3000)"
    big.select(query)  # warm up
    elapsed, ticks = _ticks_during(lambda: [big.select(query) for _ in range(5)])
    # With the GIL held for the whole call the ticker would see almost nothing.
    assert elapsed > 0.05, "the selection is too fast to observe"
    assert ticks >= 0.2 * elapsed * 1000


def test_a_selection_is_cancelled_by_its_context(big):
    context = molframe.ExecutionContext()
    timer = threading.Timer(0.01, context.cancel)
    timer.start()

    def run_until_cancelled() -> None:
        for _ in range(10_000):
            big.select("within 6 of (resid 1:3000)", context=context)

    try:
        with pytest.raises(molframe.Cancelled):
            run_until_cancelled()
    finally:
        timer.cancel()


def test_a_selection_under_a_cancelled_context_is_refused(big):
    context = molframe.ExecutionContext()
    context.cancel()
    with pytest.raises(molframe.Cancelled):
        big.select("name CA", context=context)
    with pytest.raises(molframe.Cancelled):
        molframe.Query("name CA").select(big, context=context)


def test_select_accepts_a_context_without_changing_the_answer(big):
    context = molframe.ExecutionContext(workers=1)
    plain = big.select("name CA")
    governed = big.select("name CA", context=context)
    assert np.array_equal(plain.indices, governed.indices)
    assert len(plain.indices) == RESIDUES


def test_single_name_selections_match_a_reference_scan(big):
    expected_ca = np.arange(RESIDUES, dtype=np.uint32) * len(NAMES) + 1
    assert np.array_equal(big.select("name CA").indices, expected_ca)
    both = big.select("name N O").indices
    assert np.array_equal(both, np.sort(np.concatenate([expected_ca - 1, expected_ca + 2])))
    assert len(big.select("name CA and resid 1:10").indices) == 10
    assert len(big.select("not name CA").indices) == RESIDUES * (len(NAMES) - 1)
    assert len(big.select("name C*").indices) == RESIDUES * 3


def test_a_name_selection_is_not_slower_than_a_compound_one(big):
    big.select("name CA")
    start = time.perf_counter()
    for _ in range(5):
        big.select("name CA")
    single = time.perf_counter() - start
    start = time.perf_counter()
    for _ in range(5):
        big.select("resname ALA and name CA")
    compound = time.perf_counter() - start
    assert single < 5 * compound + 0.05
