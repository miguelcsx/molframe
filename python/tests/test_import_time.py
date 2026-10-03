"""Importing the package stays under the 50 ms specification budget."""

from __future__ import annotations

import statistics
import subprocess
import sys

BUDGET_SECONDS = 0.050
RUNS = 5
SNIPPET = "import time;t=time.perf_counter();import molframe;print(time.perf_counter()-t)"


def _import_seconds() -> float:
    completed = subprocess.run(
        [sys.executable, "-c", SNIPPET],
        check=True,
        capture_output=True,
        text=True,
    )
    return float(completed.stdout.strip())


def test_import_is_under_budget() -> None:
    median = statistics.median(_import_seconds() for _ in range(RUNS))
    assert median < BUDGET_SECONDS, f"median import took {median * 1000:.1f} ms"
