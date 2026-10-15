"""The default flows run on a clean machine: no data directory, no profile, no env.

The suite's conftest points the process at the repository's tables, and the table
store reads the environment once per process, so the clean host is a subprocess
with its own empty home and every ``MOLFRAME_`` variable removed. What answers
there is only what the wheel compiles in.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

CRAMBIN = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data" / "1crn.cif"

SNIPPET = """
import json, sys
import numpy as np
import molframe
from molframe import analysis

structure = molframe.read(sys.argv[1])
read_bonds = structure.bond_count
structure.infer_bonds()
area = analysis.sasa(structure)
under_policy = analysis.sasa(structure, policy=molframe.AnalysisPolicy())
residues = analysis.residue_sasa(structure)
print(json.dumps({
    "read_bonds": read_bonds,
    "inferred_bonds": structure.bond_count,
    "selected": len(structure.select("protein and name CA")),
    "area": float(np.nansum(area.value)),
    "policy_area": float(np.nansum(under_policy.value)),
    "residues": len(residues.value),
    "dssp": len(analysis.dssp(structure)),
    "assigned": sum(1 for state in structure.secondary_structure if state is not None),
    "contacts": len(analysis.contacts(structure, 4.0).value),
}))
"""


def _clean_environment(home: Path) -> dict[str, str]:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("MOLFRAME_")
    }
    environment["HOME"] = str(home)
    environment["XDG_CACHE_HOME"] = str(home / ".cache")
    environment["USERPROFILE"] = str(home)
    environment["LOCALAPPDATA"] = str(home / "AppData" / "Local")
    return environment


def test_the_default_flows_need_no_data_directory_or_profile(tmp_path: Path) -> None:
    """Reading, bonds, selection, surface area and secondary structure all answer."""
    completed = subprocess.run(
        [sys.executable, "-c", SNIPPET, str(CRAMBIN)],
        capture_output=True,
        check=False,
        cwd=tmp_path,
        env=_clean_environment(tmp_path),
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    found = json.loads(completed.stdout)
    assert found["read_bonds"] > 0
    assert found["inferred_bonds"] > 0
    assert found["selected"] == 46
    assert found["area"] > 0
    assert found["policy_area"] == found["area"]
    assert found["residues"] > 0
    assert found["dssp"] > 0
    assert found["contacts"] > 0
    assert not (tmp_path / ".cache").exists(), "nothing was cached"
