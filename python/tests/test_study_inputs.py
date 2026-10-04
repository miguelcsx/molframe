import hashlib
import runpy
from pathlib import Path

import pytest

VERIFY_FILES = runpy.run_path(
    Path(__file__).resolve().parents[2] / "studies/semantic-stress/build_corpus.py"
)["verify_files"]


def test_frozen_study_inputs_accept_matching_bytes_and_refuse_changed_bytes(tmp_path):
    source = tmp_path / "protein.pdb"
    source.write_bytes(b"original")
    records = [{"path": source.name, "sha256": hashlib.sha256(b"original").hexdigest()}]
    VERIFY_FILES(tmp_path, records)

    source.write_bytes(b"modified")
    with pytest.raises(ValueError, match="frozen input changed"):
        VERIFY_FILES(tmp_path, records)


def test_missing_frozen_study_inputs_are_refused(tmp_path):
    with pytest.raises(FileNotFoundError):
        VERIFY_FILES(tmp_path, [{"path": "missing.pdb", "sha256": "0" * 64}])
