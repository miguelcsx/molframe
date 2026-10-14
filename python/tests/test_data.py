"""The table cache from Python: list, verify and fetch pass paths in and typed errors out."""

from pathlib import Path

import pytest

import molframe
from molframe import data

MANIFEST = """\
[[entry]]
name = "lab-radii"
version = "2"
url = "https://example.invalid/radii.radii"
licence = "CC0-1.0"
sha256 = "0000000000000000000000000000000000000000000000000000000000000000"
extension = "radii"
"""
LIMITS = data.DownloadLimits(max_bytes=1024, timeout_seconds=1, redirect_limit=0)


@pytest.fixture
def manifest(tmp_path: Path) -> Path:
    path = tmp_path / "manifest.toml"
    path.write_text(MANIFEST)
    return path


def test_list_reports_a_missing_entry_without_the_network(manifest: Path, tmp_path: Path):
    (row,) = data.list(manifest, tmp_path / "cache")
    assert (row.name, row.version, row.licence) == ("lab-radii", "2", "CC0-1.0")
    assert row.status == "missing"
    assert row.detail is None
    assert Path(row.path) == tmp_path / "cache" / "lab-radii.radii"


def test_verify_marks_a_cache_file_that_is_not_a_table_invalid(manifest: Path, tmp_path: Path):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "lab-radii.radii").write_text("not a table\n")
    (row,) = data.verify(manifest, cache, name="lab-radii")
    assert row.status == "invalid"
    assert row.detail


def test_an_unknown_entry_is_e7102_and_a_bad_manifest_is_e2201(manifest: Path, tmp_path: Path):
    with pytest.raises(molframe.MolframeError, match="E7102"):
        data.verify(manifest, tmp_path, name="nowhere")
    with pytest.raises(molframe.MolframeError, match="E7102"):
        data.fetch(manifest, "nowhere", LIMITS, tmp_path)
    bad = tmp_path / "bad.toml"
    bad.write_text(MANIFEST.replace("https://", "http://"))
    with pytest.raises(molframe.MolframeError, match="E2201"):
        data.list(bad, tmp_path)


def test_download_limits_are_explicit():
    assert (LIMITS.max_bytes, LIMITS.timeout_seconds, LIMITS.redirect_limit) == (1024, 1, 0)
    with pytest.raises(TypeError):
        data.DownloadLimits()  # type: ignore[call-arg]
