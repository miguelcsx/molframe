"""The structure-level SASA options, read through the governed Python entry point."""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import analysis

BENCH = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"

# FreeSASA 2.2.1, Shrake-Rupley, Bondi radii, 1.4 angstrom probe.
CRAMBIN = 2970.152
HEMOGLOBIN_POLYMER = 26097.315
HEMOGLOBIN_ALL = 24073.036


def _total(result) -> float:
    return float(np.nansum(result.value))


def test_the_default_reproduces_freesasa_on_crambin():
    structure = molframe.read(BENCH / "1crn.cif")
    assert _total(analysis.sasa(structure)) == pytest.approx(CRAMBIN, rel=5e-3)


def test_hemoglobin_measures_the_polymer_alone_by_default_and_every_heteroatom_on_request():
    structure = molframe.read(BENCH / "4hhb.pdb")
    default = analysis.sasa(structure)
    assert _total(default) == pytest.approx(HEMOGLOBIN_POLYMER, rel=5e-3)
    assert int(np.isnan(default.value).sum()) > 0
    everything = analysis.sasa(structure, hetero="all")
    assert _total(everything) == pytest.approx(HEMOGLOBIN_ALL, rel=5e-3)
    listed = analysis.sasa(structure, hetero="waters,ligands,ions")
    assert np.array_equal(listed.value, everything.value, equal_nan=True)


def test_lee_richards_reads_the_same_atoms_as_shrake_rupley():
    structure = molframe.read(BENCH / "1crn.cif")
    shrake = analysis.sasa(structure)
    lee = analysis.sasa(structure, method="lee-richards")
    assert np.array_equal(np.isnan(shrake.value), np.isnan(lee.value))
    assert _total(lee) == pytest.approx(_total(shrake), rel=2e-2)


def test_unknown_words_are_refused_with_the_choices():
    structure = molframe.read(BENCH / "1crn.cif")
    with pytest.raises(ValueError, match="shrake-rupley or lee-richards"):
        analysis.sasa(structure, method="fast")
    with pytest.raises(ValueError, match="waters, ligands, ions"):
        analysis.sasa(structure, hetero="lipids")
    with pytest.raises(ValueError, match="none, vdw or ionic"):
        analysis.sasa(structure, metals="copper")
