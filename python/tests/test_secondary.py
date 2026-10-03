from pathlib import Path

import pytest

import molframe


def test_explicit_dssp_requires_and_uses_caller_selected_ccd_roles():
    root = Path(__file__).resolve().parents[2]
    structure = molframe.read(root / "crates/molframe-bench/data/4hhb.cif")
    with pytest.raises(ValueError, match="polymer"):
        molframe.analysis.dssp(structure)
    rules = [
        molframe.chemistry.PolymerRoleRule(name, 1 << bit, component_kind=1)
        for bit, name in enumerate(("N", "CA", "C", "O"))
    ]
    report = molframe.chemistry.apply_polymer_role_profile(
        structure,
        root / "crates/molframe-chem/data/CCD-amino-acids.cif",
        rules,
        profile_id="wwpdb-backbone-2026-10-03",
        version="wwPDB-2026-10-03",
    )
    assert report.unresolved_components == ["HEM", "HOH", "PO4"]
    assert report.dictionary_version == "wwPDB-2026-10-03"
    assert report.profile_id == "wwpdb-backbone-2026-10-03"
    table = molframe.analysis.dssp(report.structure)
    kinds = list(table["kind"])
    assert {2, 4, 5, 6, 9}.issubset(kinds)
    empty = molframe.chemistry.apply_polymer_role_profile(
        structure,
        root / "crates/molframe-chem/data/CCD-amino-acids.cif",
        [molframe.chemistry.PolymerRoleRule("nonexistent", 1, component_kind=1)],
        profile_id="deliberately-unmatched-backbone",
    )
    assert set(molframe.analysis.dssp(empty.structure)["kind"]) == {0}


def test_secondary_codes_and_helix_membership_preserve_all_states():
    ss = molframe.SecondaryStructure
    states = [
        ss.Unknown,
        ss.Coil,
        ss.AlphaHelix,
        ss.Strand,
        ss.Turn,
        ss.ThreeTenHelix,
        ss.PiHelix,
        ss.OtherHelix,
        ss.BetaBridge,
        ss.Bend,
        ss.PolyProline,
    ]
    assert [int(state) for state in states] == list(range(11))
    assert [int(state) for state in states if state.is_helix()] == [2, 5, 6, 7, 10]
    assert [int(state) for state in states if state.is_strand()] == [3]
    assert [int(state) for state in states if state.is_sheet_like()] == [3, 8]
