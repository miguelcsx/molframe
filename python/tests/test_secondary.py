from pathlib import Path

import molframe


def test_dssp_without_a_role_profile_uses_the_standard_protein_roles():
    root = Path(__file__).resolve().parents[2]
    structure = molframe.read(root / "crates/molframe-bench/data/1crn.cif")
    table = molframe.analysis.dssp(structure)
    kinds = [int(kind) for kind in table["kind"]]
    assert len(kinds) == len(table["residue"]) > 0
    assert int(molframe.SecondaryStructure.AlphaHelix) in kinds
    assert int(molframe.SecondaryStructure.Strand) in kinds


def test_explicit_dssp_roles_agree_with_the_standard_roles():
    root = Path(__file__).resolve().parents[2]
    structure = molframe.read(root / "crates/molframe-bench/data/4hhb.cif")
    default = molframe.analysis.dssp(structure)
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
    assert list(default["kind"]) == kinds
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
