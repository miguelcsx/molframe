"""Alternate conformations retain native coverage and stable atom indices."""

import molframe


def test_native_resolution_retains_blank_atoms_and_the_best_chain_conformer():
    pdb = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA AALA A   1       1.000   0.000   0.000  0.70  0.00           C\n"
        "ATOM      3  CA BALA A   1       2.000   0.000   0.000  0.30  0.00           C\n"
        "END\n"
    )
    structure = molframe.read(pdb.encode(), format="pdb")
    result = structure.resolve_altlocs(
        policy=molframe.AnalysisPolicy(altloc="conformer_consistent")
    )
    assert result.value is not None
    assert result.value.indices.tolist() == [0, 1]
    assert structure.atom_count == 3
