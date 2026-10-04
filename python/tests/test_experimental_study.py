import runpy
from pathlib import Path
from types import SimpleNamespace

import pytest


def workflow(monkeypatch, name):
    folder = Path(__file__).resolve().parents[2] / "studies/semantic-stress"
    monkeypatch.syspath_prepend(str(folder))
    return runpy.run_path(folder / name)


def test_related_experimental_complexes_stay_in_one_transitive_holdout_group(monkeypatch):
    cluster = workflow(monkeypatch, "prepare_skempi.py")["cluster_records"]
    records = [
        {"#Pdb": "A", "Hold_out_proteins": "A,B", "Protein 1": "first", "Protein 2": "second"},
        {"#Pdb": "B", "Hold_out_proteins": "B", "Protein 1": "third", "Protein 2": "fourth"},
        {"#Pdb": "C", "Hold_out_proteins": "C", "Protein 1": "third", "Protein 2": "fifth"},
        {"#Pdb": "D", "Hold_out_proteins": "D", "Protein 1": "unrelated", "Protein 2": "other"},
    ]
    groups = cluster(records)
    assert groups["A"] == groups["B"] == groups["C"]
    assert groups["D"] != groups["A"]


def test_experimental_mapping_refuses_wrong_wild_type_and_keeps_discordant_labels(monkeypatch):
    mapped = workflow(monkeypatch, "experimental_interfaces.py")["mapped_measurements"]
    residue = SimpleNamespace(index=7, name="ASP")
    records = [
        {"Mutation(s)_PDB": "DA10A", "ddg_kcal_mol": 1.0},
        {"Mutation(s)_PDB": "DA10A", "ddg_kcal_mol": 3.0},
        {"Mutation(s)_PDB": "WA10A", "ddg_kcal_mol": 4.0},
    ]
    rows, excluded = mapped({("A", 10, ""): [residue]}, records)
    assert excluded == {"wild_type_identity_mismatch": 1}
    assert len(rows) == 1
    assert rows[0]["hotspot"] is None
    assert rows[0]["ddg_median"] == 2.0
    assert len(rows[0]["measurements"]) == 2


def test_hotspot_threshold_is_fitted_from_training_examples_and_retains_single_class_groups(
    monkeypatch,
):
    fit = workflow(monkeypatch, "experimental_interfaces.py")["fit_threshold"]
    groups = [
        [{"score": 0, "hotspot": False}, {"score": 4, "hotspot": True}],
        [{"score": 2, "hotspot": False}, {"score": 6, "hotspot": True}],
        [{"score": 100, "hotspot": True}],
    ]
    assert fit(groups) == 3.0


def test_different_reference_temperatures_and_declared_experimental_identity_conflicts_are_refused(
    monkeypatch,
):
    conflict = workflow(monkeypatch, "prepare_skempi.py")["condition_conflict"]
    note = "wild type parameters measured at 303K"
    assert conflict({"Notes": note, "Temperature": "303"}) is None
    assert (
        conflict({"Notes": note, "Temperature": "279"}) == "reference_and_mutant_temperature_differ"
    )
    assert (
        conflict({"Notes": "In paper, wild-type AA is LYS, not ARG.", "Temperature": "298"})
        == "experimental_wild_type_identity_conflict"
    )


def test_energy_reference_requires_one_finite_result_and_ignores_runtime_text(monkeypatch):
    parse = workflow(monkeypatch, "matched_energy.py")["energy"]
    assert parse("Total = -24.94\nTotal time spent: 0.03\n") == -24.94
    for invalid in ("Total = nan\n", "Total = 1e999\n", "Total = 1\nTotal = 2\n", "no result"):
        with pytest.raises(ValueError, match="exactly one finite"):
            parse(invalid)
