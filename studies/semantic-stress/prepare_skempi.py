"""Freeze an outcome-independent experimental interface pilot from SKEMPI 2.

The publisher XLS splits some fields at spaces. A pinned CSV mirror supplies
the readable schema; selected scientific values must match original XLS cells.
Only single alanine substitutions with numeric positive affinities and numeric
temperature are eligible. No censored measurement receives an invented value.
"""

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


def equivalent(first, second) -> bool:
    try:
        return math.isclose(float(first), float(second), rel_tol=1e-10)
    except ValueError:
        return str(first).replace(" ", "") == str(second).replace(" ", "")


def condition_conflict(record) -> str | None:
    note = record["Notes"]
    reference = re.search(r"wild type parameters measured at (\d+(?:\.\d+)?)K", note)
    if reference and not math.isclose(
        float(reference[1]), float(record["Temperature"]), abs_tol=0.01
    ):
        return "reference_and_mutant_temperature_differ"
    if "In paper, wild-type AA is" in note:
        return "experimental_wild_type_identity_conflict"
    return None


def cluster_records(records) -> dict:
    # Connected components use all source metadata, never observed predictions.
    parents = {record["#Pdb"]: record["#Pdb"] for record in records}

    def find(key):
        while parents[key] != key:
            key = parents[key]
        return key

    seen = {}
    for record in records:
        key = record["#Pdb"]
        tokens = record["Hold_out_proteins"].split(",")
        tokens += [
            f"protein:{record[field].strip().casefold()}"
            for field in ("Protein 1", "Protein 2")
            if record[field].strip()
        ]
        for token in tokens:
            if not token:
                continue
            if token in seen:
                a, b = sorted((find(key), find(seen[token])))
                parents[b] = a
            seen[token] = key
    return {key: find(key) for key in parents}


def select_pilot(eligible, parents) -> tuple:
    mutations = defaultdict(set)
    for record in eligible:
        mutations[record["#Pdb"]].add(record["Mutation(s)_PDB"])
    candidates = [key for key, values in mutations.items() if len(values) >= 10]

    def order(key):
        return hashlib.sha256(key.encode()).hexdigest()

    selected, clusters = [], set()
    for key in sorted(candidates, key=order):
        group = parents[key]
        if group not in clusters and len(selected) < 24:
            selected.append({"complex": key, "cluster": group, "mutations": len(mutations[key])})
            clusters.add(group)
    return candidates, selected


def prepare(root: Path) -> dict:
    import xlrd

    source = root / "skempi2-supplement.xls"
    mirror = root / "skempi_v2.csv"
    sheet = xlrd.open_workbook(source).sheet_by_index(0)
    with mirror.open() as stream:
        records = list(csv.DictReader(stream, delimiter=";"))
    if sheet.nrows != len(records) + 1:
        msg = "original supplement and CSV have different row counts"
        raise ValueError(msg)
    fields = [
        "#Pdb",
        "Mutation(s)_PDB",
        "Mutation(s)_cleaned",
        "iMutation_Location(s)",
        "Hold_out_type",
        "Hold_out_proteins",
        "Affinity_mut (M)",
        "Affinity_wt (M)",
        "Reference",
    ]
    eligible, exclusions, mismatches = [], Counter(), []
    for number, record in enumerate(records, start=1):
        cells = sheet.row_values(number)
        matches = all(equivalent(cells[i], record[field]) for i, field in enumerate(fields))
        temperature_column = 9 + sum(
            max(1, len(record[k].split())) for k in ("Protein 1", "Protein 2")
        )
        matches = matches and equivalent(cells[temperature_column], record["Temperature"])
        if not matches:
            mismatches.append(number + 1)
        if not re.fullmatch(r"[A-Z][A-Za-z0-9]-?\d+[A-Za-z]?A", record["Mutation(s)_PDB"]):
            exclusions["not_single_alanine"] += 1
            continue
        try:
            mutant, wild, temperature = (
                float(record[field])
                for field in ("Affinity_mut (M)", "Affinity_wt (M)", "Temperature")
            )
        except ValueError:
            exclusions["non_numeric_or_nonpositive_affinity_or_temperature"] += 1
            continue
        if not all(math.isfinite(value) and value > 0 for value in (mutant, wild, temperature)):
            exclusions["non_numeric_or_nonpositive_affinity_or_temperature"] += 1
            continue
        if not matches:
            exclusions["original_supplement_mismatch"] += 1
            continue
        conflict = condition_conflict(record)
        if conflict is not None:
            exclusions[conflict] += 1
            continue
        eligible.append(
            {
                **record,
                "source_row": number + 1,
                "ddg_kcal_mol": 0.00198720425864083 * temperature * math.log(mutant / wild),
            }
        )

    parents = cluster_records(records)
    candidates, selected = select_pilot(eligible, parents)
    return {
        "sources": [
            {"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in (source, mirror)
        ],
        "source_rows": len(records),
        "concordant_rows": len(records) - len(mismatches),
        "mismatched_source_rows": mismatches,
        "exclusions": dict(exclusions),
        "eligible_measurements": len(eligible),
        "candidate_complexes": len(candidates),
        "selection": "At least ten distinct eligible alanine substitutions; SHA256(complex) order; at most one per metadata/shared-protein component; cap 24. No score-dependent selection.",
        "cluster_limit": "Author holdout links and identical protein names, not independently validated sequence-family clusters.",
        "selected": selected,
        "measurements": [
            record
            for record in eligible
            if record["#Pdb"] in {entry["complex"] for entry in selected}
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    arguments = parser.parse_args()
    print(json.dumps(prepare(arguments.inputs), indent=1, allow_nan=False))


if __name__ == "__main__":
    main()
