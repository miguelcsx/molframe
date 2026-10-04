"""Evaluate structural contact degree against independent SKEMPI measurements.

All policies use the same experimentally mutated author residue and the same
specified protein partners. Replicates and conditions are retained; a measured
residue is one statistical unit. Discordant hotspot labels remain ambiguous.
This exploratory pilot does not introduce or validate a binding-energy model.
"""

import argparse
import hashlib
import json
import re
import statistics
import warnings
from collections import Counter, defaultdict
from itertools import pairwise
from pathlib import Path

import numpy as np
from matched_interactions import AMINO_ACIDS, entity_kind

import molframe
from molframe import analysis

POLICIES = [
    (radii, slack)
    for radii in ("bondi", "charmm", "alvarez", "amber_united")
    for slack in (0.0, 0.25, 0.5)
]
BASELINE = "bondi:0.5"
ONE_TO_THREE = dict(
    zip(
        "ARNDCQEGHILKMFPSTWYV",
        [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ],
        strict=True,
    )
)


def balanced_accuracy(labels, predictions) -> float | None:
    positives, negatives = labels, ~labels
    if not positives.any() or not negatives.any():
        return None
    return float((predictions[positives].mean() + (~predictions[negatives]).mean()) / 2)


def auc(labels, scores) -> float | None:
    positive, negative = scores[labels], scores[~labels]
    if not len(positive) or not len(negative):
        return None
    differences = positive[:, None] - negative[None, :]
    return float(np.mean((differences > 0) + 0.5 * (differences == 0)))


def fit_threshold(groups) -> float | None:
    values = sorted({row["score"] for group in groups for row in group})
    if not values:
        return None
    thresholds = [
        values[0] - 1,
        *[(a + b) / 2 for a, b in pairwise(values)],
        values[-1] + 1,
    ]
    best, chosen = -1.0, None
    for threshold in thresholds:
        accuracies = []
        for group in groups:
            labels = np.array([row["hotspot"] for row in group], dtype=bool)
            predictions = np.array([row["score"] >= threshold for row in group])
            measured = balanced_accuracy(labels, predictions)
            if measured is not None:
                accuracies.append(measured)
        if accuracies and statistics.fmean(accuracies) > best:
            best, chosen = statistics.fmean(accuracies), threshold
    return chosen


def mapped_measurements(residue_map, measurements) -> tuple:
    mutations = defaultdict(list)
    for measurement in measurements:
        mutations[measurement["Mutation(s)_PDB"]].append(measurement)
    rows, exclusions = [], Counter()
    for mutation, replicates in sorted(mutations.items()):
        match = re.fullmatch(r"([A-Z])([A-Za-z0-9])(-?\d+)([A-Za-z]?)A", mutation)
        if match is None:
            exclusions["mutation_notation"] += 1
            continue
        wild, chain, number, insertion = match.groups()
        candidates = residue_map[(chain, int(number), insertion)]
        if len(candidates) != 1:
            exclusions["absent_or_ambiguous_author_residue"] += 1
            continue
        residue = candidates[0]
        if residue.name != ONE_TO_THREE[wild]:
            exclusions["wild_type_identity_mismatch"] += 1
            continue
        measured = [replicate["ddg_kcal_mol"] for replicate in replicates]
        labels = {value >= 2.0 for value in measured}
        rows.append(
            {
                "mutation": mutation,
                "residue": residue.index,
                "ddg_median": statistics.median(measured),
                "ddg_min": min(measured),
                "ddg_max": max(measured),
                "hotspot": next(iter(labels)) if len(labels) == 1 else None,
                "measurements": replicates,
                "scores": {},
            }
        )
    return rows, exclusions


def partner_residues(structure, groups) -> tuple:
    residue_map, partner = defaultdict(list), {}
    for residue in structure.residues:
        if residue.name not in AMINO_ACIDS or entity_kind(structure, residue) != "polymer":
            continue
        chain = residue.chain.auth_label
        if chain not in groups[0] | groups[1]:
            continue
        residue_map[(chain, residue.auth_number, residue.insertion_code or "")].append(residue)
        partner[residue.index] = 0 if chain in groups[0] else 1
    if not all(any(value == index for value in partner.values()) for index in (0, 1)):
        msg = "specified author-chain partners absent"
        raise ValueError(msg)
    return residue_map, partner


def namespace_ablation(structure, rows) -> dict:
    """Measure an actual wrong-namespace join, without assigning it to a runner."""
    labels = defaultdict(list)
    for residue in structure.residues:
        if residue.name in AMINO_ACIDS and entity_kind(structure, residue) == "polymer":
            labels[(residue.chain.label, residue.number, residue.insertion_code or "")].append(
                residue.index
            )
    counts = Counter()
    for row in rows:
        match = re.fullmatch(r"([A-Z])([A-Za-z0-9])(-?\d+)([A-Za-z]?)A", row["mutation"])
        if match is None:
            continue
        _, chain, number, insertion = match.groups()
        candidates = labels[(chain, int(number), insertion)]
        outcome = "same_residue" if candidates == [row["residue"]] else "wrong_or_missing_residue"
        counts[outcome] += 1
    return dict(counts)


def describe(path: Path, entry: dict, measurements: list[dict]) -> dict:
    structure = molframe.read(path, options=molframe.ReadOptions(digest_input=True))
    _, first_chain, second_chain = entry["complex"].split("_")
    groups = (set(first_chain), set(second_chain))
    residue_map, partner = partner_residues(structure, groups)
    rows, exclusions = mapped_measurements(residue_map, measurements)
    runs = []
    for radii, slack in POLICIES:
        key = f"{radii}:{slack}"
        policy = molframe.AnalysisPolicy(
            hydrogens="exclude", vdw_radii=radii, contact_def=f"distance:{slack}"
        )
        result = analysis.contacts_by_definition(structure, policy=policy)
        if result.value is None:
            runs.append({"policy": key, "refused": "indeterminate"})
            continue
        origins = np.asarray(result.atom_origin, dtype=np.int64)
        atom_residue = np.array([atom.residue.index for atom in structure.atoms])
        first = atom_residue[origins[np.asarray(result.value["first"], dtype=np.int64)]]
        second = atom_residue[origins[np.asarray(result.value["second"], dtype=np.int64)]]
        neighbours = defaultdict(set)
        for a, b in zip(first, second, strict=True):
            if a in partner and b in partner and partner[a] != partner[b]:
                neighbours[int(a)].add(int(b))
                neighbours[int(b)].add(int(a))
        for row in rows:
            row["scores"][key] = len(neighbours[row["residue"]])
        runs.append(
            {
                "policy": key,
                "estimand": "distinct opposite-partner residues in vdW contact with the experimentally mutated author residue",
                "analysed_atoms": len(origins),
            }
        )
    for row in rows:
        if len(row["scores"]) == len(POLICIES):
            row["scores"]["policy_median"] = statistics.median(row["scores"].values())
    return {
        **entry,
        "mapped_residues": len(rows),
        "mapping_exclusions": dict(exclusions),
        "label_namespace_ablation": namespace_ablation(structure, rows),
        "runs": runs,
        "residues": rows,
    }


def evaluate(entries: list[dict]) -> dict:
    from scipy.stats import spearmanr

    methods = [f"{radii}:{slack}" for radii, slack in POLICIES] + ["policy_median"]
    metrics = []
    for entry in entries:
        for method in methods:
            residues = [row for row in entry["residues"] if method in row["scores"]]
            if len(residues) < 3:
                continue
            scores = np.array([row["scores"][method] for row in residues], dtype=float)
            targets = np.array([row["ddg_median"] for row in residues])
            rho = (
                float(spearmanr(scores, targets).statistic)
                if np.ptp(scores) and np.ptp(targets)
                else None
            )
            classified = [row for row in residues if row["hotspot"] is not None]
            labels = np.array([row["hotspot"] for row in classified], dtype=bool)
            test_scores = np.array([row["scores"][method] for row in classified], dtype=float)
            training = [
                [
                    {"score": row["scores"][method], "hotspot": row["hotspot"]}
                    for row in other["residues"]
                    if method in row["scores"] and row["hotspot"] is not None
                ]
                for other in entries
                if other["cluster"] != entry["cluster"]
            ]
            threshold = fit_threshold(training)
            metrics.append(
                {
                    "complex": entry["complex"],
                    "cluster": entry["cluster"],
                    "method": method,
                    "residues": len(residues),
                    "ambiguous_labels": len(residues) - len(classified),
                    "hotspots": int(labels.sum()),
                    "spearman": rho,
                    "auc": auc(labels, test_scores),
                    "training_threshold": threshold,
                    "held_out_balanced_accuracy": balanced_accuracy(
                        labels, test_scores >= threshold
                    )
                    if threshold is not None
                    else None,
                }
            )
    summary = {}
    for method in methods:
        selected = [row for row in metrics if row["method"] == method]
        summary[method] = {}
        for metric in ("spearman", "auc", "held_out_balanced_accuracy"):
            values = [row[metric] for row in selected if row[metric] is not None]
            summary[method][metric] = {
                "complexes": len(values),
                "mean": statistics.fmean(values) if values else None,
            }
    paired = {}
    for metric in ("spearman", "auc", "held_out_balanced_accuracy"):
        baseline = {
            row["cluster"]: row[metric]
            for row in metrics
            if row["method"] == BASELINE and row[metric] is not None
        }
        alternative = {
            row["cluster"]: row[metric]
            for row in metrics
            if row["method"] == "policy_median" and row[metric] is not None
        }
        differences = np.array(
            [
                alternative[key] - baseline[key]
                for key in sorted(baseline.keys() & alternative.keys())
            ]
        )
        if len(differences):
            generator = np.random.default_rng(20261004)
            bootstrap = differences[
                generator.integers(0, len(differences), size=(10000, len(differences)))
            ].mean(axis=1)
            paired[metric] = {
                "clusters": len(differences),
                "mean_difference": float(differences.mean()),
                "paired_cluster_bootstrap_95_percent": list(
                    map(float, np.quantile(bootstrap, [0.025, 0.975]))
                ),
            }
    return {
        "methods": summary,
        "uncertainty_limit": "Bootstrap resamples paired test-cluster metrics conditional on fitted leave-cluster-out predictions; does not refit calibration or establish confirmatory significance.",
        "paired_policy_median_minus_baseline": paired,
        "by_complex": metrics,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    arguments = parser.parse_args()
    root = arguments.protocol.parent
    protocol = json.loads(arguments.protocol.read_text())
    acquisition = json.loads((root / "acquisition.json").read_text())
    digest = hashlib.sha256(arguments.protocol.read_bytes()).hexdigest()
    if digest != acquisition["protocol_sha256"]:
        msg = "sampling protocol changed since acquisition"
        raise ValueError(msg)
    for record in protocol["sources"] + acquisition["files"]:
        if "error" in record:
            continue
        if hashlib.sha256((root / record["path"]).read_bytes()).hexdigest() != record["sha256"]:
            msg = f"frozen input changed: {record['path']}"
            raise ValueError(msg)
    entries, refusals = [], []
    with warnings.catch_warnings(record=True) as diagnostics:
        warnings.simplefilter("always")
        for entry in protocol["selected"]:
            path = root / "deposited" / f"{entry['complex'].split('_')[0].lower()}.cif"
            measurements = [
                row for row in protocol["measurements"] if row["#Pdb"] == entry["complex"]
            ]
            try:
                entries.append(describe(path, entry, measurements))
            except (ValueError, RuntimeError) as error:
                refusals.append({**entry, "refused": str(error)})
    print(
        json.dumps(
            {
                "molframe": molframe.__version__,
                "protocol_sha256": digest,
                "endpoint": "Measured single-alanine binding disruption; geometric descriptor, not an energy model",
                "baseline": BASELINE,
                "cluster_limit": protocol["cluster_limit"],
                "entries": entries,
                "refusals": refusals,
                "diagnostics": sorted({str(item.message) for item in diagnostics}),
                "evaluation": evaluate(entries),
            },
            indent=1,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
