"""How buried a ligand is, under the decisions that surface area depends on.

For every complex of the frozen corpus, a ligand atom is *buried* when its solvent-accessible
area in the complex is under 1 square angstrom, and the ligand is *mostly buried* when more
than half of its atoms are. Both are computed under every defensible combination of:

- whether the hydrogens the file carries occlude (``hydrogens``: interpretive),
- which van der Waals radii the atoms are given (``vdw_radii``: algorithmic).

United-atom radii already contain the hydrogens, so pairing them with modelled hydrogens
counts the hydrogens twice. That combination is forbidden, with the reason recorded, and the
audit runs a constrained plan: its attribution is Shapley's, not the product's.

The threshold (1 square angstrom) and the fraction (one half) are parameters of the
question, not decisions of the policy, and are fixed here rather than varied.

Run from the repository root:

    PYTHONPATH=python uv run --no-sync python studies/semantic-stress/burial.py \
        --refined-set /path/to/refined-set > studies/semantic-stress/burial.json
"""

import argparse
import json
import statistics
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
from build_corpus import sha256, verify_files
from pockets import combine

import molframe
from molframe import analysis, audit

warnings.simplefilter("ignore")

BURIED_BELOW = 1.0
MOSTLY = 0.5

SPACE = (
    audit.PolicySpace()
    .vary(
        "hydrogens",
        ["explicit_only", "exclude"],
        rationale="prepared protein files carry added hydrogens; deposited entries do not",
        evidence="PDBbind v2020 protein files are about half hydrogen",
    )
    .vary(
        "vdw_radii",
        ["bondi", "charmm", "alvarez", "amber_united"],
        rationale="four radius sets ship with the library; none is claimed to be correct",
        evidence="molframe-chem radius tables",
    )
    .forbid(
        ("hydrogens", "explicit_only"),
        ("vdw_radii", "amber_united"),
    )
)
PLAN = SPACE.plan(constrained=True)


def run_entry(refined: Path, entry: dict, scratch: Path) -> dict:
    code = entry["id"]
    combined = scratch / f"{code}.pdb"
    keys = combine(
        refined / entry["files"]["protein"]["path"],
        refined / entry["files"]["ligand"]["path"],
        combined,
    )
    n_protein = len(keys)
    structure = molframe.read(combined, options=molframe.ReadOptions(digest_input=True))
    cache: dict[str, object] = {}

    def analyse(policy):
        # Both audits below ask the same analysis of the same universe; run it once.
        key = repr(policy)
        if key not in cache:
            cache[key] = analysis.sasa(structure, policy=policy)
        return cache[key]

    def ligand_areas(analysed) -> np.ndarray:
        origin = np.asarray(analysed.atom_origin, dtype=np.int64)
        return np.asarray(analysed.value)[origin >= n_protein]

    def buried_fraction(analysed) -> float:
        # An atom with no radius has no area and is not counted as buried or exposed.
        areas = ligand_areas(analysed)
        areas = areas[np.isfinite(areas)]
        return float(np.mean(areas < BURIED_BELOW)) if areas.size else float("nan")

    fraction = audit.run(PLAN, analyse, metric="absolute", project=buried_fraction)
    verdict = audit.run(
        PLAN,
        analyse,
        metric="flip",
        project=lambda analysed: buried_fraction(analysed) > MOSTLY,
    )
    exposed = audit.run(
        PLAN,
        analyse,
        metric="relative",
        project=lambda analysed: float(np.nansum(ligand_areas(analysed))),
    )
    fractions = [buried_fraction(run) for run in fraction.runs]
    unmeasured = [int(np.sum(~np.isfinite(ligand_areas(run)))) for run in fraction.runs]
    return {
        "id": code,
        "strata": entry["strata"],
        "ligand_atoms": structure.atom_count - n_protein,
        "universes": len(fraction.runs),
        "ligand_atoms_without_area_by_universe": unmeasured,
        "buried_fraction_min": min(fractions),
        "buried_fraction_max": max(fractions),
        "buried_fraction_range": max(fractions) - min(fractions),
        "mostly_buried_in": sum(f > MOSTLY for f in fractions),
        "conclusion_flips": 0 < sum(f > MOSTLY for f in fractions) < len(fractions),
        "agreement_of_verdict": verdict.agreement_with_first,
        "exposed_area_relative_change": {
            effect.field: effect.mean_change for effect in exposed.effects or []
        }
        if exposed.effects
        else {},
        "shapley_buried_fraction": {share.name: share.share for share in fraction.shapley or []},
        "shapley_exposed_area": {share.name: share.share for share in exposed.shapley or []},
        "balanced": fraction.balanced,
    }


def summarise(rows: list[dict]) -> dict:
    kept = [row for row in rows if "refused" not in row]

    def mean(values) -> float:
        values = list(values)
        return statistics.fmean(values) if values else float("nan")

    ranges = sorted(row["buried_fraction_range"] for row in kept)
    return {
        "complexes_in_corpus": len(rows),
        "refused_by_molframe": len(rows) - len(kept),
        "complexes": len(kept),
        "universes_per_complex": kept[0]["universes"] if kept else 0,
        "conclusion_flips": sum(row["conclusion_flips"] for row in kept) / len(kept),
        "buried_fraction_range_median": statistics.median(ranges),
        "buried_fraction_range_max": ranges[-1],
        "mean_shapley_exposed_area": {
            key: mean(row["shapley_exposed_area"].get(key, 0.0) for row in kept)
            for key in ("hydrogens", "vdw_radii")
        },
        "mean_shapley_buried_fraction": {
            key: mean(row["shapley_buried_fraction"].get(key, 0.0) for row in kept)
            for key in ("hydrogens", "vdw_radii")
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refined-set", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, default=Path(__file__).with_name("corpus.json"))
    parser.add_argument("--limit", type=int, default=0)
    arguments = parser.parse_args()
    corpus = json.loads(arguments.corpus.read_text())
    entries = corpus["entries"][: arguments.limit or None]
    verify_files(arguments.refined_set, [f for e in entries for f in e["files"].values()])
    rows = []
    with tempfile.TemporaryDirectory() as scratch:
        for number, entry in enumerate(entries, start=1):
            try:
                row = run_entry(arguments.refined_set, entry, Path(scratch))
            except molframe.MolframeError as error:
                row = {"id": entry["id"], "refused": error.code, "message": error.message}
            rows.append(row)
            sys.stderr.write(f"{number}/{len(entries)} {entry['id']}\n")
    json.dump(
        {
            "corpus_sha256": sha256(arguments.corpus),
            "molframe_version": molframe.__version__,
            "plan": {
                "universes": PLAN.cost,
                "skipped": PLAN.skipped,
                "balanced": PLAN.balanced,
                "decisions": [
                    {
                        "field": decision.field,
                        "uncertainty": decision.uncertainty,
                        "rationale": decision.rationale,
                        "evidence": decision.evidence,
                    }
                    for decision in PLAN.decisions
                ],
            },
            "buried_below": BURIED_BELOW,
            "mostly": MOSTLY,
            "summary": summarise(rows),
            "entries": rows,
        },
        sys.stdout,
        indent=1,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
