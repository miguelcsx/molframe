"""Freeze the semantic-stress corpus: which entries, from which files, with which hashes.

Run from the repository root (the dataset path is the only argument that varies):

    uv run python studies/semantic-stress/build_corpus.py \
        --refined-set /path/to/refined-set > studies/semantic-stress/corpus.json

The selection is deterministic: strata are computed from the files themselves,
entries inside a stratum are ordered by a SHA-256 of their identifier, and the
first ``--per-stratum`` are kept. Nothing about the selection depends on a
random seed or on the order of a directory listing.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

BENCH = Path("crates/molframe-bench/data")
METALS = {
    "ZN",
    "MG",
    "CA",
    "MN",
    "FE",
    "FE2",
    "CU",
    "CU1",
    "NI",
    "CO",
    "CD",
    "HG",
    "NA",
    "K",
    "SR",
    "BA",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def index(refined: Path) -> dict[str, dict[str, object]]:
    """Resolution, release year and ligand name for every entry of the dataset's own index."""
    entries: dict[str, dict[str, object]] = {}
    for line in (refined / "index" / "INDEX_refined_data.2020").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        ligand = line.rsplit("(", 1)[-1].rstrip(") ")
        entries[fields[0]] = {
            "resolution": float(fields[1]) if fields[1] != "NMR" else None,
            "year": int(fields[2]),
            "ligand": ligand,
        }
    return entries


def profile(protein: Path) -> dict[str, object]:
    """What the protein file contains, read from the fixed PDB columns."""
    chains: set[str] = set()
    residues: set[tuple[str, str, str]] = set()
    atoms = hydrogens = insertion_residues = hetero_atoms = altlocs = 0
    insertions: set[tuple[str, str, str]] = set()
    for line in protein.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms += 1
        chains.add(line[21])
        key = (line[21], line[22:26], line[26])
        residues.add(key)
        if line[26] != " ":
            insertions.add(key)
        if line.startswith("HETATM"):
            hetero_atoms += 1
        if line[76:78].strip() == "H":
            hydrogens += 1
        if line[16] != " ":
            altlocs += 1
    insertion_residues = len(insertions)
    return {
        "atoms": atoms,
        "residues": len(residues),
        "chains": len(chains),
        "hydrogen_fraction": round(hydrogens / atoms, 4) if atoms else 0.0,
        "insertion_code_residues": insertion_residues,
        "hetero_atoms": hetero_atoms,
        "altloc_atoms": altlocs,
    }


def strata(entry: dict[str, object], facts: dict[str, object]) -> list[str]:
    tags = []
    if facts["insertion_code_residues"]:
        tags.append("insertion-codes")
    if facts["chains"] > 2:
        tags.append("multi-chain")
    if entry["ligand"] in METALS:
        tags.append("metal-ligand")
    resolution = entry["resolution"]
    if resolution is not None:
        tags.append(
            "resolution-high"
            if resolution <= 1.6
            else "resolution-low"
            if resolution >= 2.6
            else "resolution-mid"
        )
    if not tags or tags == ["resolution-mid"]:
        tags.append("control")
    return tags


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refined-set", type=Path, required=True)
    parser.add_argument("--per-stratum", type=int, default=12)
    arguments = parser.parse_args()
    refined = arguments.refined_set
    meta = index(refined)
    buckets: dict[str, list[str]] = {}
    facts_by_id: dict[str, dict[str, object]] = {}
    for code in sorted(meta):
        protein = refined / code / f"{code}_protein.pdb"
        if not protein.exists():
            continue
        facts = profile(protein)
        facts_by_id[code] = facts
        for tag in strata(meta[code], facts):
            buckets.setdefault(tag, []).append(code)
    chosen: dict[str, list[str]] = {}
    for tag, codes in sorted(buckets.items()):
        ordered = sorted(
            codes, key=lambda code: hashlib.sha256(f"{tag}:{code}".encode()).hexdigest()
        )
        for code in ordered[: arguments.per_stratum]:
            chosen.setdefault(code, []).append(tag)

    entries = []
    for code in sorted(chosen):
        base = refined / code
        files = {
            kind: {
                "path": f"{code}/{code}_{kind}.{ext}",
                "sha256": sha256(base / f"{code}_{kind}.{ext}"),
            }
            for kind, ext in (("protein", "pdb"), ("ligand", "sdf"), ("pocket", "pdb"))
        }
        entries.append(
            {
                "id": code,
                "strata": sorted(chosen[code]),
                **meta[code],
                "facts": facts_by_id[code],
                "files": files,
            }
        )

    bench = []
    for code, name in (
        ("1crn", "1crn.cif"),
        ("1ubq", "1ubq.cif"),
        ("4hhb", "4hhb.cif"),
        ("1aon", "1aon.cif.gz"),
        ("2m7c", "2m7c.pdb"),
    ):
        path = BENCH / name
        bench.append(
            {"id": code, "path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
        )

    json.dump(
        {
            "name": "molframe semantic-stress corpus",
            "version": 1,
            "refined_set": {
                "dataset": "PDBbind v2020 refined set (protein files prepared by X-Tool)",
                "root": "<refined-set>",
                "entries": len(meta),
                "strata_sizes": {tag: len(codes) for tag, codes in sorted(buckets.items())},
                "per_stratum": arguments.per_stratum,
                "selection": "sha256(tag:id) ascending within each stratum",
            },
            "entries": entries,
            "bench": bench,
        },
        sys.stdout,
        indent=1,
    )
    print()


if __name__ == "__main__":
    main()
