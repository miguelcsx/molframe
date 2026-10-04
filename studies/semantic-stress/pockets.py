"""How much a ligand pocket depends on defensible decisions.

For every complex of the frozen corpus, the pocket is the set of non-HOH residues in
the receptor file that touch the ligand (including ions and cofactors). It is computed
under the valid combinations of three decisions an analyst
makes without usually noticing:

- whether the hydrogens the file carries count (``hydrogens``: interpretive),
- which van der Waals radii the distance is measured against (``vdw_radii``: algorithmic),
- how much slack is added to the sum of radii (``contact_def``: algorithmic).

The audit reports how often the pocket changes, how much, which decision does it, and
whether the decisions act together. The default pockets of three other tools are then
placed among MolFrame's universes: a default is one point in this space, and the
question is which one.

Run from the repository root:

    PYTHONPATH=python uv run --no-sync \
        --with biopython==1.88 --with biotite==1.7.1 --with gemmi==0.7.5 \
        python studies/semantic-stress/pockets.py \
            --refined-set /path/to/refined-set > studies/semantic-stress/pockets.json
"""

import argparse
import json
import statistics
import sys
import tempfile
import warnings
from importlib import metadata
from pathlib import Path

import numpy as np
from build_corpus import sha256, verify_files

import molframe
from molframe import analysis, audit

warnings.simplefilter("ignore")

TOOL_CUTOFF = 4.0  # the distance at which each tool's neighbour search is asked

SPACE = (
    audit.PolicySpace()
    .vary(
        "hydrogens",
        ["explicit_only", "exclude"],
        rationale="prepared protein files carry added hydrogens; deposited entries do not",
        evidence="PDBbind v2020 protein files are written by X-Tool and are about half hydrogen",
    )
    .vary(
        "vdw_radii",
        ["bondi", "charmm", "alvarez", "amber_united"],
        rationale="four radius sets ship with the library; none is claimed to be correct",
        evidence="molframe-chem radius tables",
    )
    .vary(
        "contact_def",
        ["distance:0.0", "distance:0.25", "distance:0.5"],
        rationale="no slack, a quarter angstrom, and the library's own default of half",
        evidence="ContactDefinition::default",
    )
    .forbid(("hydrogens", "explicit_only"), ("vdw_radii", "amber_united"))
)
PLAN = SPACE.plan(constrained=True)


def parse_sdf(path: Path) -> list[tuple[str, float, float, float]]:
    """Read the atoms of the first molecule of a V2000 SDF file."""
    lines = path.read_text().splitlines()
    count = int(lines[3][0:3])
    return [
        (line[31:34].strip(), float(line[0:10]), float(line[10:20]), float(line[20:30]))
        for line in lines[4 : 4 + count]
    ]


def combine(protein: Path, ligand: Path, destination: Path) -> list[dict[str, str]]:
    """Write the protein and the ligand as one PDB file; return the protein's atom keys."""
    keys = []
    out = []
    for line in protein.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            out.append(line)
            keys.append(
                {
                    "residue": "|".join((line[21], line[22:26].strip(), line[26])),
                    "name": line[17:20].strip(),
                    "element": line[76:78].strip(),
                    "serial": line[6:11].strip(),
                }
            )
    serial = len(out)
    for index, (element, x, y, z) in enumerate(parse_sdf(ligand), start=1):
        out.append(
            f"HETATM{serial + index:5d} {element + str(index):<4s} LIG Z   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
    out.append("END")
    destination.write_text("\n".join(out) + "\n")
    return keys


def jaccard(first: set[str], second: set[str]) -> float:
    union = first | second
    return 1.0 - len(first & second) / len(union) if union else 0.0


def tool_pockets(path: Path, n_protein: int, keys: list[dict[str, str]]) -> dict[str, set[str]]:
    """The pocket each other tool's neighbour search gives at one distance, all atoms."""
    found: dict[str, set[str]] = {}
    polymer = np.array([atom["name"] != "HOH" for atom in keys])

    def residues(indices) -> set[str]:
        return {keys[i]["residue"] for i in indices if i < n_protein and polymer[i]}

    from Bio.PDB import NeighborSearch, PDBParser

    # Biopython groups residues by chain, so its iteration order is not the file's; the
    # PDB serial number ties one of its atoms back to the file.
    index_of_serial = {atom["serial"]: i for i, atom in enumerate(keys)}
    structure = PDBParser(QUIET=True).get_structure("x", str(path))
    atoms = list(structure[0].get_atoms())
    ligand = [atom for atom in atoms if atom.get_parent().get_resname() == "LIG"]
    others = [atom for atom in atoms if atom.get_parent().get_resname() != "LIG"]
    search = NeighborSearch(others)
    hits = {
        index_of_serial[str(near.get_serial_number())]
        for atom in ligand
        for near in search.search(atom.coord, TOOL_CUTOFF, "A")
    }
    found["biopython"] = residues(hits)

    import biotite.structure as bs
    import biotite.structure.io.pdb as pdb

    array = pdb.PDBFile.read(str(path)).get_structure(model=1)
    cells = bs.CellList(array[:n_protein], TOOL_CUTOFF)
    hits = set()
    for coordinate in array.coord[n_protein:]:
        hits.update(int(i) for i in cells.get_atoms(coordinate, TOOL_CUTOFF))
    found["biotite"] = residues(hits)

    import gemmi

    # Gemmi groups atoms by chain name, so its iteration order is not the file's; the PDB
    # serial number is what ties one of its atoms back to the file.
    index_of_serial = {atom["serial"]: i for i, atom in enumerate(keys)}
    model = gemmi.read_structure(str(path))[0]
    neighbours = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate()
    hits = set()
    for chain in model:
        if chain.name != "Z":
            continue
        for residue in chain:
            for atom in residue:
                for mark in neighbours.find_atoms(atom.pos, "\0", radius=TOOL_CUTOFF):
                    near = mark.to_cra(model)
                    if near.chain.name != "Z":
                        hits.add(index_of_serial[str(near.atom.serial)])
    found["gemmi"] = residues(hits)
    return found


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
    if structure.atom_count <= n_protein:
        return {"id": code, "skipped": "no ligand atoms"}
    polymer = np.array([atom["name"] != "HOH" for atom in keys])

    def pocket(analysed) -> list[str]:
        origin = np.asarray(analysed.atom_origin, dtype=np.int64)
        first = origin[np.asarray(analysed.value["first"], dtype=np.int64)]
        second = origin[np.asarray(analysed.value["second"], dtype=np.int64)]
        atoms = first[(first < n_protein) & (second >= n_protein)]
        return sorted({keys[i]["residue"] for i in atoms[polymer[atoms]]})

    result = audit.run(
        PLAN,
        lambda policy: analysis.contacts_by_definition(structure, policy=policy),
        metric="set",
        project=pocket,
    )
    pockets = [set(pocket(run)) for run in result.runs]
    union = set().union(*pockets)
    inter = set.intersection(*pockets)
    effects = {effect.field: effect.mean_change for effect in result.effects or []}
    shapley = {share.name: share.share for share in result.shapley or []}
    classes = {share.name: share.share for share in result.by_class or []}
    tools = {
        name: {
            "size": len(found),
            "equals_a_universe": any(found == universe for universe in pockets),
            "jaccard_to_first": 1.0 - jaccard(found, pockets[0]),
            "nearest_universe_jaccard": 1.0 - min(jaccard(found, universe) for universe in pockets),
        }
        for name, found in tool_pockets(combined, n_protein, keys).items()
    }
    return {
        "id": code,
        "strata": entry["strata"],
        "protein_atoms": n_protein,
        "ligand_atoms": structure.atom_count - n_protein,
        "universes": result.runs.__len__(),
        "balanced": result.balanced,
        "distinct_pockets": len({tuple(sorted(p)) for p in pockets}),
        "union_residues": len(union),
        "invariant_residues": len(inter),
        "sensitive_residues": len(union - inter),
        "agreement_with_first": result.agreement_with_first,
        "mean_distance": result.mean_distance,
        "max_distance": result.max_distance,
        "effects": effects,
        "shapley": shapley,
        "by_class": classes,
        "interactions": {f"{i.first}x{i.second}": i.share for i in result.interactions or []}
        if result.balanced
        else None,
        "higher_order": result.higher_order if result.balanced else None,
        "tools": tools,
    }


def summarise(rows: list[dict]) -> dict:
    refused = [row for row in rows if "refused" in row]
    kept = [row for row in rows if "skipped" not in row and "refused" not in row]

    def mean(values) -> float:
        values = list(values)
        return statistics.fmean(values) if values else float("nan")

    changed = [row for row in kept if row["distinct_pockets"] > 1]
    return {
        "complexes_in_corpus": len(rows),
        "refused_by_molframe": len(refused),
        "refusal_codes": sorted({row["refused"] for row in refused}),
        "complexes": len(kept),
        "pocket_changes_under_some_universe": len(changed) / len(kept),
        "mean_distinct_pockets": mean(row["distinct_pockets"] for row in kept),
        "mean_fraction_of_pocket_residues_sensitive": mean(
            row["sensitive_residues"] / row["union_residues"]
            for row in kept
            if row["union_residues"]
        ),
        "mean_pairwise_jaccard_distance": mean(row["mean_distance"] for row in kept),
        "max_jaccard_distance": max(row["max_distance"] for row in kept),
        "entries_where_hydrogens_alone_change_the_pocket": sum(
            row["effects"]["hydrogens"] > 0 for row in kept
        )
        / len(kept),
        "mean_shapley": {
            key: mean(row["shapley"].get(key, 0.0) for row in kept)
            for key in ("hydrogens", "vdw_radii", "contact_def")
        },
        "mean_by_class": {
            key: mean(row["by_class"].get(key, 0.0) for row in kept)
            for key in ("interpretive", "algorithmic")
        },
        "tools": {
            tool: {
                "default_equals_a_universe": mean(
                    row["tools"][tool]["equals_a_universe"] for row in kept
                ),
                "mean_nearest_universe_jaccard": mean(
                    row["tools"][tool]["nearest_universe_jaccard"] for row in kept
                ),
                "mean_jaccard_to_first_universe": mean(
                    row["tools"][tool]["jaccard_to_first"] for row in kept
                ),
            }
            for tool in ("biopython", "biotite", "gemmi")
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
                # A refusal is a result: the library declined where a lenient tool would not.
                row = {"id": entry["id"], "refused": error.code, "message": error.message}
            rows.append(row)
            sys.stderr.write(f"{number}/{len(entries)} {entry['id']}\n")
    json.dump(
        {
            "corpus_sha256": sha256(arguments.corpus),
            "versions": {
                name: metadata.version(name)
                for name in ("molframe", "biopython", "biotite", "gemmi")
            },
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
            "tool_cutoff": TOOL_CUTOFF,
            "summary": summarise(rows),
            "entries": rows,
        },
        sys.stdout,
        indent=1,
    )
    print()


if __name__ == "__main__":
    main()
