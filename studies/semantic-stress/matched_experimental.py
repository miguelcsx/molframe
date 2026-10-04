"""Check experimental residue joins with Gemmi and interface degree with ProLIF.

The independent reader verifies all measured author identifiers and wild types.
The domain reference verifies twelve matched vdW kernels on three complexes;
shared atoms and radii isolate geometry rather than chemical interpretation.
"""

import argparse
import hashlib
import json
from collections import defaultdict
from importlib import metadata
from itertools import combinations
from pathlib import Path

import numpy as np
from experimental_interfaces import ONE_TO_THREE, POLICIES, partner_residues
from matched_interactions import prolif_transfer

import molframe
from molframe import analysis, chemistry


def chain_sequences(root: Path, entries: list[dict]) -> list[dict]:
    import gemmi

    sequences = []
    for entry in entries:
        code, first, second = entry["complex"].split("_")
        model = gemmi.read_structure(str(root / "deposited" / f"{code.lower()}.cif"))[0]
        for chain in model:
            if chain.name not in set(first + second):
                continue
            sequence = "".join(
                gemmi.find_tabulated_residue(residue.name).one_letter_code.upper()
                for residue in chain
                if gemmi.find_tabulated_residue(residue.name).is_amino_acid()
            )
            sequence = "".join(letter if letter in ONE_TO_THREE else "X" for letter in sequence)
            if sequence:
                sequences.append(
                    {
                        "complex": entry["complex"],
                        "cluster": entry["cluster"],
                        "chain": chain.name,
                        "sequence": sequence,
                    }
                )
    return sequences


def sequence_leakage(root: Path, entries: list[dict]) -> dict:
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    sequences = chain_sequences(root, entries)
    flags, compared = [], 0
    for first, second in combinations(sequences, 2):
        if first["cluster"] == second["cluster"]:
            continue
        alignment = aligner.align(first["sequence"], second["sequence"])[0]
        lengths = [sum(int(end - start) for start, end in blocks) for blocks in alignment.aligned]
        identity = alignment.counts().identities / min(lengths) if min(lengths) else 0
        coverage = [
            length / len(entry["sequence"])
            for length, entry in zip(lengths, (first, second), strict=True)
        ]
        compared += 1
        if identity >= 0.3 and min(coverage) >= 0.8:
            flags.append(
                {
                    "first": [first["complex"], first["chain"]],
                    "second": [second["complex"], second["chain"]],
                    "identity": identity,
                    "coverage": coverage,
                }
            )
    return {
        "design": "Post-run leakage screen: global BLOSUM62, gaps -10/-0.5; flag identity >=30% with >=80% coverage on both chains. Passing is not exhaustive family validation.",
        "chains": len(sequences),
        "cross_cluster_pairs": compared,
        "flagged": flags,
        "sequences": sequences,
    }


def verify_reader(root: Path, entries: list[dict]) -> dict:
    import gemmi

    verified, mismatches = 0, []
    for entry in entries:
        model = gemmi.read_structure(
            str(root / "deposited" / f"{entry['complex'].split('_')[0].lower()}.cif")
        )[0]
        residues = defaultdict(list)
        for chain in model:
            for residue in chain:
                residues[(chain.name, str(residue.seqid))].append(residue.name)
        for row in entry["residues"]:
            mutation = row["mutation"]
            expected = ONE_TO_THREE[mutation[0]]
            found = residues[(mutation[1], mutation[2:-1])]
            if found == [expected]:
                verified += 1
            else:
                mismatches.append(
                    {
                        "complex": entry["complex"],
                        "mutation": mutation,
                        "gemmi": found,
                        "expected": expected,
                    }
                )
    return {"verified_author_residues": verified, "mismatches": mismatches}


def verify_degree(root: Path, entry: dict) -> list[dict]:
    import prolif

    path = root / "deposited" / f"{entry['complex'].split('_')[0].lower()}.cif"
    structure = molframe.read(path)
    _, first, second = entry["complex"].split("_")
    _, partners = partner_residues(structure, (set(first), set(second)))
    rows = []
    for radii, slack in POLICIES:
        policy = molframe.AnalysisPolicy(
            hydrogens="exclude", vdw_radii=radii, contact_def=f"distance:{slack}"
        )
        analysed = analysis.contacts_by_definition(structure, policy=policy)
        indices = list(map(int, analysed.atom_origin))
        atoms = [
            [i for i in indices if partners.get(structure.atoms[i].residue.index) == group]
            for group in (0, 1)
        ]
        values = np.asarray(chemistry.vdw_radii(structure, radii=radii))
        supplied = {structure.atoms[i].element: float(values[i]) for group in atoms for i in group}
        fingerprint = prolif.Fingerprint(
            ["VdWContact"], parameters={"VdWContact": {"tolerance": slack, "vdwradii": supplied}}
        )
        contacts = fingerprint.generate(
            prolif_transfer(structure, atoms[0]),
            prolif_transfer(structure, atoms[1]),
            residues="all",
        )
        neighbours = defaultdict(set)
        for (a, b), bits in contacts.items():
            if bits.any():
                neighbours[a.number - 1].add(b.number - 1)
                neighbours[b.number - 1].add(a.number - 1)
        key = f"{radii}:{slack}"
        mismatches = [
            {
                "mutation": row["mutation"],
                "native": row["scores"][key],
                "prolif": len(neighbours[row["residue"]]),
            }
            for row in entry["residues"]
            if row["scores"][key] != len(neighbours[row["residue"]])
        ]
        rows.append(
            {
                "complex": entry["complex"],
                "policy": key,
                "comparisons": len(entry["residues"]),
                "mismatches": mismatches,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--ids", nargs="+", default=["1BRS_A_D", "1GUA_A_B", "4PWX_AB_CD"])
    arguments = parser.parse_args()
    root = arguments.results.parent
    results = json.loads(arguments.results.read_text())
    protocol = json.loads((root / "protocol.json").read_text())
    acquisition = json.loads((root / "acquisition.json").read_text())
    for source in protocol["sources"] + acquisition["files"]:
        if hashlib.sha256((root / source["path"]).read_bytes()).hexdigest() != source["sha256"]:
            msg = f"frozen input changed: {source['path']}"
            raise ValueError(msg)
    selected = [entry for entry in results["entries"] if entry["complex"] in arguments.ids]
    if {entry["complex"] for entry in selected} != set(arguments.ids):
        msg = "requested reference complex missing from results"
        raise ValueError(msg)
    print(
        json.dumps(
            {
                "results_sha256": hashlib.sha256(arguments.results.read_bytes()).hexdigest(),
                "versions": {
                    name: metadata.version(name)
                    for name in ("gemmi", "prolif", "rdkit", "biopython")
                },
                "reader": verify_reader(root, results["entries"]),
                "sequence_leakage": sequence_leakage(root, results["entries"]),
                "domain_reference": [
                    row for entry in selected for row in verify_degree(root, entry)
                ],
            },
            indent=1,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
