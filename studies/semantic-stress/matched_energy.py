"""Descriptive EvoEF binding-change comparison on the same experimental targets.

Use the author's unmodified executable and libraries outside the checkout.
Training overlap with SKEMPI is unresolved: these scores are not independent
generalization evidence. Contact counts and energies are compared only as ranks.
"""

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from build_corpus import verify_files
from experimental_interfaces import BASELINE, auc, partner_residues

import molframe
from molframe import analysis


def energy(text: str) -> float:
    values = re.findall(r"^Total\s*=\s*([-+\d.eE]+)\s*$", text, re.MULTILINE)
    if len(values) != 1 or not np.isfinite(float(values[0])):
        msg = "reference did not return exactly one finite binding energy"
        raise ValueError(msg)
    return float(values[0])


def residue_identity(path: Path) -> dict:
    import gemmi

    normal = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS"}
    structure = gemmi.read_structure(str(path))
    residues = [
        (
            chain.name,
            residue.seqid.num,
            residue.seqid.icode.strip(),
            normal.get(residue.name, residue.name),
        )
        for chain in structure[0]
        for residue in chain
    ]
    keys = [(chain, number, insertion) for chain, number, insertion, _ in residues]
    if len(keys) != len(set(keys)):
        msg = "reference model has ambiguous residue identities"
        raise ValueError(msg)
    return {key: residue[3] for key, residue in zip(keys, residues, strict=True)}


def transfer(path: Path, entry: dict, destination: Path) -> None:
    import gemmi

    structure = molframe.read(path)
    _, first, second = entry["complex"].split("_")
    _, partners = partner_residues(structure, (set(first), set(second)))
    policy = molframe.AnalysisPolicy(
        hydrogens="exclude", vdw_radii="bondi", contact_def="distance:0.5"
    )
    analysed = analysis.contacts_by_definition(structure, policy=policy)
    output, model = gemmi.Structure(), gemmi.Model("1")
    residues, chains = {}, {}
    for origin in analysed.atom_origin:
        atom = structure.atoms[int(origin)]
        residue = atom.residue
        if residue.index not in partners:
            continue
        if residue.insertion_code:
            msg = "EvoEF transfer refuses insertion codes"
            raise ValueError(msg)
        if residue.index not in residues:
            fresh = gemmi.Residue()
            fresh.name, fresh.het_flag = residue.name, "A"
            fresh.seqid = gemmi.SeqId(residue.auth_number, " ")
            residues[residue.index] = fresh
        fresh_atom = gemmi.Atom()
        fresh_atom.name, fresh_atom.element = atom.name, gemmi.Element(atom.element)
        fresh_atom.pos = gemmi.Position(*atom.coordinate)
        fresh_atom.occ = atom.occupancy if atom.occupancy is not None else 1.0
        fresh_atom.b_iso = atom.b_factor if atom.b_factor is not None else 0.0
        residues[residue.index].add_atom(fresh_atom)
    for residue in structure.residues:
        if residue.index not in residues:
            continue
        chain = residue.chain.auth_label
        if chain not in chains:
            chains[chain] = gemmi.Chain(chain)
        chains[chain].add_residue(residues[residue.index])
    for chain in chains.values():
        model.add_chain(chain)
    output.add_model(model)
    output.write_pdb(str(destination))


def compare(reference: Path, root: Path, entry: dict, scratch: Path) -> dict:
    from scipy.stats import spearmanr

    scratch.mkdir(parents=True, exist_ok=True)
    code, first, second = entry["complex"].split("_")
    transfer(root / "deposited" / f"{code.lower()}.cif", entry, scratch / "input.pdb")
    original = residue_identity(scratch / "input.pdb")

    def execute(label, *arguments: str):
        # The supplied academic reference receives constructed argv without a shell.
        completed = subprocess.run(  # noqa: S603
            [str(reference), *arguments],
            cwd=scratch,
            capture_output=True,
            text=True,
            check=True,
            timeout=300,
        )
        (scratch / f"{label}.log").write_text(completed.stdout + completed.stderr)
        return completed.stdout

    execute("repair", "--command=RepairStructure", "--pdb=input.pdb", "--num_of_runs=3")
    if residue_identity(scratch / "input_Repair.pdb") != original:
        msg = "repair changed residue identity"
        raise ValueError(msg)
    wild = energy(
        execute(
            "reference-energy",
            "--command=ComputeBinding",
            "--pdb=input_Repair.pdb",
            f"--split={first},{second}",
        )
    )
    (scratch / "individual_list.txt").write_text(
        "\n".join(row["mutation"] + ";" for row in entry["residues"]) + "\n"
    )
    execute(
        "mutants",
        "--command=BuildMutant",
        "--pdb=input_Repair.pdb",
        "--mutant_file=individual_list.txt",
        "--num_of_runs=10",
    )
    rows = []
    for number, row in enumerate(entry["residues"], start=1):
        model = f"input_Repair_Model_{number:04d}.pdb"
        expected = dict(original)
        expected[(row["mutation"][1], int(row["mutation"][2:-1]), "")] = "ALA"
        if residue_identity(scratch / model) != expected:
            msg = f"reference mutant changed unintended residue identity: {row['mutation']}"
            raise ValueError(msg)
        mutant = energy(
            execute(
                f"energy-{number:04d}",
                "--command=ComputeBinding",
                f"--pdb={model}",
                f"--split={first},{second}",
            )
        )
        rows.append(
            {
                "mutation": row["mutation"],
                "experimental_ddg": row["ddg_median"],
                "hotspot": row["hotspot"],
                "evoef_ddg": mutant - wild,
                "baseline_contact_degree": row["scores"][BASELINE],
                "median_contact_degree": row["scores"]["policy_median"],
            }
        )
    labels = np.array([row["hotspot"] for row in rows if row["hotspot"] is not None], dtype=bool)
    metrics = {}
    for method in ("evoef_ddg", "baseline_contact_degree", "median_contact_degree"):
        values = np.array([row[method] for row in rows if row["hotspot"] is not None])
        metrics[method] = {
            "auc": auc(labels, values),
            "spearman": float(
                spearmanr(
                    [row[method] for row in rows], [row["experimental_ddg"] for row in rows]
                ).statistic
            ),
        }
    files = [
        {"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for path in sorted(scratch.iterdir())
        if path.is_file()
    ]
    return {
        "complex": entry["complex"],
        "reference_binding_energy": wild,
        "residues": rows,
        "metrics": metrics,
        "files": files,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    reference, root = arguments.reference.resolve(), arguments.results.parent
    source = json.loads(arguments.results.read_text())
    protocol = json.loads((root / "protocol.json").read_text())
    acquisition = json.loads((root / "acquisition.json").read_text())
    if (
        source["protocol_sha256"]
        != hashlib.sha256((root / "protocol.json").read_bytes()).hexdigest()
    ):
        msg = "experimental result belongs to a different protocol"
        raise ValueError(msg)
    verify_files(root, protocol["sources"] + acquisition["files"])
    selected = [
        entry
        for entry in source["entries"]
        if entry["complex"] in {"1BRS_A_D", "1GUA_A_B", "4PWX_AB_CD"}
    ]
    if len(selected) != 3:
        msg = "three specified energy-reference complexes required"
        raise ValueError(msg)
    print(
        json.dumps(
            {
                "reference": "EvoEF 1.1; repair 3 passes; mutants 10 passes; explicit partner split",
                "training_limit": "Parameters optimized using SKEMPI-derived binding changes; pilot training membership unresolved. Descriptive comparison only.",
                "reference_binary_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
                "reference_libraries": [
                    {
                        "path": str(path.relative_to(reference.parent)),
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    }
                    for path in sorted((reference.parent / "library").iterdir())
                    if path.is_file()
                ],
                "results_sha256": hashlib.sha256(arguments.results.read_bytes()).hexdigest(),
                "entries": [
                    compare(reference, root, entry, arguments.output / entry["complex"])
                    for entry in selected
                ],
            },
            indent=1,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
